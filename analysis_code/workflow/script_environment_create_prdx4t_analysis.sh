#!/usr/bin/env bash
set -Eeuo pipefail

readonly DEFAULT_ENV_NAME="prdx4t-analysis"
readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

env_name="${PRDX4T_ENV_NAME:-$DEFAULT_ENV_NAME}"
platform=""
dry_run=0

usage() {
    cat <<EOF
Usage: $(basename "$0") [--name ENV_NAME] [--platform PLATFORM] [--dry-run]

Create or update a Mamba environment containing the software required to
run the PRDX4t RNA-seq analysis pipeline.

Options:
  --name NAME          Environment name (default: ${DEFAULT_ENV_NAME})
  --platform PLATFORM  Conda platform (auto-detected by default)
  --dry-run            Resolve packages without creating/updating the environment
  -h, --help           Show this help

Environment:
  PRDX4T_ENV_NAME      Alternative environment name
  MAMBA_EXE            Explicit path to the mamba executable

On Apple silicon, osx-64 is selected intentionally because the required
Bioconda command-line tools are distributed for Intel macOS and run through
Rosetta 2.
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

while (($#)); do
    case "$1" in
        --name)
            (($# >= 2)) || die "--name requires a value"
            env_name="$2"
            shift 2
            ;;
        --platform)
            (($# >= 2)) || die "--platform requires a value"
            platform="$2"
            shift 2
            ;;
        --dry-run)
            dry_run=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "unknown argument: $1 (use --help)"
            ;;
    esac
done

[[ "$env_name" =~ ^[A-Za-z0-9._-]+$ ]] || die "environment name contains unsupported characters: $env_name"

mamba_exe="${MAMBA_EXE:-}"
if [[ -n "$mamba_exe" ]]; then
    :
elif command -v mamba >/dev/null 2>&1; then
    mamba_exe="$(command -v mamba)"
elif [[ -x "${HOME}/miniforge3/bin/mamba" ]]; then
    mamba_exe="${HOME}/miniforge3/bin/mamba"
else
    die "mamba was not found. Install Miniforge or set MAMBA_EXE."
fi
[[ -x "$mamba_exe" ]] || die "mamba is not executable: $mamba_exe"

if [[ -z "$platform" ]]; then
    case "$(uname -s)-$(uname -m)" in
        Darwin-arm64)
            platform="osx-64"
            if ! /usr/bin/arch -x86_64 /usr/bin/true >/dev/null 2>&1; then
                die "Rosetta 2 is required. Install it with: softwareupdate --install-rosetta"
            fi
            ;;
        Darwin-x86_64)
            platform="osx-64"
            ;;
        Linux-x86_64)
            platform="linux-64"
            ;;
        *)
            die "unsupported host; specify a Conda platform with --platform"
            ;;
    esac
fi

# Command-line tool versions are pinned for reproducibility. Python libraries
# use the versions used for the figure.
packages=(
    "python=3.11"
    "star=2.7.11b"
    "sra-tools=3.1.1"
    "bowtie2=2.5.4"
    "samtools=1.20"
    "fastp=0.23.4"
    "pigz"
    "pysam=0.22.1"
    "numpy=1.26.4"
    "matplotlib=3.8.4"
)

if "$mamba_exe" env list | awk -v wanted="$env_name" '$1 == wanted { found=1 } END { exit !found }'; then
    action="install"
    printf 'Updating existing environment: %s\n' "$env_name"
else
    action="create"
    printf 'Creating environment: %s\n' "$env_name"
fi
printf 'Conda platform: %s\n' "$platform"
printf 'Mamba: %s\n' "$mamba_exe"

mamba_args=(
    "$action"
    --yes
    --name "$env_name"
    --platform "$platform"
    --override-channels
    --channel conda-forge
    --channel bioconda
    --strict-channel-priority
)
if ((dry_run)); then
    mamba_args+=(--dry-run)
fi

"$mamba_exe" "${mamba_args[@]}" "${packages[@]}"

if ((dry_run)); then
    printf '\nDry run completed; no environment changes were made.\n'
    exit 0
fi

printf '\nVerifying executables and Python imports...\n'
"$mamba_exe" run --name "$env_name" bash -c '
set -Eeuo pipefail
printf "STAR:         %s\n" "$(STAR --version)"
printf "Bowtie2:      %s\n" "$(bowtie2 --version | sed -n "1p")"
printf "samtools:     %s\n" "$(samtools --version | sed -n "1p")"
printf "fastp:        %s\n" "$(fastp --version 2>&1 | sed -n "1p")"
printf "SRA Toolkit:  %s\n" "$(fasterq-dump --version 2>&1 | awk "NF {print; exit}")"
python -c "import pysam, numpy, matplotlib; print(\"Python stack: pysam\", pysam.__version__, \"numpy\", numpy.__version__, \"matplotlib\", matplotlib.__version__)"
'


cat <<EOF

Environment ready. Activate it with:

  mamba activate ${env_name}

If shell activation is unavailable, run commands with:

  ${mamba_exe} run --name ${env_name} COMMAND
EOF
