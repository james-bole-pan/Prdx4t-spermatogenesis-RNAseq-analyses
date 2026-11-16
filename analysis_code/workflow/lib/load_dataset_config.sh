#!/usr/bin/env bash

# Shared dataset-config loader. This file is sourced by workflow steps 03-05.
# Dataset configuration is intentionally data-only; processing logic belongs in
# the shared workflow scripts.

load_dataset_config() {
    local selector="$1"
    local script_root="$2"
    local config_path

    if [[ -f "$selector" ]]; then
        config_path="$selector"
    elif [[ -f "${script_root}/datasets/${selector}/dataset.conf" ]]; then
        config_path="${script_root}/datasets/${selector}/dataset.conf"
    else
        printf 'ERROR: unknown dataset %s; expected a config file or one of:\n' "$selector" >&2
        find "${script_root}/datasets" -mindepth 2 -maxdepth 2 -name dataset.conf \
            -exec dirname {} \; 2>/dev/null | sed 's|.*/|  |' >&2
        return 1
    fi

    config_path="$(cd -- "$(dirname -- "$config_path")" && pwd)/$(basename -- "$config_path")"
    # shellcheck disable=SC1090
    source "$config_path"

    local required
    for required in DATASET_ID DATASET_LABEL DEFAULT_DATASET_ROOT DEFAULT_REFERENCE_DIR RUN_MANIFEST_FILE PLOT_DATASET_LABEL MAIN_STAGE_GROUPS; do
        [[ -n "${!required:-}" ]] || {
            printf 'ERROR: %s is not defined in %s\n' "$required" "$config_path" >&2
            return 1
        }
    done

    DATASET_CONFIG="$config_path"
    DATASET_CONFIG_DIR="$(dirname -- "$config_path")"
    RUN_MANIFEST="${DATASET_CONFIG_DIR}/${RUN_MANIFEST_FILE}"
    [[ -r "$RUN_MANIFEST" ]] || {
        printf 'ERROR: run manifest is not readable: %s\n' "$RUN_MANIFEST" >&2
        return 1
    }
}

dataset_selector_from_args() {
    local selector="${PRDX4T_DATASET:-PRJNA548952}"
    local expect_value=0
    local argument
    for argument in "$@"; do
        if ((expect_value)); then
            selector="$argument"
            expect_value=0
        elif [[ "$argument" == "--dataset" ]]; then
            expect_value=1
        fi
    done
    ((expect_value == 0)) || {
        printf 'ERROR: --dataset requires an ID or config file\n' >&2
        return 1
    }
    printf '%s\n' "$selector"
}
