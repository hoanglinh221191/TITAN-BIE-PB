#!/usr/bin/env bash
set -euo pipefail

# GROMACS slab-density workflow for the phosphorylated/WT N-protein-RNA runs.
# Override defaults with environment variables, for example:
#   ROOT=/path/to/slab START_PS=1000000 STRIDE_PS=2000 bash gmx_slab_density_temperature.sh

GMX=${GMX:-gmx_mpi}
ROOT=${ROOT:-"$PWD"}
START_PS=${START_PS:-1000000}
END_PS=${END_PS:-auto}
STRIDE_PS=${STRIDE_PS:-2000}
SLICES=${SLICES:-300}
BLOCKS=${BLOCKS:-5}
TEMPERATURES=${TEMPERATURES:-"290 300 310 320 330 340 350 360 370"}
OUTPUT_NAME=${OUTPUT_NAME:-gmx_density_1000ns}

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

command -v "$GMX" >/dev/null 2>&1 || die "$GMX was not found in PATH"
[[ "$BLOCKS" =~ ^[1-9][0-9]*$ ]] || die "BLOCKS must be a positive integer"
[[ "$SLICES" =~ ^[1-9][0-9]*$ ]] || die "SLICES must be a positive integer"

get_end_time_ps() {
    local trajectory=$1
    local check_log=$2
    local check_output
    local detected

    check_output=$("$GMX" check -f "$trajectory" 2>&1 || true)
    printf '%s\n' "$check_output" >"$check_log"
    detected=$(printf '%s\n' "$check_output" | awk '
        /Last frame/ {
            for (i = 1; i <= NF; ++i) {
                token = $i
                gsub(/[:,]/, "", token)
                if (token == "time" && i < NF) {
                    last_frame_time = $(i + 1)
                    break
                }
            }
        }
        /Reading frame/ {
            for (i = 1; i <= NF; ++i) {
                token = $i
                gsub(/[:,]/, "", token)
                if (token == "time" && i < NF) {
                    reading_time = $(i + 1)
                    if (!have_first_time) {
                        first_time = reading_time
                        have_first_time = 1
                    }
                    break
                }
            }
        }
        /^[[:space:]]*Time[[:space:]]+[0-9]+[[:space:]]+/ {
            frame_count = $2
            frame_timestep = $3
        }
        END {
            if (last_frame_time != "") {
                value = last_frame_time
            } else if (have_first_time && frame_count > 0 && frame_timestep >= 0) {
                value = first_time + (frame_count - 1) * frame_timestep
            } else {
                value = reading_time
            }
            gsub(/[,;()]/, "", value)
            if (value ~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) {
                printf "%.10g\n", value
            }
        }
    ')
    if [[ -z "$detected" ]]; then
        printf 'Could not parse the last frame time. Raw gmx check output:\n%s\n' \
            "$check_output" >&2
        return 1
    fi
    printf '%s\n' "$detected"
}

require_index_group() {
    local index_file=$1
    local group_name=$2
    grep -Eq "^[[:space:]]*\[[[:space:]]*${group_name}[[:space:]]*\]" "$index_file" \
        || die "Group '${group_name}' is missing from ${index_file}"
}

run_density() {
    local tpr=$1
    local trajectory=$2
    local index_file=$3
    local group_name=$4
    local begin_ps=$5
    local end_ps=$6
    local output=$7

    printf '%s\n' "$group_name" \
        | "$GMX" density \
            -s "$tpr" \
            -f "$trajectory" \
            -n "$index_file" \
            -b "$begin_ps" \
            -e "$end_ps" \
            -d Z \
            -sl "$SLICES" \
            -dens mass \
            -o "$output"
}

for temperature in $TEMPERATURES; do
    source_dir="${ROOT}/${temperature}"
    trajectory="${source_dir}/md1_1.xtc"
    full_tpr="${source_dir}/md1_1.tpr"
    full_index="${source_dir}/index.ndx"

    if [[ ! -s "$trajectory" || ! -s "$full_tpr" || ! -s "$full_index" ]]; then
        log "Skipping ${temperature} K: md1_1.xtc, md1_1.tpr, or index.ndx is missing"
        continue
    fi

    output_dir="${source_dir}/${OUTPUT_NAME}"
    mkdir -p "$output_dir/blocks"
    log "Starting ${temperature} K -> ${output_dir}"

    subset_tpr="${output_dir}/protein_rna.tpr"
    subset_index="${output_dir}/protein_rna.ndx"
    centered_xtc="${output_dir}/protein_rna_centered.xtc"

    # Build a topology whose atom order exactly matches the selected trajectory.
    printf 'Protein_RNA\n' \
        | "$GMX" convert-tpr \
            -s "$full_tpr" \
            -n "$full_index" \
            -o "$subset_tpr"

    # Make every molecule whole, center the common Protein_RNA assembly, and
    # output the same Protein_RNA selection. Do not center Protein/RNA separately.
    printf 'Protein_RNA\nProtein_RNA\n' \
        | "$GMX" trjconv \
            -s "$full_tpr" \
            -f "$trajectory" \
            -n "$full_index" \
            -b "$START_PS" \
            -dt "$STRIDE_PS" \
            -pbc mol \
            -ur compact \
            -center \
            -boxcenter rect \
            -o "$centered_xtc"

    # Regenerate indices against the subset TPR. In this topology, System is
    # exactly Protein_RNA, while Protein and RNA remain separate selections.
    printf 'q\n' | "$GMX" make_ndx -f "$subset_tpr" -o "$subset_index"
    require_index_group "$subset_index" Protein
    require_index_group "$subset_index" RNA
    require_index_group "$subset_index" System

    if [[ "$END_PS" == "auto" ]]; then
        end_ps=$(get_end_time_ps "$centered_xtc" "${output_dir}/gmx_check_centered.log") \
            || die "Could not determine the final time of ${centered_xtc}"
    else
        end_ps=$END_PS
    fi

    awk -v start="$START_PS" -v end="$end_ps" 'BEGIN {exit !(end > start)}' \
        || die "End time ${end_ps} ps is not later than START_PS=${START_PS} ps"
    printf '%s\n' "$end_ps" >"${output_dir}/trajectory_end_ps.txt"

    run_density "$subset_tpr" "$centered_xtc" "$subset_index" Protein \
        "$START_PS" "$end_ps" "${output_dir}/density_mass_protein.xvg"
    run_density "$subset_tpr" "$centered_xtc" "$subset_index" RNA \
        "$START_PS" "$end_ps" "${output_dir}/density_mass_rna.xvg"
    run_density "$subset_tpr" "$centered_xtc" "$subset_index" System \
        "$START_PS" "$end_ps" "${output_dir}/density_mass_protein_rna.xvg"

    # Independent time blocks for uncertainty estimation. Boundaries may share
    # one frame; this has negligible impact for long blocks and preserves all data.
    for ((block=1; block<=BLOCKS; ++block)); do
        block_begin=$(awk -v s="$START_PS" -v e="$end_ps" -v i="$((block-1))" -v n="$BLOCKS" \
            'BEGIN {printf "%.3f", s + (e-s)*i/n}')
        block_end=$(awk -v s="$START_PS" -v e="$end_ps" -v i="$block" -v n="$BLOCKS" \
            'BEGIN {printf "%.3f", s + (e-s)*i/n}')
        block_tag=$(printf 'block%02d' "$block")

        run_density "$subset_tpr" "$centered_xtc" "$subset_index" Protein \
            "$block_begin" "$block_end" \
            "${output_dir}/blocks/${block_tag}_density_mass_protein.xvg"
        run_density "$subset_tpr" "$centered_xtc" "$subset_index" RNA \
            "$block_begin" "$block_end" \
            "${output_dir}/blocks/${block_tag}_density_mass_rna.xvg"
        run_density "$subset_tpr" "$centered_xtc" "$subset_index" System \
            "$block_begin" "$block_end" \
            "${output_dir}/blocks/${block_tag}_density_mass_protein_rna.xvg"
    done

    cat >"${output_dir}/density_run_metadata.txt" <<EOF
temperature_K=${temperature}
start_ps=${START_PS}
end_ps=${end_ps}
stride_ps=${STRIDE_PS}
slices=${SLICES}
blocks=${BLOCKS}
centering_group=Protein_RNA
density_type=mass
density_axis=Z
EOF
    log "Finished ${temperature} K"
done

log "All available temperatures completed"
