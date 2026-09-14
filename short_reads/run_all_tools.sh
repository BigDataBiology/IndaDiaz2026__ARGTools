#!/usr/bin/env bash
# Usage:
#   ./run_all_tools.sh <r1.fastq[.gz]> <r2.fastq[.gz]> <output_dir> [threads]
#
# ONLY_TOOL=fargene|rgi|deeparg runs just one tool instead of all three.
# RGI_LOCALDB_DIR and DEEPARG_HF_DIR must be set (see pixi.toml's
# prepare-rgi-db / download-deeparg-db tasks).
set -euo pipefail

R1="$1"
R2="$2"
OUTDIR="$3"
THREADS="${4:-4}"
ONLY_TOOL="${ONLY_TOOL:-all}"

RGI_LOCALDB_DIR="${RGI_LOCALDB_DIR:-/EDIT/ME/rgi_work_dir}"
DEEPARG_HF_DIR="${DEEPARG_HF_DIR:-/EDIT/ME/deeparg_hf_bundle}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FARGENE_BIN="$SCRIPT_DIR/.pixi/envs/fargene/bin/fargene"
RGI_BIN="$SCRIPT_DIR/.pixi/envs/rgi/bin/rgi"
DEEPARG_BIN="$SCRIPT_DIR/.pixi/envs/deeparg/bin/deeparg"

run_fargene() {
    declare -A models=(
        [class_a]=class_a               [class_b1_b2]=class_b_1_2        [class_b3]=class_b_3
        [class_c]=class_c               [class_d1]=class_d_1             [class_d2]=class_d_2
        [mph]=mph                       [erm_1]=erm_type_a               [erm_2]=erm_type_f
        [tet_enzyme]=tet_enzyme         [tet_rpg]=tet_rpg
        [aph2b]=aminoglycoside_model_g  [aph3p]=aminoglycoside_model_h   [aph6]=aminoglycoside_model_i
        [aac2p]=aminoglycoside_model_a  [aac3_1]=aminoglycoside_model_b  [aac3_2]=aminoglycoside_model_c
        [aac6p_1]=aminoglycoside_model_d [aac6p_2]=aminoglycoside_model_e [aac6p_3]=aminoglycoside_model_f
        [tet_efflux]=tet_efflux         [qnr]=qnr
    )
    mkdir -p "$OUTDIR/fargene"

    # fargene can't read .gz directly -- decompress once into a scratch dir
    # under outdir, reuse across all 22 classes, delete when done.
    local tmp_dir="$OUTDIR/fargene/.decompressed"
    mkdir -p "$tmp_dir"
    local r1_plain="$tmp_dir/r1.fastq"
    local r2_plain="$tmp_dir/r2.fastq"
    gunzip -c "$R1" > "$r1_plain" 2>/dev/null || cp "$R1" "$r1_plain"
    gunzip -c "$R2" > "$r2_plain" 2>/dev/null || cp "$R2" "$r2_plain"

    for class in "${!models[@]}"; do
        model="${models[$class]}"
        echo ">>> fargene: $class ($model)"
        "$FARGENE_BIN" \
            -i "$r1_plain" "$r2_plain" \
            --hmm-model "$model" \
            --meta \
            -o "$OUTDIR/fargene/$class" \
            -p "$THREADS" \
            --force
    done

    rm -rf "$tmp_dir"
}

run_rgi() {
    if [ ! -d "$RGI_LOCALDB_DIR/localDB" ]; then
        echo "error: no localDB/ under RGI_LOCALDB_DIR=$RGI_LOCALDB_DIR -- run" >&2
        echo "  pixi run prepare-rgi-db $RGI_LOCALDB_DIR" >&2
        exit 1
    fi
    mkdir -p "$OUTDIR/rgi"
    echo ">>> rgi bwt"
    ( cd "$RGI_LOCALDB_DIR" && \
      "$RGI_BIN" bwt \
          --read_one "$R1" \
          --read_two "$R2" \
          --output_file "$OUTDIR/rgi/sample.bwt" \
          --local \
          --threads "$THREADS" )
}

run_deeparg() {
    mkdir -p "$OUTDIR/deeparg"

    # deeparg derives intermediate filenames (.paired/.unpaired/.merged/
    # .unmerged) by appending suffixes directly onto whatever path it's
    # given -- symlink inputs into outdir/deeparg/ first so those land there
    # instead of next to the original reads, then delete them when done.
    local r1_link="$OUTDIR/deeparg/$(basename "$R1")"
    local r2_link="$OUTDIR/deeparg/$(basename "$R2")"
    ln -sf "$(cd "$(dirname "$R1")" && pwd)/$(basename "$R1")" "$r1_link"
    ln -sf "$(cd "$(dirname "$R2")" && pwd)/$(basename "$R2")" "$r2_link"

    local hf_flag=()
    if [ -d "$DEEPARG_HF_DIR" ]; then
        hf_flag=(--hf-model-path "$DEEPARG_HF_DIR")
    else
        echo "warning: DEEPARG_HF_DIR=$DEEPARG_HF_DIR doesn't exist yet -- run" >&2
        echo "  pixi run download-deeparg-db $DEEPARG_HF_DIR" >&2
        echo "falling back to a live Hugging Face download instead." >&2
    fi
    echo ">>> deeparg short_reads_pipeline"
    "$DEEPARG_BIN" short_reads_pipeline \
        --forward_pe_file "$r1_link" \
        --reverse_pe_file "$r2_link" \
        --output_file "$OUTDIR/deeparg/sample" \
        "${hf_flag[@]}"

    rm -f "$OUTDIR"/deeparg/*.paired "$OUTDIR"/deeparg/*.unpaired \
          "$OUTDIR"/deeparg/*.merged "$OUTDIR"/deeparg/*.unmerged \
          "$r1_link" "$r2_link"
}

case "$ONLY_TOOL" in
    all)     run_fargene; run_rgi; run_deeparg ;;
    fargene) run_fargene ;;
    rgi)     run_rgi ;;
    deeparg) run_deeparg ;;
    *) echo "ONLY_TOOL must be all|fargene|rgi|deeparg, got: $ONLY_TOOL" >&2; exit 1 ;;
esac

echo ">>> done -- outputs in $OUTDIR/{fargene,rgi,deeparg}"
