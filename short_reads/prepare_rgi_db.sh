#!/usr/bin/env bash
# Usage:
#   ./prepare_rgi_db.sh <work_dir> [card.json]
#
# Downloads CARD v4.0.0 automatically if [card.json] is omitted.
# Set CARD_VERSION to fetch a different version.
# <work_dir> ends up with a localDB/ -- pass it as RGI_LOCALDB_DIR to run_all_tools.sh.
set -euo pipefail

WORK_DIR="$1"
CARD_JSON_ARG="${2:-}"
CARD_VERSION="${CARD_VERSION:-4.0.0}"
CARD_URL="https://card.mcmaster.ca/download/0/broadstreet-v${CARD_VERSION}.tar.bz2"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RGI_BIN="$SCRIPT_DIR/.pixi/envs/rgi/bin/rgi"

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

if [ -n "$CARD_JSON_ARG" ]; then
    echo ">>> using provided card.json: $CARD_JSON_ARG"
    cp "$CARD_JSON_ARG" card.json
else
    echo ">>> downloading CARD v$CARD_VERSION from $CARD_URL"
    curl -fsSL -o card_data.tar.bz2 "$CARD_URL"
    tar xjf card_data.tar.bz2 ./card.json
    rm -f card_data.tar.bz2
fi

echo ">>> rgi card_annotation"
"$RGI_BIN" card_annotation -i card.json > card_annotation.log 2>&1
ver=$(ls card_database_v*.fasta | grep -v _all | sed -E 's/card_database_v(.*)\.fasta/\1/')

echo ">>> rgi load (CARD v$ver)"
"$RGI_BIN" clean --local
"$RGI_BIN" load --card_json card.json \
    --card_annotation "card_database_v${ver}.fasta" \
    --card_annotation_all_models "card_database_v${ver}_all.fasta" \
    --local

echo ">>> materializing localDB/ with a throwaway seed run"
printf ">seed\nATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGT\n" > _seed.fna
"$RGI_BIN" main -a DIAMOND -i _seed.fna -o _seed_out --local --clean -t contig -n 1
rm -f _seed*

echo ">>> done -- localDB/ ready at $WORK_DIR/localDB (CARD v$ver)"
