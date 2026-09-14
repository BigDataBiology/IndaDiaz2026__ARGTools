#!/usr/bin/env bash
# Usage:
#   ./download_deeparg_db.sh <dir>
#
# <dir> is what you then pass as DEEPARG_HF_DIR to run_all_tools.sh.
set -euo pipefail

DIR="$1"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="$SCRIPT_DIR/.pixi/envs/deeparg/bin/python"

mkdir -p "$DIR"
echo ">>> downloading gaarangoa/deeparg (DeepARG v2 model + database bundle) into $DIR"
"$PYTHON_BIN" -c "
from huggingface_hub import snapshot_download
path = snapshot_download(repo_id='gaarangoa/deeparg', repo_type='model', local_dir='$DIR')
print('done:', path)
"
echo ">>> done -- pass DEEPARG_HF_DIR=$DIR to run_all_tools.sh"
