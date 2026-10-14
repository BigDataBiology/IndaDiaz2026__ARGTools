#!/usr/bin/env python
"""
Downloads the gaarangoa/deeparg model+database bundle (DeepARG v2) into
<dir>, and writes STAMP_FILE there once the download is complete.

jugfile.py's prepare_deeparg_db task normally calls this (and every
run_deeparg task waits for it); it can also be run by hand with

    pixi run download-deeparg-db <dir>

This has to be a separate script rather than a function in jugfile.py
because it needs huggingface_hub, which only exists in the 'deeparg' pixi
environment (it comes in as a dependency of deeparg-modern) -- jugfile.py
itself runs under the 'jug' environment. The pixi task is defined under
[feature.deeparg.tasks] so `pixi run download-deeparg-db` always picks that
environment too.

Only the stamp-file name is shared with jugfile.py, which imports it from
here; the huggingface_hub import is deliberately kept inside main() so that
importing this module from the 'jug' environment works.
"""
import argparse
import os

# Written last, so that its presence means "this directory holds a complete
# bundle" -- an interrupted download leaves the files but not the stamp, and
# is retried rather than mistaken for a finished one.
STAMP_FILE = ".download-complete"

REPO_ID = "gaarangoa/deeparg"


def download(dest):
    """Download the bundle into `dest` (created if needed); returns `dest`."""
    from huggingface_hub import snapshot_download

    dest = os.path.abspath(dest)
    os.makedirs(dest, exist_ok=True)
    print(f">>> downloading {REPO_ID} (DeepARG v2 model + database bundle) into {dest}")
    snapshot_download(repo_id=REPO_ID, repo_type="model", local_dir=dest)
    with open(os.path.join(dest, STAMP_FILE), "w") as fh:
        fh.write(f"{REPO_ID}\n")
    return dest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dir", help="destination directory")
    args = parser.parse_args()

    dest = download(args.dir)
    print(f">>> done -- set DEEPARG_HF_DIR = {dest!r} at the top of jugfile.py")


if __name__ == "__main__":
    main()
