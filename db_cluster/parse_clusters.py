#!/usr/bin/env python3
"""Parses the 4 .clstr files from cluster_cdhit.sh into:
  - one long-format TSV per threshold: protein_id, cluster_id,
    is_representative, pct_identity_to_rep, seq_length
  - one combined wide-format TSV: protein_id -> its cluster_id at each of
    the 4 thresholds, for easy cross-referencing of the same protein
    across thresholds

protein_id is the full merged-fasta header (tag@@@original_header), i.e.
exactly what merge_and_tag.py wrote - so these tables join straight back
onto the merged FASTA.

Usage: python3 parse_clusters.py <cdhit_outdir> [-o OUTDIR]
"""
import argparse
import re
from pathlib import Path

THRESHOLDS = {"90": "clusters_90.clstr", "95": "clusters_95.clstr",
              "975": "clusters_975.clstr", "99": "clusters_99.clstr"}

LINE_RE = re.compile(
    r"^\d+\t\d+aa, >(?P<id>.+)\.\.\. (?:\*|at (?:[+-]/)?(?P<pct>[\d.]+)%)$"
)


def parse_clstr(path):
    """Yields (cluster_id, list of (protein_id, is_representative, pct))."""
    cluster_id = None
    members = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">Cluster"):
                if cluster_id is not None:
                    yield cluster_id, members
                cluster_id = int(line.split()[-1])
                members = []
                continue
            m = LINE_RE.match(line)
            if not m:
                raise ValueError(f"unparsed .clstr line in {path}: {line!r}")
            pid = m.group("id")
            is_rep = line.rstrip().endswith("*")
            pct = float(m.group("pct")) if m.group("pct") else 100.0
            members.append((pid, is_rep, pct))
    if cluster_id is not None:
        yield cluster_id, members


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cdhit_outdir")
    ap.add_argument("-o", "--outdir", default=None,
                     help="defaults to the same directory as cdhit_outdir")
    args = ap.parse_args()

    indir = Path(args.cdhit_outdir)
    outdir = Path(args.outdir) if args.outdir else indir
    outdir.mkdir(parents=True, exist_ok=True)

    # protein_id -> {threshold_tag: cluster_id}
    wide = {}

    for tag, fname in THRESHOLDS.items():
        clstr_path = indir / fname
        if not clstr_path.exists():
            print(f"WARNING: {clstr_path} not found, skipping {tag}%")
            continue

        long_path = outdir / f"cluster_membership_{tag}.tsv"
        n_clusters = 0
        n_members = 0
        long_cols = [
            "cluster_id",
            "protein_id",
            "is_representative",
            "pct_identity_to_rep",
        ]
        with open(long_path, "w") as out:
            out.write("\t".join(long_cols) + "\n")
            for cluster_id, members in parse_clstr(clstr_path):
                n_clusters += 1
                for pid, is_rep, pct in members:
                    row = [cluster_id, pid, is_rep, pct]
                    out.write("\t".join(str(x) for x in row) + "\n")
                    n_members += 1
                    wide.setdefault(pid, {})[tag] = cluster_id

        print(f"{tag}%: {n_clusters} clusters, "
              f"{n_members} proteins -> {long_path}")

    wide_path = outdir / "cluster_membership_all_thresholds.tsv"
    wide_cols = [
        "protein_id",
        "cluster_90",
        "cluster_95",
        "cluster_975",
        "cluster_99",
    ]
    with open(wide_path, "w") as out:
        out.write("\t".join(wide_cols) + "\n")
        for pid in sorted(wide):
            row = wide[pid]
            values = [pid]
            for t in ("90", "95", "975", "99"):
                values.append(str(row.get(t, "")))
            out.write("\t".join(values) + "\n")
    print(f"\nCombined wide table ({len(wide)} proteins) -> {wide_path}")


if __name__ == "__main__":
    main()
