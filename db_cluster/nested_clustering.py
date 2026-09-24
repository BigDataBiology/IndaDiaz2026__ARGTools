#!/usr/bin/env python3
"""Nested/hierarchical CD-HIT clustering: first splits the full input at
90% identity, then re-clusters each resulting group's members
INDEPENDENTLY (i.e. only against each other, not the whole dataset) at
each subsequent, tighter threshold - 95%, 97.5%, 98%, 99% - producing a
full nested cluster-ID lineage per protein. This is different from
running 4 separate all-vs-all clusterings: a 99% cluster ID here is only
meaningful within its parent 98% group, etc.

Usage: python3 nested_cluster.py <input.faa> <outdir>
Set CDHIT_BIN env var if cd-hit isn't on PATH.
"""
import os
import re
import subprocess
import sys
from pathlib import Path

CDHIT_BIN = os.environ.get("CDHIT_BIN", "cd-hit")
LEVELS = [90, 95, 97.5, 98, 99]

LINE_RE = re.compile(r"^\d+\t\d+aa, >(?P<id>.+)\.\.\. ")

# cd-hit truncates each .clstr header at its FIRST SPACE regardless of -d,
# since FASTA-convention "ID" is everything before the first whitespace -
# our headers routinely have spaces (description text after the ~~~
# fields), so sequences get written to cd-hit with spaces swapped for this
# placeholder, and swapped back when reading cluster membership back out.
SPACE_PLACEHOLDER = "\x01"


def hide_spaces(s):
    return s.replace(" ", SPACE_PLACEHOLDER)


def restore_spaces(s):
    return s.replace(SPACE_PLACEHOLDER, " ")


def read_fasta(path):
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)


def run_cdhit(in_fasta, out_prefix, pct):
    c = pct / 100
    cmd = [CDHIT_BIN, "-i", str(in_fasta), "-o", str(out_prefix),
           "-c", str(c), "-n", "5", "-d", "0", "-M", "0", "-T", "0", "-g", "1"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"cd-hit failed on {in_fasta}:\n{res.stdout}\n{res.stderr}", file=sys.stderr)
        raise SystemExit(1)
    return Path(str(out_prefix) + ".clstr")


def parse_clstr_members(clstr_path):
    """Yields a list of member IDs, one list per cluster, in file order."""
    members = []
    started = False
    for line in open(clstr_path):
        line = line.rstrip("\n")
        if line.startswith(">Cluster"):
            if started:
                yield members
            members = []
            started = True
            continue
        m = LINE_RE.match(line)
        if m:
            members.append(restore_spaces(m.group("id")))
    if started:
        yield members


def clean(pct):
    return str(pct).replace(".", "")


def main():
    in_fasta = Path(sys.argv[1])
    outdir = Path(sys.argv[2])
    outdir.mkdir(parents=True, exist_ok=True)
    tmp_dir = outdir / "tmp_nested"
    tmp_dir.mkdir(exist_ok=True)

    seqs = dict(read_fasta(in_fasta))
    assignment = {pid: [] for pid in seqs}

    to_process = [(list(seqs.keys()), "L")]

    for level_idx, pct in enumerate(LEVELS):
        next_to_process = []
        for group_ids, group_label in to_process:
            if len(group_ids) == 1:
                pid = group_ids[0]
                sub_label = f"{group_label}.0"
                assignment[pid].append(sub_label)
                next_to_process.append(([pid], sub_label))
                continue

            tag = group_label.replace(".", "_")
            sub_fasta = tmp_dir / f"lvl{level_idx}_{tag}.faa"
            with open(sub_fasta, "w") as out:
                for pid in group_ids:
                    out.write(f">{hide_spaces(pid)}\n{seqs[pid]}\n")

            out_prefix = tmp_dir / f"lvl{level_idx}_{tag}_out"
            clstr_path = run_cdhit(sub_fasta, out_prefix, pct)

            for cluster_num, members in enumerate(parse_clstr_members(clstr_path)):
                sub_label = f"{group_label}.{cluster_num}"
                for pid in members:
                    assignment[pid].append(sub_label)
                next_to_process.append((members, sub_label))

        to_process = next_to_process
        print(f"after {pct}% level: {len(to_process)} groups", file=sys.stderr)

    out_tsv = outdir / "nested_cluster_membership.tsv"
    cols = [f"cluster_{clean(pct)}" for pct in LEVELS]
    with open(out_tsv, "w") as out:
        out.write("protein_id\t" + "\t".join(cols) + "\n")
        for pid, labels in assignment.items():
            out.write(pid + "\t" + "\t".join(labels) + "\n")
    print(f"\nWrote {out_tsv}", file=sys.stderr)

    for f in tmp_dir.glob("*"):
        f.unlink()
    tmp_dir.rmdir()


if __name__ == "__main__":
    main()
