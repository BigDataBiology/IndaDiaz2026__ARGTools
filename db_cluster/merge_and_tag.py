#!/usr/bin/env python3
"""Merges 9 ARG reference sources into one amino-acid FASTA.

Six sources are nucleotide and get translated by checking all 6 reading
frames (3 forward, 3 reverse-complement) and keeping the single longest
open reading frame found anywhere. Two sources are already protein and are
copied through unchanged. The RGI/CARD source is pulled directly out of
card.json's protein_sequence fields (protein homolog model + protein
overexpression model entries - the two model types RGI-DIAMOND actually
reports hits against) rather than from card_database_v4.0.0_all.fasta -
that file turned out to be nucleotide, not protein, and RGI's own DIAMOND
search is itself built from a protein FASTA derived from card.json this
same way, so this is also the more authoritative source, not just a
workaround.

Every header becomes: >{tag}@@@{original header, unmodified}

Usage: python3 merge_and_tag.py [-o OUT_FASTA]
Edit SOURCES / CARD_JSON_PATH below if any path changes.
"""
import argparse
import json
import sys
from pathlib import Path

# (path, tag, needs_translation)
SOURCES = [
    ("/work/microbiome/users/juan/resfinder_databases/resfinder_db/all.fsa", "resfinder", True),
    ("/work/microbiome/users/juan/e/abricate2026/db/resfinder/sequences", "abricate-resfinder", True),
    ("/work/microbiome/users/juan/e/abricate2026/db/argannot/sequences", "abricate-argannot", True),
    ("/work/microbiome/users/juan/e/abricate2026/db/megares/sequences", "abricate-megares", True),
    ("/work/microbiome/users/juan/e/abricate2026/db/card/sequences", "abricate-card", True),
    ("/work/microbiome/users/juan/e/abricate2026/db/ncbi/sequences", "abricate-ncbi", True),
    ("/work/microbiome/users/juan/deeparg_data/database/v2/features.fasta", "deeparg", False),
    ("/work/microbiome/users/juan/e/amrfinder/share/amrfinderplus/data/2024-12-18.1/AMRProt.fa", "amrfinderplus", False),
]

CARD_JSON_PATH = "/work/microbiome/users/juan/rgi/card.json"
CARD_JSON_TAG = "rgi-card"


def extract_card_json_proteins(card_json_path):
    """Yields (header, protein_sequence) for every protein-homolog-model
    entry in card.json, with the header built as ARO:X|ID:Y|Name:Z|NCBI:W -
    the same convention card_database_v4.0.0_all.fasta used, so every piece
    of downstream code that parses that header format (rgi_model_id
    extraction etc.) keeps working unchanged."""
    INCLUDED_MODEL_TYPES = {"protein homolog model", "protein overexpression model"}
    card = json.load(open(card_json_path))
    for key, model in card.items():
        if not isinstance(model, dict):
            continue
        if model.get("model_type") not in INCLUDED_MODEL_TYPES:
            continue
        aro_accession = model.get("ARO_accession", "")
        model_name = (model.get("model_name") or "").replace(" ", "_")
        seqs = model.get("model_sequences", {}).get("sequence", {})
        for seqid, info in seqs.items():
            prot = info.get("protein_sequence")
            if not prot or not prot.get("sequence"):
                continue
            dna = info.get("dna_sequence", {})
            dna_accession = dna.get("accession", "")
            header = f"ARO:{aro_accession}|ID:{key}|Name:{model_name}|NCBI:{dna_accession}"
            yield header, prot["sequence"]

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(nt_seq):
    return nt_seq.translate(COMPLEMENT)[::-1]


def translate_frame(nt_seq):
    return "".join(CODON_TABLE.get(nt_seq[i:i+3], "X") for i in range(0, len(nt_seq) - 2, 3))


def longest_orf(nt_seq):
    """All 6 reading frames, longest stop-delimited segment wins. No ATG
    requirement (bacterial genes often start GTG/TTG; some reference
    entries carry a little flanking sequence)."""
    nt_seq = nt_seq.upper().replace("U", "T")
    rc = reverse_complement(nt_seq)
    best = ""
    for strand_seq in (nt_seq, rc):
        for frame in range(3):
            for segment in translate_frame(strand_seq[frame:]).split("*"):
                if len(segment) > len(best):
                    best = segment
    return best


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default="merged_reference_proteins.faa")
    args = ap.parse_args()

    total = 0
    short_orf_rows = []
    with open(args.out, "w") as out:
        for path, tag, needs_translation in SOURCES:
            p = Path(path)
            if not p.exists():
                print(f"ERROR: not found: {path}", file=sys.stderr)
                sys.exit(1)
            n = 0
            for header, seq in read_fasta(p):
                if not seq:
                    continue
                if needs_translation:
                    if len(seq) < 3:
                        continue
                    aa = longest_orf(seq)
                    if not aa:
                        continue
                    if len(aa) * 3 < len(seq) * 0.9:
                        short_orf_rows.append((tag, header, len(seq), len(aa), round(len(aa) * 3 / len(seq), 3)))
                else:
                    aa = seq
                out.write(f">{tag}@@@{header}\n{aa}\n")
                n += 1
                total += 1
            kind = "translated (best-of-6-frames)" if needs_translation else "protein, copied as-is"
            print(f"{tag}: {n} sequences [{kind}]", file=sys.stderr)

        p = Path(CARD_JSON_PATH)
        if not p.exists():
            print(f"ERROR: not found: {CARD_JSON_PATH}", file=sys.stderr)
            sys.exit(1)
        n = 0
        for header, aa in extract_card_json_proteins(p):
            out.write(f">{CARD_JSON_TAG}@@@{header}\n{aa}\n")
            n += 1
            total += 1
        print(f"{CARD_JSON_TAG}: {n} sequences [protein, straight from card.json]", file=sys.stderr)

    print(f"\nTotal merged sequences: {total}", file=sys.stderr)
    print(f"Output: {args.out}", file=sys.stderr)
    if short_orf_rows:
        report_path = Path(args.out).with_suffix("").as_posix() + ".short_orf_report.tsv"
        with open(report_path, "w") as rep:
            rep.write("source\theader\tnt_length\tbest_orf_aa_length\tcoverage_fraction\n")
            for row in short_orf_rows:
                rep.write("\t".join(str(x) for x in row) + "\n")
        print(f"NOTE: {len(short_orf_rows)} translated entries came up short "
              f"(<90% of expected length even in the best of 6 frames) - "
              f"likely fragmented/partial/non-coding entries.\n"
              f"      Full list written to: {report_path}", file=sys.stderr)
        by_source = {}
        for tag, *_ in short_orf_rows:
            by_source[tag] = by_source.get(tag, 0) + 1
        for tag, n in sorted(by_source.items(), key=lambda kv: -kv[1]):
            print(f"        {tag}: {n}", file=sys.stderr)


if __name__ == "__main__":
    main()
