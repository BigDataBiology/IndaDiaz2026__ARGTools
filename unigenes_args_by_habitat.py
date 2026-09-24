#!/usr/bin/env python3
"""
For a set of wanted genes, find which habitats they show up in, based on
per-sample abundance and a sample -> habitat mapping.

Usage:
    python3 genes_by_habitat.py genes_prot_dna.txt GMGC10.sample-abundance.tsv.xz \
        metadata_GMGC10.sample.meta.tsv output.tsv

Output is a long-format tsv: gene<TAB>habitat, one row per (gene, habitat)
pair where the gene had a present (nonzero, if a value column exists)
abundance in at least one sample belonging to that habitat.
"""
import sys
import argparse
import gzip
import bz2
import lzma
import time


def open_maybe_compressed(path, mode='rt'):
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    elif path.endswith('.xz'):
        return lzma.open(path, mode)
    elif path.endswith('.bz2'):
        return bz2.open(path, mode)
    else:
        return open(path, mode)


CHUNK_SIZE = 10_000_000

# Candidate header names to auto-detect the gene/sample/value columns in the
# big abundance file. Override with --gene-col / --sample-col / --value-col
# if none of these match (the script will print the real header and exit
# rather than silently guess wrong).
GENE_COL_CANDIDATES = ['gene', 'gene_id', 'GMGC10_id', 'gmgc10_id', 'unigene', 'orf_id', 'id']
SAMPLE_COL_CANDIDATES = ['sample', 'sample_id', 'Sample', 'SampleID']
VALUE_COL_CANDIDATES = ['abundance', 'value', 'count', 'norm_abundance', 'raw_count', 'abund']


def find_col(header_fields, candidates, forced=None, what=''):
    if forced is not None:
        if forced not in header_fields:
            sys.exit(f'--{what}-col {forced!r} not found in header: {header_fields}')
        return header_fields.index(forced)
    lower = [h.lower() for h in header_fields]
    for cand in candidates:
        if cand.lower() in lower:
            return lower.index(cand.lower())
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('genes_wanted_file', help='one gene id per line')
    ap.add_argument('abundance_file', help='GMGC10.sample-abundance.tsv(.xz/.gz/.bz2)')
    ap.add_argument('metadata_file', help='metadata_GMGC10.sample.meta.tsv')
    ap.add_argument('output_file', help='output tsv: gene<TAB>habitat')
    ap.add_argument('--gene-col', default=None, help='abundance file column name holding the gene id (auto-detected if omitted)')
    ap.add_argument('--sample-col', default=None, help='abundance file column name holding the sample id (auto-detected if omitted)')
    ap.add_argument('--value-col', default=None,
                     help='abundance file column name holding the abundance value (auto-detected if present; '
                          'if the file has no such column, every row is treated as a present/nonzero observation)')
    ap.add_argument('--meta-sample-col', default='sample_id', help='metadata column holding the sample id (default: sample_id)')
    ap.add_argument('--meta-habitat-col', default='habitat', help='metadata column holding the habitat (default: habitat)')
    args = ap.parse_args()

    log_file = args.output_file + '.log'
    log = open(log_file, 'w')

    def log_msg(msg):
        ts = time.strftime('%Y-%m-%d %H:%M:%S')
        line = f'[{ts}] {msg}'
        log.write(line + '\n')
        log.flush()
        print(line, file=sys.stderr)

    # --- wanted genes ---
    with open_maybe_compressed(args.genes_wanted_file) as f:
        wanted = {line.strip() for line in f if line.strip()}
    log_msg(f'Loaded {len(wanted):,} wanted genes from {args.genes_wanted_file}')

    # --- metadata: sample_id -> habitat ---
    sample_to_habitat = {}
    with open_maybe_compressed(args.metadata_file) as f:
        header = f.readline().rstrip('\n').split('\t')
        try:
            sid_idx = header.index(args.meta_sample_col)
            hab_idx = header.index(args.meta_habitat_col)
        except ValueError:
            sys.exit(f'Could not find {args.meta_sample_col!r}/{args.meta_habitat_col!r} '
                      f'in metadata header: {header}')
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            sample_to_habitat[fields[sid_idx]] = fields[hab_idx]
    habitats = sorted(set(sample_to_habitat.values()))
    log_msg(f'Loaded {len(sample_to_habitat):,} samples across {len(habitats)} habitats: {habitats}')

    # --- scan the big abundance file ---
    with open_maybe_compressed(args.abundance_file) as fin:
        header_line = fin.readline().rstrip('\n')
        header_fields = header_line.split('\t')

        gene_idx = find_col(header_fields, GENE_COL_CANDIDATES, args.gene_col, 'gene')
        sample_idx = find_col(header_fields, SAMPLE_COL_CANDIDATES, args.sample_col, 'sample')
        value_idx = find_col(header_fields, VALUE_COL_CANDIDATES, args.value_col, 'value')

        if gene_idx is None or sample_idx is None:
            sys.exit(
                'Could not auto-detect the gene/sample columns in the abundance file.\n'
                f'Header was: {header_fields}\n'
                'Re-run with --gene-col and --sample-col set to the exact column names above.'
            )
        log_msg(
            f'Abundance file columns -> gene: {header_fields[gene_idx]!r} (col {gene_idx}), '
            f'sample: {header_fields[sample_idx]!r} (col {sample_idx}), '
            + (f'value: {header_fields[value_idx]!r} (col {value_idx})'
               if value_idx is not None else
               'value: none detected (every row treated as present)')
        )

        found_pairs = set()          # (gene, habitat)
        genes_seen = set()           # any wanted gene observed at all, any sample
        samples_missing_meta = set()
        line_count = 0
        chunks_done = 0

        for line in fin:
            line_count += 1
            fields = line.rstrip('\n').split('\t')

            gene = fields[gene_idx]
            if gene in wanted:
                if value_idx is not None:
                    try:
                        present = float(fields[value_idx]) > 0
                    except ValueError:
                        present = False
                else:
                    present = True

                if present:
                    sample = fields[sample_idx]
                    genes_seen.add(gene)
                    habitat = sample_to_habitat.get(sample)
                    if habitat is None:
                        samples_missing_meta.add(sample)
                    else:
                        found_pairs.add((gene, habitat))

            if line_count % CHUNK_SIZE == 0:
                chunks_done += 1
                log_msg(f'Chunk {chunks_done} done ({line_count:,} lines read, '
                        f'{len(found_pairs):,} gene/habitat pairs so far)')

    missing_genes = wanted - genes_seen
    log_msg(f'Finished scanning: {line_count:,} lines read, {chunks_done} full chunks')
    log_msg(f'{len(genes_seen):,}/{len(wanted):,} wanted genes observed at all in the abundance file '
            f'({len(missing_genes):,} never seen)')
    if samples_missing_meta:
        log_msg(f'WARNING: {len(samples_missing_meta):,} distinct sample ids in the abundance file '
                 'had no matching entry in the metadata (skipped): '
                 f'{sorted(samples_missing_meta)[:10]}{" ..." if len(samples_missing_meta) > 10 else ""}')

    with open_maybe_compressed(args.output_file, 'wt') as fout:
        fout.write('gene\thabitat\n')
        for gene, habitat in sorted(found_pairs):
            fout.write(f'{gene}\t{habitat}\n')
    log_msg(f'Wrote {len(found_pairs):,} gene/habitat rows to {args.output_file}')

    log.close()


if __name__ == '__main__':
    main()
