#!/usr/bin/env python3
"""
Count unigenes present in samples outside a set of excluded habitats, for two
different habitat-exclusion sets, in a single pass over the abundance file.

Usage:
    python3 count-unigenes_in_interesting_habitats.py \
        GMGC10.sample.meta.tsv.gz GMGC10.sample-abundance.tsv.xz counts.txt
"""
import argparse
import pandas as pd

CHUNK_SIZE = 10_000_000

NOT_WANTED_SETS = {
    'not_wanted_habitats': {'amplicon', 'built-environment', 'isolate'},
    'not_wanted_habitats_2': {'amplicon', 'isolate'},
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('meta_file', help='GMGC10.sample.meta.tsv.gz')
    ap.add_argument('abundance_file', help='GMGC10.sample-abundance.tsv.xz')
    ap.add_argument('counts_out', help='output path for the counts summary')
    args = ap.parse_args()

    meta = pd.read_csv(args.meta_file, sep='\t', index_col=0)

    exclude_samples = {
        name: set(meta.index[meta['habitat'].map(habitats.__contains__)])
        for name, habitats in NOT_WANTED_SETS.items()
    }

    keep = {name: set() for name in NOT_WANTED_SETS}

    chks = pd.read_csv(
        args.abundance_file,
        chunksize=CHUNK_SIZE, sep='\t', usecols=['Unnamed: 0', 'sample'])

    for i, ch in enumerate(chks, 1):
        for name, not_wanted in exclude_samples.items():
            sub = ch[~ch['sample'].map(not_wanted.__contains__)]
            keep[name].update(sub['Unnamed: 0'])
        print(f'chunk {i} done')

    with open(args.counts_out, 'w') as out:
        for name in NOT_WANTED_SETS:
            line = f'{name} (excluding habitats {NOT_WANTED_SETS[name]}): {len(keep[name]):,} unigenes'
            print(line)
            out.write(line + '\n')

    print(f'Wrote counts to {args.counts_out}')


if __name__ == '__main__':
    main()
