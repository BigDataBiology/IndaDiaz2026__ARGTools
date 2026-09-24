#!/usr/bin/env python3
"""
List unigenes present in samples outside a set of excluded habitats, for two
different habitat-exclusion sets, in a single pass over the abundance file.

Usage:
    python3 count-unigenes_in_interesting_habitats.py \
        GMGC10.sample.meta.tsv.gz GMGC10.sample-abundance.tsv.xz \
        output_not_wanted_habitats.txt.gz output_not_wanted_habitats_2.txt.gz
"""
import argparse
import gzip
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
    ap.add_argument('output_not_wanted_habitats', help="output path for the 'not_wanted_habitats' set (amplicon, built-environment, isolate)")
    ap.add_argument('output_not_wanted_habitats_2', help="output path for the 'not_wanted_habitats_2' set (amplicon, isolate)")
    args = ap.parse_args()

    out_paths = {
        'not_wanted_habitats': args.output_not_wanted_habitats,
        'not_wanted_habitats_2': args.output_not_wanted_habitats_2,
    }

    meta = pd.read_csv(args.meta_file, sep='\t', index_col=0)

    # samples to exclude, per habitat-exclusion set
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

    for name in NOT_WANTED_SETS:
        out_path = out_paths[name]
        with gzip.open(out_path, 'wt') as out:
            out.write(f'{len(keep[name])} unigenes considered\n')
            out.write('\nFull list of unigenes considered:\n')
            for g in sorted(keep[name]):
                out.write(f'{g}\n')
        print(f'Wrote {out_path} ({len(keep[name]):,} unigenes, excluding habitats {NOT_WANTED_SETS[name]})')


if __name__ == '__main__':
    main()
