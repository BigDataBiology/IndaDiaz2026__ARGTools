#!/usr/bin/env python3
import sys
import gzip
import bz2
import lzma


def open_maybe_compressed(path, mode='rt'):
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    elif path.endswith('.xz'):
        return lzma.open(path, mode)
    elif path.endswith('.bz2'):
        return bz2.open(path, mode)
    else:
        return open(path, mode)


def main():
    if len(sys.argv) != 4:
        sys.exit(f'Usage: {sys.argv[0]} genes_wanted.txt big_file output_file')

    wanted_file, big_file, output_file = sys.argv[1:4]

    with open_maybe_compressed(wanted_file) as f:
        wanted = {line.strip() for line in f if line.strip()}

    with open_maybe_compressed(big_file) as fin, open_maybe_compressed(output_file, 'wt') as fout:
        header = fin.readline()
        fout.write(header)
        for line in fin:
            first_col = line.split('\t', 1)[0]
            if first_col in wanted:
                fout.write(line)


if __name__ == '__main__':
    main()