import itertools
import subprocess
from os import makedirs, path
import gzip 
import sys


def fasta_iter(fname, full_header=False):
    header = None
    chunks = []
    if fname.endswith('.gz'):
        import gzip
        op = gzip.open
    else:
        op = open
    with op(fname, 'rt') as f:
        for line in f:
            if line[0] == '>':
                if header is not None:
                    yield header,''.join(chunks)
                line = line[1:].strip()
                if not line:
                    header = ''
                elif full_header:
                    header = line.strip()
                else:
                    header = line.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, ''.join(chunks)


def split_seq_file(fa, MAX_SEQ_CHUNKS):
    makedirs('data', exist_ok=True)
    basename = path.basename(fa)
    cur_n = MAX_SEQ_CHUNKS + 1
    ix = 0
    out = None
    partials = []
    for h, seq in fasta_iter(fa):
        if cur_n >= MAX_SEQ_CHUNKS:
            if fa.endswith('.faa'):
                makedirs('data/faa', exist_ok=True)
                partials.append(f'data/faa/{basename}_block_{ix:04}.faa.gz')
            elif fa.endswith('.fna'):
                makedirs('data/fna', exist_ok=True)
                partials.append(f'data/fna/{basename}_block_{ix:04}.fna.gz')
            else:
                raise ValueError(f'Unexpected file name: {fa}')
            out = gzip.open(partials[-1], compresslevel = 1, mode = 'wt')
            cur_n = 0
            ix += 1
        out.write(f'>{h}\n{seq}\n')
        cur_n +=1
    out.close()
    return partials

data_faa = sys.argv[1] 
data_fna = sys.argv[2] 
MAX_SEQ_CHUNKS = int(sys.argv[3])
#1000000


splits_faa = split_seq_file(data_faa, MAX_SEQ_CHUNKS)
splits_fna = split_seq_file(data_fna, MAX_SEQ_CHUNKS)

