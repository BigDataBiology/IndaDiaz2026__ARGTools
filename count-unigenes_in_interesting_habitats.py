import pandas as pd
import gzip
meta = pd.read_csv('GMGC10.data/metadata/GMGC10.sample.meta.tsv.gz', sep='\t', index_col=0)
not_wanted_habitats = {'amplicon', 'built-environment', 'isolate'}
not_wanted = set(meta.index[meta['habitat'].map(not_wanted_habitats.__contains__)])

CHUNK_SIZE = 10_000_000

keep = set()
chks = pd.read_csv(
        'GMGC10.data/GMGC10.sample-abundance.tsv.xz',
        chunksize=CHUNK_SIZE, sep='\t', usecols=['Unnamed: 0', 'sample'])

for ch in chks:
    ch = ch[~ch['sample'].map(not_wanted.__contains__)]
    keep.update( ch['Unnamed: 0'] )

with gzip.open('Statement__unigenes-considered.txt', 'wt') as out:
    out.write(f'{len(keep)} unigenes considered\n')
    out.write(f'\nFull list of unigenes considered:\n')
    for g in sorted(keep):
        out.write(f'{g}\n')
