import pandas as pd

df = pd.read_csv('faa_args_vs_SARG_filtered_id70_cov60.ranked.tsv', sep='\t')


def count_rankings_per_unigene(pident_th, cov_th):
    filtered = df[(df.pident >= pident_th) & (df.qcovhsp >= cov_th)]
    rank_counts = filtered.groupby('qseqid')['rank'].nunique()
    dist = rank_counts.value_counts().to_dict()
    return dist

thresholds = [(70, 80), (80, 80), (90, 80)]
results = []

for p, c in thresholds:
    dist = count_rankings_per_unigene(p, c)
    results.append({
        'Identity (%)': p,
        'Coverage (%)': c,
        '1 Ranking': dist.get(1, 0),
        '2 Rankings': dist.get(2, 0),
        '3 Rankings': dist.get(3, 0),
        '4 or more': dist.get(4, 0) + dist.get(5, 0)
    })

results_df = pd.DataFrame(results)
print(results_df.to_string(index=False))
