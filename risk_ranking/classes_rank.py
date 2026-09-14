import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from argnorm.lib import get_aro_mapping_table

MIN_NR_GENES_FOR_FRACTION = 20

mapping_table = get_aro_mapping_table('sarg')
arg_rank = pd.read_csv("ARG_rank.txt", sep="\t", encoding="utf-8-sig", index_col=0)
arg_rank['ARO'] = mapping_table.reindex(arg_rank.index)['ARO']
classes_df = pd.read_csv("ARG_classes.csv")
aro2class = classes_df[['ARO_Term_ID', 'gene_class']].set_index('ARO_Term_ID')
arg_rank = arg_rank.merge(aro2class, left_on='ARO', right_index=True)
cross = pd.crosstab(arg_rank['gene_class'], arg_rank['Rank'])
cross4 = cross[['I', 'II', 'III', 'IV']]
cross4 = cross4[cross4.sum(1) > MIN_NR_GENES_FOR_FRACTION]
frac = pd.DataFrame((cross.values.T / cross.sum(1).values).T,
                    index=cross.index,
                    columns=cross.columns)

fig, axes = plt.subplots(4,1, sharex=True)
for r,ax in zip(['I', 'II', 'III', 'IV'], axes):
    sns.kdeplot(frac[r]*100, ax=ax, clip=(0,100))
    ax.set_yticks([])
    ax.set_title(r)
sns.despine(fig, trim=True)
fig.tight_layout()
ax.set_xlabel("Distribution of ranks within a gene class (%)")
fig.savefig('rank_prob_distribution.png')
