import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from argnorm.lib import get_aro_mapping_table

df = pd.read_csv('faa_args_vs_SARG_filtered_id70_cov60.ranked.tsv', sep='\t')
df_conv = pd.read_csv('../code_R_analysis/output_abundance_diversity_resistome/conversion_ARO_parent_new_level.csv')
mapping_table = get_aro_mapping_table('sarg')

df['aro'] = df['sseqid'].str.extract(r'(ARO:\d+)')
df['mapped_aro'] = df['sseqid'].map(mapping_table['ARO'])
df['aro'] = df['aro'].fillna(df['mapped_aro'])


df = df.merge(df_conv[['Term_ID', 'new_level']], left_on='aro', right_on='Term_ID', how='left')
df['arg_class'] = df['new_level'].fillna('Unmapped')

thresholds = [(70, 60), (80, 70), (90, 70)]

fig1, axes1 = plt.subplots(1, 3, figsize=(18, 5))  # rank distribution
fig2, axes2 = plt.subplots(1, 3, figsize=(22, 8))  # ARG classes rank

for i, (p, c) in enumerate(thresholds):
    filtered = df[(df.pident >= p) & (df.qcovhsp >= c)]
    top = filtered.sort_values('bitscore', ascending=False).drop_duplicates('qseqid')
    
    sns.countplot(data=top, x='rank', hue='rank', 
                  order=['I', 'II', 'III', 'IV', 'notassessed'], 
                  ax=axes1[i], palette='viridis', legend=False)
    axes1[i].set_title(f'Id >= {p}%, Cov >= {c}%\nTotal Unigenes: {len(top)}')
    axes1[i].set_xlabel('Rank')
    axes1[i].set_ylabel('Count')

    mapped_top = top[top['arg_class'] != 'Unmapped']
    top_15_classes = mapped_top['arg_class'].value_counts().head(15).index
    plot_data = mapped_top[mapped_top['arg_class'].isin(top_15_classes)]
    counts = plot_data.groupby(['rank', 'arg_class']).size().unstack(fill_value=0)
    counts = counts.reindex(['I', 'II', 'III', 'IV', 'notassessed'], fill_value=0)
    
  
    counts.plot(kind='bar', stacked=True, ax=axes2[i], colormap='tab20')  
    axes2[i].set_title(f'Top 15 Gene Classes\nId >= {p}%, Cov >= {c}%')
    axes2[i].set_xlabel('Rank')
    axes2[i].set_ylabel('Count')
    axes2[i].tick_params(axis='x', rotation=0) 
    
    if i == 2:
        axes2[i].legend(title='Gene Class', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    else:
        axes2[i].get_legend().remove()

fig1.tight_layout()
fig1.savefig('rank_dist_new_level.png')

fig2.tight_layout()
fig2.savefig('stacked_rank_class_dist_top15.png')
