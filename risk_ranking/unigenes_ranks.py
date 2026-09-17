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

thresholds = [(70, 80), (80, 80), (90, 80), (95, 80)]
rank_order = ['I', 'II', 'III', 'IV', 'notassessed']

global_mapped = df[df['arg_class'] != 'Unmapped']
global_top_15 = global_mapped['arg_class'].value_counts().head(15).index

# Rank Distribution and ARG Classes Rank 
fig1, axes1 = plt.subplots(2, 2, figsize=(22, 14))  # rank distribution
axes1 = axes1.flatten()  
fig2, axes2 = plt.subplots(2, 2, figsize=(22, 14))  # ARG classes rank
axes2 = axes2.flatten()  

for i, (p, c) in enumerate(thresholds):
    filtered = df[(df.pident >= p) & (df.qcovhsp >= c)]
    #top = filtered.sort_values('bitscore', ascending=False).drop_duplicates('qseqid')
    
    sns.countplot(data=filtered, x='rank', hue='rank', 
                  order=rank_order, hue_order=rank_order,
                  ax=axes1[i], palette='viridis', legend=False)
    axes1[i].set_title(f'Id >= {p}%, Cov >= {c}%\nTotal Unigenes: {len(filtered)}', fontsize=14, fontweight='bold')
    axes1[i].set_xlabel('Rank', fontsize=12, fontweight='bold')
    axes1[i].set_ylabel('Count', fontsize=12, fontweight='bold')

    mapped_top = filtered[filtered['arg_class'] != 'Unmapped']
    plot_data = mapped_top[mapped_top['arg_class'].isin(global_top_15)]
    counts = (plot_data.groupby(['rank', 'arg_class']).size().unstack(fill_value=0))
    counts = counts.reindex(index=rank_order, columns=global_top_15, fill_value=0)
    
    counts.plot(kind='bar', stacked=True, ax=axes2[i], colormap='tab20')  
    axes2[i].set_title(f'Top 15 Gene Classes\nId >= {p}%, Cov >= {c}%', fontsize=14, fontweight='bold')
    axes2[i].set_xlabel('Rank', fontsize=12, fontweight='bold')
    axes2[i].set_ylabel('Count', fontsize=12, fontweight='bold')
    axes2[i].tick_params(axis='x', rotation=0) 
    
    if i == 3:
        axes2[i].legend(title='Gene Class', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    else:
        legend = axes2[i].get_legend()
        if legend is not None:
            legend.remove()
        

fig1.tight_layout()
fig1.savefig('rank_dist_new_level.png', dpi=300, bbox_inches='tight')

fig2.tight_layout()
fig2.savefig('stacked_rank_class_dist_top15.png', dpi=300, bbox_inches='tight')

# Rank Overlap Matrix
fig3, axes3 = plt.subplots(2, 2, figsize=(22, 14))
axes3 = axes3.flatten()  

for i, (p, c) in enumerate(thresholds):
    filtered = df[(df.pident >= p) & (df.qcovhsp >= c)]
    
    if not filtered.empty:
        rank_dummies = pd.get_dummies(filtered['rank'])
        unigene_ranks = rank_dummies.groupby(filtered['qseqid']).max().astype(int)
        
        overlap_matrix = unigene_ranks.T.dot(unigene_ranks)
        
        overlap_matrix = overlap_matrix.reindex(index=rank_order, columns=rank_order, fill_value=0)
        
        sns.heatmap(overlap_matrix, annot=True, fmt='d', cmap='Purples', ax=axes3[i], cbar=(i==3), annot_kws={"size": 12, "weight": "bold"})
        axes3[i].set_title(f'All-Hits Rank Overlap\nId >= {p}%, Cov >= {c}%', fontsize=14, fontweight='bold')
        axes3[i].set_ylabel('Rank A', fontsize=12, fontweight='bold')
        axes3[i].set_xlabel('Rank B', fontsize=12, fontweight='bold')

fig3.tight_layout()
fig3.savefig('rank_overlap_matrix.png', dpi=300, bbox_inches='tight')


# Normalized Bitscore Diffs and Riskiest Hit Shift
risk_map = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'notassessed': 5}
df['risk_score'] = df['rank'].map(risk_map)

high_thresholds = [(90, 80), (95, 80)] 
rank_order = ['I', 'II', 'III', 'IV', 'notassessed']

fig4, axes4 = plt.subplots(1, 2, figsize=(14, 5))  # For the Histograms
fig5, axes5 = plt.subplots(1, 2, figsize=(15, 6))  # For the Heatmaps

for i, (p, c) in enumerate(high_thresholds):
    filtered = df[(df.pident >= p) & (df.qcovhsp >= c)]
    
    results = []
    for name, group in filtered.groupby('qseqid'):
        if len(group['rank'].unique()) > 1:
            max_b = group['bitscore'].max()
            
            top_rank = group.loc[group['bitscore'].idxmax(), 'rank']
            diff_hits = group[group['rank'] != top_rank]
            
            if not diff_hits.empty:
                # max_norm_diff = (max_b - diff_hits['bitscore'].min()) / max_b #using normalized difference, i.e., (max bitscore - alternative rank bitscore) / max bitscore
                abs_diff = max_b - diff_hits['bitscore'].max() #using absolute difference, i.e., max bitscore - alternative rank bitscore
                
                riskiest_rank_val = group['risk_score'].min()
                riskiest_rank = {v: k for k, v in risk_map.items()}[riskiest_rank_val]
                
                results.append({'abs_diff': abs_diff,
                    # 'norm_diff': max_norm_diff, 
                                'top_rank': top_rank, 
                                'risk_rank': riskiest_rank})
                
    res_df = pd.DataFrame(results)
    
    if not res_df.empty:
        # sns.histplot(res_df['norm_diff'], ax=axes4[i], bins=30, kde=True, color='blue' if i==0 else 'purple')
        sns.histplot(res_df['abs_diff'], ax=axes4[i], bins=30, kde=True, color='blue' if i==0 else 'purple')
        axes4[i].set_title(f'Absolute Bit-score Diff (Id >= {p}%)', fontsize=14, fontweight='bold')
        # axes4[i].set_title(f'Max Norm Bit-score Diff (Id >= {p}%)', fontsize=14, fontweight='bold')
        
        # axes4[i].set_xlabel('(Max Bitscore - Smaller) / Max Bitscore', fontsize=12, fontweight='bold')
        axes4[i].set_xlabel('Max Bitscore - Alternative Rank Bitscore', fontsize=12, fontweight='bold')
        axes4[i].set_ylabel('Count of Unigenes', fontsize=12, fontweight='bold')

        cm = pd.crosstab(res_df['top_rank'], res_df['risk_rank']).reindex(index=rank_order, columns=rank_order, fill_value=0)

        # Calculate Same vs Changed for the title
        import numpy as np
        same_count = np.trace(cm.values)
        changed_count = cm.values.sum() - same_count

        sns.heatmap(cm, annot=True, fmt='d', cmap='Reds', ax=axes5[i], cbar=(i==1), annot_kws={"size": 12, "weight": "bold"})

        # axes5[i].set_title(f'Rank Shift: Max Bit-score vs Riskiest Hit\n(Id >= {p}%)', fontsize=14, fontweight='bold')
        axes5[i].set_title(f'Rank Shift Heatmap (Id >= {p}%, Cov >= {c}%)\nSame: {same_count} | Changed: {changed_count}', fontsize=14, fontweight='bold')
        # axes5[i].set_ylabel('Rank Assigned by Max Bit-score', fontsize=12, fontweight='bold')
        axes5[i].set_ylabel('Rank by Max Bit-score (Original)', fontweight='bold')
        # axes5[i].set_xlabel('Rank Assigned by Riskiest Hit', fontsize=12, fontweight='bold')
        axes5[i].set_xlabel('Rank by Riskiest Hit (Alternative)', fontweight='bold')
        axes5[i].set_xticklabels(axes5[i].get_xticklabels(), fontweight='bold')
        axes5[i].set_yticklabels(axes5[i].get_yticklabels(), fontweight='bold')

fig4.tight_layout()
# fig4.savefig('normalized_bitscore_diffs.png', dpi=300, bbox_inches='tight')
fig4.savefig('absolute_bitscore_diffs.png', dpi=300, bbox_inches='tight')

fig5.tight_layout()
fig5.savefig('riskiest_hit_rank_shifts.png', dpi=300, bbox_inches='tight')