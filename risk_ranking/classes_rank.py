import pandas as pd
import re
import matplotlib.pyplot as plt

fasta_path = "SARG.db.fasta" 
fasta_data = []

with open(fasta_path, "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            header = line[1:]
            header_parts = header.split()
            
            if not header_parts:
                continue
                
            seq_id = header_parts[0] 
            
            aro_match = re.search(r"ARO:\d+", header)
            aro_id = aro_match.group(0) if aro_match else None
            
            fasta_data.append({"seq_id": seq_id, "aro_id": aro_id})

fasta_df = pd.DataFrame(fasta_data)

# Read argRanker output table 
arg_rank = pd.read_csv("ARG_rank.txt", sep="\t", encoding="utf-8-sig")

# Clean headers and extract embedded ARO IDs
cols = list(arg_rank.columns)
cols[0] = 'ARG'
arg_rank.columns = cols
arg_rank['aro_from_rank'] = arg_rank['ARG'].astype(str).str.extract(r"(ARO:\d+)")

# Merge FASTA mapping with argRanker data
merged_arg = pd.merge(arg_rank, fasta_df, left_on="ARG", right_on="seq_id", how="left")
merged_arg['ARO'] = merged_arg['aro_from_rank'].combine_first(merged_arg['aro_id'])

# Load the mapping file (ARG_classes.csv)
classes_df = pd.read_csv("ARG_classes.csv")
classes_df = classes_df[['ARO_Term_ID', 'abbreviation']]

# Final Annotation Join
# Merge the argRanker data with the custom file based on ARO ID
final_mapped_df = pd.merge(merged_arg, classes_df, left_on="ARO", right_on="ARO_Term_ID", how="left")

# Assign the abbreviation as the Gene_Class. If there is no match in mapping file, label it 'Unknown'
final_mapped_df['Gene_Class'] = final_mapped_df['abbreviation'].fillna("Unknown")

class_rank_summary = pd.crosstab(final_mapped_df['Gene_Class'], final_mapped_df['Rank'])
class_rank_summary['Total Genes'] = class_rank_summary.sum(axis=1)
class_rank_summary = class_rank_summary.sort_values(by='Total Genes', ascending=False)

final_mapped_df.to_csv("ARG_rank_gene_classes.csv", index=False)
class_rank_summary.to_csv("Gene_Class_Rank_Summary.csv")

# Save any unmapped sequences
unmapped_entries = final_mapped_df[final_mapped_df['Gene_Class'] == 'Unknown']
unmapped_entries.to_csv("Missing_ARO_Entries.csv", index=False)

# Generate Stacked Bar Chart
plot_df = class_rank_summary[class_rank_summary.index != 'Unknown']

# Select the top 20 gene classes 
top_plot_df = plot_df.head(20).drop(columns=['Total Genes'], errors='ignore')

# Transpose the dataframe so Ranks are on the X-axis (index) and Gene Classes are the stacked columns
transposed_df = top_plot_df.T

plt.figure(figsize=(12, 8))
ax = transposed_df.plot(kind='bar', stacked=True, colormap='tab20', figsize=(12, 8), edgecolor='black', linewidth=0.5)

# Formatting
plt.title('Gene Classes Distributed Across Risk Ranks', fontsize=16, fontweight='bold', pad=15)
plt.xlabel('Risk Rank', fontsize=12, fontweight='bold')
plt.ylabel('Number of Genes', fontsize=12, fontweight='bold')
plt.xticks(rotation=0, fontsize=11)
plt.yticks(fontsize=10)
plt.legend(title='Gene Class', title_fontsize='11', fontsize='9', bbox_to_anchor=(1.02, 1), loc='upper left')
plt.tight_layout()
plt.savefig('Rank_GeneClass_StackedBar.png', dpi=300)
print("Processing complete. Stacked bar chart saved as 'Rank_GeneClass_StackedBar.png'")
