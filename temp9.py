import polars as pl
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import pandas as pd

# Load data
covid = pl.read_parquet("Results/Covid/node_mutations.parquet")
influenza = pl.read_parquet("Results/Influenza/node_mutations.parquet")

print("=== DATASET OVERVIEW ===")
print(f"Covid: {covid.height} mutations")
print(f"Influenza: {influenza.height} mutations")
print(f"Covid time range: {covid['time_iso'].drop_nulls().min()} to {covid['time_iso'].drop_nulls().max()}")
print(f"Influenza has time data: {influenza['time_iso'].null_count() == 0}")

# === COVID TIME SERIES ANALYSIS ===
print("\n=== COVID TEMPORAL TRENDS ===")

# Filter to coding mutations with valid time (coding = gene is not null)
covid_coding = covid.filter((pl.col("gene").is_not_null()) & (pl.col("time_iso").is_not_null()))
covid_coding = covid_coding.with_columns(
    time_period=pl.col("time_iso").str.slice(0, 7)  # YYYY-MM
)

# 1. Trinucleotide context distribution over time
print("\n1. Trinucleotide Context (96) Distribution:")
ctx_by_time = (
    covid_coding
    .filter(pl.col("ctx96").is_not_null())
    .group_by("time_period", "ctx96")
    .agg(pl.count().alias("count"))
    .sort("time_period")
)
print(f"Total context-annotated mutations: {ctx_by_time['count'].sum()}")
print(f"Unique contexts observed: {ctx_by_time['ctx96'].n_unique()}")

# 2. Codon usage bias (from->to shifts)
print("\n2. Codon Bias Analysis:")
codon_pairs = (
    covid_coding
    .filter((pl.col("codon_from").is_not_null()) & (pl.col("codon_to").is_not_null()) 
            & (pl.col("alt").is_not_null()))  # Exclude deletions (alt=null)
)
print(f"Mutations with codon info: {codon_pairs.height}")

codon_pairs[0].drop("time_iso")

# Count most common codon substitutions
codon_subs = codon_pairs.select(
    (pl.col("codon_from") + "->" + pl.col("codon_to")).alias("substitution")
).to_series().to_list()
codon_counter = Counter(codon_subs)
print("Top 15 codon substitutions (excluding deletions):")
for sub, count in codon_counter.most_common(15):
    print(f"  {sub}: {count}")

# 3. Amino acid substitution patterns
print("\n3. Amino Acid Substitution Patterns:")
aa_pairs = (
    covid_coding
    .filter((pl.col("aa_from").is_not_null()) & (pl.col("aa_to").is_not_null()))
)
print(f"Synonymous mutations: {(covid_coding.filter((pl.col('codon_from') == pl.col('codon_to'))).height)}")
print(f"Non-synonymous mutations: {aa_pairs.height}")

aa_subs = aa_pairs.select(
    (pl.col("aa_from") + "->" + pl.col("aa_to")).alias("substitution")
).to_series().to_list()
aa_counter = Counter(aa_subs)
print("Top 15 amino acid substitutions:")
for sub, count in aa_counter.most_common(15):
    print(f"  {sub}: {count}")

# 4. Mutation rate over time
print("\n4. Mutation Rate Over Time:")
mutation_rate = (
    covid_coding
    .group_by("time_period")
    .agg(pl.count().alias("mutations"))
    .sort("time_period")
)
print(mutation_rate.to_pandas().to_string())

# === VISUALIZATIONS ===
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Plot 1: Mutation count over time
ax = axes[0, 0]
rate_df = mutation_rate.to_pandas()
rate_df['time_period'] = pd.to_datetime(rate_df['time_period'])
ax.plot(rate_df['time_period'], rate_df['mutations'], marker='o', linewidth=2, markersize=4)
ax.set_xlabel('Time')
ax.set_ylabel('Mutations per Month')
ax.set_title('COVID: Mutation Rate Over Time')
ax.grid(True, alpha=0.3)

# Plot 2: Top contexts by time (stacked bar for first/last periods)
ax = axes[0, 1]
ctx_time_df = ctx_by_time.to_pandas()
early = ctx_time_df[ctx_time_df['time_period'] <= '2020-06'].groupby('ctx96')['count'].sum().nlargest(10)
late = ctx_time_df[ctx_time_df['time_period'] >= '2021-01'].groupby('ctx96')['count'].sum().nlargest(10)
x = np.arange(len(early))
width = 0.35
ax.bar(x - width/2, early.values, width, label='Early (≤2020-06)', alpha=0.7)
ax.bar(x + width/2, late.values[:len(early)], width, label='Late (≥2021-01)', alpha=0.7)
ax.set_xlabel('Context')
ax.set_ylabel('Count')
ax.set_title('COVID: Top 10 Trinucleotide Contexts (Early vs Late)')
ax.set_xticks(x)
ax.set_xticklabels(early.index, rotation=45, fontsize=8)
ax.legend()


codon_counter

# Plot 3: Top codon substitutions
ax = axes[1, 0]
top_codons = dict(codon_counter.most_common(15))
subs = list(top_codons.keys())
counts = list(top_codons.values())
ax.barh(range(len(subs)), counts, color='steelblue')
ax.set_yticks(range(len(subs)))
ax.set_yticklabels(subs, fontsize=9)
ax.set_xlabel('Count')
ax.set_title('COVID: Top 15 Codon Substitutions (excl. deletions)')
ax.invert_yaxis()

# Plot 4: Synonymous vs Non-synonymous over time
ax = axes[1, 1]
syn_nonsyn = (
    covid_coding
    .filter(pl.col("time_iso").is_not_null())
    .with_columns(
        time_period=pl.col("time_iso").str.slice(0, 7),
        is_synonymous=(pl.col("codon_from") == pl.col("codon_to"))
    )
    .group_by("time_period", "is_synonymous")
    .agg(pl.count().alias("count"))
    .sort("time_period")
)
syn_nonsyn_df = syn_nonsyn.to_pandas()
pivot_df = syn_nonsyn_df.pivot(index='time_period', columns='is_synonymous', values='count').fillna(0)
pivot_df.index = pd.to_datetime(pivot_df.index)
ax.plot(pivot_df.index, pivot_df[True], marker='o', label='Synonymous', linewidth=2, markersize=4)
ax.plot(pivot_df.index, pivot_df[False], marker='s', label='Non-synonymous', linewidth=2, markersize=4)
ax.set_xlabel('Time')
ax.set_ylabel('Count')
ax.set_title('COVID: Synonymous vs Non-synonymous Mutations Over Time')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('Results/temporal_analysis_covid.png', dpi=300, bbox_inches='tight')
print("\n✓ Saved: Results/temporal_analysis_covid.png")
plt.close()

# === COMPARISON SUMMARY ===
print("\n=== COVID vs INFLUENZA COMPARISON ===")
print(f"Covid coding mutations: {covid.filter(pl.col('gene').is_not_null()).height}")
print(f"Influenza coding mutations: {influenza.filter(pl.col('gene').is_not_null()).height}")
print(f"Covid genes: {covid.filter(pl.col('gene').is_not_null())['gene'].n_unique()}")
print(f"Influenza genes: {influenza.filter(pl.col('gene').is_not_null())['gene'].n_unique()}")

# Gene distribution comparison
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

covid_genes = covid.filter(pl.col('gene').is_not_null()).group_by('gene').agg(pl.count().alias('count')).sort('count', descending=True)
influenza_genes = influenza.filter(pl.col('gene').is_not_null()).group_by('gene').agg(pl.count().alias('count')).sort('count', descending=True)

ax = axes[0]
cg = covid_genes.to_pandas()
ax.barh(range(len(cg)), cg['count'].values, color='coral')
ax.set_yticks(range(len(cg)))
ax.set_yticklabels(cg['gene'].values)
ax.set_xlabel('Mutation Count')
ax.set_title(f'COVID: Mutations by Gene ({len(cg)} genes)')
ax.invert_yaxis()

ax = axes[1]
ig = influenza_genes.to_pandas()
ax.barh(range(len(ig)), ig['count'].values, color='skyblue')
ax.set_yticks(range(len(ig)))
ax.set_yticklabels(ig['gene'].values)
ax.set_xlabel('Mutation Count')
ax.set_title(f'Influenza: Mutations by Gene ({len(ig)} genes)')
ax.invert_yaxis()

plt.tight_layout()
plt.savefig('Results/gene_comparison.png', dpi=300, bbox_inches='tight')
print("✓ Saved: Results/gene_comparison.png")
plt.close()

print("\n=== Analysis complete ===")
print("Visualizations saved to Results/temporal_analysis_covid.png and Results/gene_comparison.png")
