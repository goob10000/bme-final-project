import os
import polars as pl
import numpy as np
import json
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from AnalysisFunctions import CODON_ORDER, adjustTPM, joinMatrixGenesAndCounts, joinMatrixGenesAndCountsWithMT, weightByExpression, weightByAdjustedExpression, diff, compare, JSD
from codon_frequency import GENETIC_CODE
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime
from Bio import SeqIO
from numpy.typing import NDArray
from collections import Counter
import matplotlib.patches as mpatches

InfluenzaDFPath = "influenzaAlignment.csv"

mtGenesList = pl.read_csv("mtGenes.tsv", separator="\t").filter(pl.col("gene_biotype") == "protein_coding")

viruses = ["Covid", "Influenza", "Measles", "EAV", "HIV1", "HIV2", "ChickenPox_VZV", "HSV1", "HSV2", "T4Virus"]
controlViruses = ["Measles", "EAV", "HIV1", "HIV2", "ChickenPox_VZV", "HSV1", "HSV2", "T4Virus"]
controlVirusCBs = []

for v in controlViruses:
    if not os.path.exists(f"Results2/{v}/stuff.npz"):
        raise ValueError(f"Data for {v} not found. Please ensure the stuff.npz file exists at the specified path.")
    s = np.load(f"Results2/{v}/stuff.npz")
    m = s["matrix"]
    controlVirusCBs.append(m.mean(axis=0))

species = [d for d in os.listdir("Results2/") if (d not in viruses)]
tab20b_cmap = plt.colormaps.get_cmap("tab20b")
species_colors = tab20b_cmap(np.linspace(0, 1, max(1, len(species))))
stuff = [np.load(f"Results2/{d}/stuff.npz") for d in species]
exp = [pl.read_csv(f"Results2/{d}/ExpressionData.csv") for d in species]
expression = []
for ex in exp:
    ex.columns = ["gene_id"] + ex.columns[1:]
    ex = ex.select(ex.columns[:2])
    ex.columns = ["gene_id", "expression"]
    expression.append(ex)

dfs = joinMatrixGenesAndCounts(species, stuff, mtGenesList)
dfs_withMT = joinMatrixGenesAndCountsWithMT(species, stuff)

# (dates, codons, genes) = covidCodonBiasMatrix.shape

# Codon biases for all strains. Averaged across all genes in covid.

dfsWithAdjustmentTPM = adjustTPM(dfs, species, expression)
dfsWithAdjustmentAndMT = adjustTPM(dfs_withMT, species, expression)

humanCovid = 4
humanHealthy = 5
humanInfluenza = 6
pig = 11
pigInfluenza = 12

dfs_expressionAdjusted = [weightByAdjustedExpression(df) for df in dfsWithAdjustmentTPM]
dfs_expressionAdjustedWithMT = [weightByAdjustedExpression(df) for df in dfsWithAdjustmentAndMT]
dfs_expressionUnadjusted = [weightByExpression(df) for df in dfsWithAdjustmentTPM]
dfs_expressionUnadjustedWithMT = [weightByExpression(df) for df in dfsWithAdjustmentAndMT]

## Define TPM adjusted expresison weighted codon biases
humanWithCovidCodonBias_TPMAdjusted_expressionWeighted = dfs_expressionAdjusted[humanCovid]
healthyHumanCodonBias_TPMAdjusted_expressionWeighted = dfs_expressionAdjusted[humanHealthy]
humanWithInfluenzaCodonBias_TPMAdjusted_expressionWeighted = dfs_expressionAdjusted[humanInfluenza]
pigCodonBias_TPMAdjusted_expressionWeighted = dfs_expressionAdjusted[pig]
pigWithInfluenzaCodonBias_TPMAdjusted_expressionWeighted = dfs_expressionAdjusted[pigInfluenza]
##Define expression weighted codon biases (no TPM adjustment)
humanWithCovidCodonBias_expressionWeighted = dfs_expressionUnadjusted[humanCovid]
healthyHumanCodonBias_expressionWeighted = dfs_expressionUnadjusted[humanHealthy]
humanWithInfluenzaCodonBias_expressionWeighted = dfs_expressionUnadjusted[humanInfluenza]
pigCodonBias_expressionWeighted = dfs_expressionUnadjusted[pig]
pigWithInfluenzaCodonBias_expressionWeighted = dfs_expressionUnadjusted[pigInfluenza]

n = 100
humanWithCovidCodonBias_expressionOrderedTop100 = dfsWithAdjustmentTPM[humanCovid].sort("adjustedExpression", descending=True).head(n).select(CODON_ORDER).mean()
healthyHumanCodonBias_expressionOrderedTop100 = dfsWithAdjustmentTPM[humanHealthy].sort("adjustedExpression", descending=True).head(n).select(CODON_ORDER).mean()
humanWithInfluenzaCodonBias_expressionOrderedTop100 = dfsWithAdjustmentTPM[humanInfluenza].sort("adjustedExpression", descending=True).head(n).select(CODON_ORDER).mean()
pigCodonBias_expressionOrderedTop100 = dfsWithAdjustmentTPM[pig].sort("adjustedExpression", descending=True).head(n).select(CODON_ORDER).mean()
pigWithInfluenzaCodonBias_expressionOrderedTop100 = dfsWithAdjustmentTPM[pigInfluenza].sort("adjustedExpression", descending=True).head(n).select(CODON_ORDER).mean()

# influenzaDF = pl.read_csv(InfluenzaDFPath)
# a = influenzaDF["Collection_Date"]
# ## Date formats
# # YYYY-MM-DD
# # YYYY-MM
# # YYYY

# ## Forget you know anything about the pandas library, only use polars.

# newDates = []

# for i in influenzaDF["Collection_Date"]:
#     if i == None:
#         newDates.append(None)
#         continue
#     if len(i) == 7:
#         i = i + "-01"
#     elif len(i) == 4:
#         i = i + "-01-01"
#     newDates.append(i)

# influenzaDF = influenzaDF.with_columns(
#     pl.Series("Collection_Date", newDates)
# )

# reference_influenza_genome_path = "Genomes/Influenza/genbank reference sequence/ncbi_dataset/data/GCA_046377675.1/GCA_046377675.1_ASM4637767v1_genomic.fna"
# reference_influenza_gtf_path = "Genomes/Influenza/genbank reference sequence/ncbi_dataset/data/GCA_046377675.1/genomic.gtf"

# gtf = pl.read_csv(reference_influenza_gtf_path, separator="\t", skip_lines=3, new_columns=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
# gtf_cds = gtf.filter(pl.col("feature") == "CDS")
# gtf_cds = gtf_cds.with_columns(
#     gene_id=pl.col("attribute").str.extract(r'gene_id "([^"]+)"')
# )


# ## Example fasta header:
# # >LC660656.1 Influenza A virus (H1N1) A/Iwate/1130/2009 PB2 gene for polymerase PB2, complete cds
# # ATGGAGAGAATAAAAGAACTGAGAGATCTAATGTCGCAGTCCCGCACTCGCGAGATACTCACTAAGACCACTGTGGACCA
# # TATGGCCATAATCAAAAAGTACACATCAGGAAGGCAAGAGAAGAACCCCGCACTCAGAATGAAGTGGATGATGGCAATGA
# # GATACCCAATTACAGCAGACAAGAGAATAATGGACATGATTCCAGAGAGGAATGAACAAGGACAAACCCTCTGGAGCAAA

# # Take cds from gtf_cds. If gene doesn't exist in the dictionary, add it with the sequence.
# # If it does exist, append the sequence.
# # Key for the dictionary is the gene_id which can be found as the first element in the attribute column, which is a string of key-value pairs separated by ";". The gene_id is in the format "gene_id "GENE_NAME"".
# gene_ids = gtf_cds["gene_id"].unique()

# gene_sequences = {}
# for gene_id in gene_ids:
#     gene_sequences[gene_id] = ""

# # Iterate through gtf_cds and append sequences to the corresponding gene_id in gene_sequences
# for row in gtf_cds.iter_rows():
#     gene_id = row[9]  # gene_id is the 10th column (index 9)
#     seqname = row[0]  # seqname is the 1st column (index 0)
#     start = row[3]  # start is the 4th column (index 3)
#     end = row[4]  # end is the 5th column (index 4)
#     strand = row[6]  # strand is the 7th column (index 6)

#     # Extract the sequence from the reference genome
#     for record in SeqIO.parse(reference_influenza_genome_path, "fasta"):
#         if record.id == seqname:
#             if strand == "+":
#                 gene_sequences[gene_id] += str(record.seq[start-1:end+3])  # -1 because of 0-based indexing
#             else:
#                 gene_sequences[gene_id] += str(record.seq[start-1:end+3].reverse_complement())

# gene_lengths = {gene_id: len(seq) for gene_id, seq in gene_sequences.items()}
# length_to_gene = {length: gene_id for gene_id, length in gene_lengths.items()}

# codonPossibilities = json.load(open("codon_possibilities.json", "r"))
# # json.dump(codonPossibilities, open("codon_possibilities.json", "w"))

# influenzaDF = influenzaDF.with_columns(
#     pl.Series([gene_lengths[gene] for gene in influenzaDF["geneID"]]).alias("gene_length")
# )

# codonCounts = influenzaDF[CODON_ORDER].to_numpy() * influenzaDF["gene_length"].to_numpy()[:, None]
# codonProportions = np.zeros_like(codonCounts)
# for i, row in enumerate(codonCounts):
#     if i%100 == 0:
#         print(f"Processing gene {i} of {len(codonCounts)}")
#     c = Counter()
#     for j, codon in enumerate(CODON_ORDER):
#         c[GENETIC_CODE[codon]] += row[j]
#     for j, codon in enumerate(CODON_ORDER):
#         if c[GENETIC_CODE[codon]] > 0:
#             codonProportions[i, j] = row[j] / c[GENETIC_CODE[codon]] / codonPossibilities[codon]

# codonProportions

# influenzaDFNew = influenzaDF.drop(CODON_ORDER).hstack(pl.DataFrame(codonProportions, schema=CODON_ORDER))

# influenzaDFNew.write_parquet("influenzaAlignmentWithCodonProportions.parquet")
influenzaDFNew = pl.read_parquet("influenzaAlignmentWithCodonProportions.parquet").filter(pl.col("geneID") != "M2").filter(pl.col("geneID") != "NEP")

influenzaDFM2 = pl.read_csv("influenzaAlignmentM2Valid.csv")
influenzaDFM2 = influenzaDFM2.insert_column(influenzaDFNew.columns.index("gene_length"),
    pl.Series("gene_length", [len(x) for x in influenzaDFM2["orf"]])
)
influenzaDFM2 = influenzaDFM2.select(influenzaDFNew.columns)

influenzaDFNEP = pl.read_csv("influenzaAlignmentNEPValid.csv")
influenzaDFNEP = influenzaDFNEP.insert_column(influenzaDFNew.columns.index("gene_length"),
    pl.Series("gene_length", [len(x) for x in influenzaDFNEP["orf"]])
)
influenzaDFNEP = influenzaDFNEP.select(influenzaDFNew.columns)

influenzaDFNewWhole = influenzaDFNew.vstack(influenzaDFM2).vstack(influenzaDFNEP)

influenzaDFNewWhole = influenzaDFNewWhole.with_columns(
    pl.Series(np.array(influenzaDFNewWhole["Collection_Date"], dtype="datetime64")).alias("Collection_Date")
)

H1N1 = influenzaDFNewWhole.filter(pl.col("Genotype") == "H1N1")
H1N1_div_humanSick = JSD(H1N1.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type:ignore

covidCodonBiasMatrix = np.load("Results2/Covid/codon_proportions.npy")
covidGenomeDF = pl.read_csv("Genomes/Covid/MT291829_cds.tsv", separator="\t")
covidStuff = np.load("Results2/Covid/stuff.npz")
covidGenes = covidGenomeDF["gene_id"].to_list()
covidMutations = pl.read_parquet("Genomes/Covid/mutations2.parquet")
covidMutations = covidMutations.filter(pl.col("mutations") != covidMutations["mutations"][1]) # Remove null (no mutations?)
covidDates = covidMutations["weeks"]

covidGeneIDS = covidStuff["gene_ids"]
covidMatrix = covidStuff["matrix"]
covidCounts = covidStuff["counts"]
covidGeneOrder = ['ORF1ab', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'N', 'ORF10']


covidAverageDF = pl.DataFrame({"date": covidDates}).hstack(pl.DataFrame({x:covidCodonBiasMatrix.mean(axis=2)[:,i] for i, x in enumerate(CODON_ORDER)}))
# covidDF = covidDF.vstack(pl.DataFrame({gene: }))
covidDFPt1 = pl.DataFrame(schema={"date":pl.Int64, "gene":str})
covidDFPt2 = pl.DataFrame(schema={x: pl.Float64 for x in CODON_ORDER})
covidDF = covidDFPt1.hstack(covidDFPt2)
for i, gene in enumerate(covidGeneOrder):
    geneSubset = covidCodonBiasMatrix[:,:,i]
    dateDFCopyCovid = pl.DataFrame({"date": covidDates})
    geneDF = pl.DataFrame({"gene": [gene] * len(covidDates)})
    codonDFPart = pl.DataFrame({x: geneSubset[:,j] for j, x in enumerate(CODON_ORDER)})
    codonDFPart = dateDFCopyCovid.hstack(geneDF).hstack(codonDFPart)
    covidDF = covidDF.vstack(codonDFPart)


fig0 = plt.figure()
ax = fig0.add_subplot(111)
ax.scatter(covidDates, JSD(covidCodonBiasMatrix.mean(axis=2).T, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]), label="Average Codon Bias Across Genes", s=0.5) #type: ignore
ax.set_xlabel("Weeks Since Week 1 of 2020")
ax.set_ylabel("JSD from Human Infected with Codon Bias")
ax.legend()
fig0.show()

fig01 = plt.figure()
ax = fig01.add_subplot(111)
for gene in covidGeneOrder:
    subset = covidDF.filter(pl.col("gene") == gene)
    jsd_values = JSD(subset.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type: ignore
    ax.scatter(subset["date"], jsd_values, s=0.5, label=gene)
ax.set_xlabel("Weeks Since Week 1 of 2020")
ax.set_ylabel("JSD from Human Infected with Codon Bias")
ax.legend()
fig01.show()

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value

fig02 = plt.figure()
ax = fig02.add_subplot(111)
colors = tab20b_cmap(np.linspace(0, 1, len(covidGeneOrder)))
dateRange = np.arange(0, np.array(covidDates).max()+1, 12) # 4 week intervals
for k, gene in enumerate(covidGeneOrder):
    jsd_values_array = []
    for i, dr1 in enumerate(zip(dateRange[:-1], dateRange[1:])):
        subset = covidDF.filter(
            (pl.col("gene") == gene) &
            (pl.col("date") >= dr1[0]) &
            (pl.col("date") < dr1[1])
        )
        if len(subset) == 0:
            jsd_values_array.append(np.array([np.nan,np.nan,np.nan]))
            continue
        jsd_values = JSD(subset.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type: ignore
        jsd_values_array.append(jsd_values)
    parts = ax.violinplot(jsd_values_array, positions=range(1, len(dateRange)), showmeans=False, showmedians=False, showextrema=False)
    qmqs = [np.percentile(j, [25, 50, 75], axis=0) for j in jsd_values_array]
    quartile1, medians, quartile3 = zip(*qmqs)
    whiskers = np.array([
        adjacent_values(sorted_array, q1, q3)
        for sorted_array, q1, q3 in zip(jsd_values_array, quartile1, quartile3)])
    whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]
    inds = np.arange(1, len(medians) + 1)
    ax.vlines(inds, whiskers_min, whiskers_max, color=colors[k], linestyle='-', lw=2)
    for pc in parts['bodies']: #type: ignore
            pc.set_facecolor(colors[k])
            pc.set_edgecolor(colors[k])
            pc.set_alpha(1)
ax.set_xticks(range(1, len(dateRange)))
ax.set_xticklabels([f"{dr[0]}-{dr[1]}" for dr in zip(dateRange[:-1], dateRange[1:])], rotation=45)
ax.set_xlabel("Weeks Since 2020 by 12 Week Intervals")
ax.set_ylabel("JSD from Human Infected with Codon Bias")
ax.set_title("Distribution of JSD from Human Infected with Codon Bias Over Time for Covid Genes")
# colored legened
handles = [mpatches.Patch(color=colors[i], label=covidGeneOrder[i]) for i in range(len(covidGeneOrder))]
ax.legend(handles=handles, bbox_to_anchor=(1, 1), loc='upper left')
fig02.show()

fig03 = plt.figure()
ax = fig03.add_subplot(111)
dateRange = np.arange(0, np.array(covidDates).max()+1, 12) # 4 week intervals
jsd_values_array = []
jsd_values_array2 = []
for i, dr1 in enumerate(zip(dateRange[:-1], dateRange[1:])):
    subset = covidAverageDF.filter(
        (pl.col("date") >= dr1[0]) &
        (pl.col("date") < dr1[1])
    )
    if len(subset) == 0:
        jsd_values_array.append(np.array([np.nan,np.nan,np.nan]))
        jsd_values_array2.append(np.array([np.nan,np.nan,np.nan]))
        continue
    jsd_values = JSD(subset.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type: ignore
    jsd_values2 = JSD(subset.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_expressionWeighted[:, np.newaxis]) #type: ignore
    jsd_values_array.append(jsd_values)
    jsd_values_array2.append(jsd_values2)
parts = ax.violinplot(jsd_values_array, positions=range(1, len(dateRange)), showmeans=False, showmedians=False, showextrema=False)
qmqs = [np.percentile(j, [25, 50, 75], axis=0) for j in jsd_values_array]
quartile1, medians, quartile3 = zip(*qmqs)  
whiskers = np.array([
    adjacent_values(sorted_array, q1, q3)
    for sorted_array, q1, q3 in zip(jsd_values_array, quartile1, quartile3)])
whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]
inds = np.arange(1, len(medians) + 1)
ax.vlines(inds, whiskers_min, whiskers_max, color='blue', linestyle='-', lw=2)

parts2 = ax.violinplot(jsd_values_array2, positions=range(1, len(dateRange)), showmeans=False, showmedians=False, showextrema=False)
qmqs2 = [np.percentile(j, [25, 50, 75], axis=0) for j in jsd_values_array2]
quartile1_2, medians_2, quartile3_2 = zip(*qmqs2)
whiskers2 = np.array([
    adjacent_values(sorted_array, q1, q3)
    for sorted_array, q1, q3 in zip(jsd_values_array2, quartile1_2, quartile3_2)])
whiskers_min_2, whiskers_max_2 = whiskers2[:, 0], whiskers2[:, 1]
inds_2 = np.arange(1, len(medians_2) + 1)
ax.vlines(inds_2, whiskers_min_2, whiskers_max_2, color='red', linestyle='-', lw=2)

for pc in parts['bodies']: #type: ignore
        pc.set_facecolor('blue')
        pc.set_edgecolor('blue')
        pc.set_alpha(1)
for pc in parts2['bodies']: #type: ignore
        pc.set_facecolor('red')
        pc.set_edgecolor('red')
        pc.set_alpha(1)
ax.set_xticks(range(1, len(dateRange)))
ax.set_xticklabels([f"{dr[0]}-{dr[1]}" for dr in zip(dateRange[:-1], dateRange[1:])], rotation=45)
ax.set_xlabel("Weeks Since 2020 by 12 Week Intervals")
ax.set_ylabel("JSD from Human Infected with Codon Bias")
ax.set_title("Distribution of JSD from Human Infected with Codon Bias Over Time for Covid Genes")
ax.legend(handles=[mpatches.Patch(color='blue', label='All Covid Genes with TPM Adjustment'), mpatches.Patch(color='red', label='Without TPM Adjustment')])
fig03.show()

fig013 = plt.figure()
ax = fig013.add_subplot(111)
dateRange = np.arange(0, np.array(covidDates).max()+1, 12) # 4 week intervals
jsd_values_array = []
jsd_values_array2 = []
for i, dr1 in enumerate(zip(dateRange[:-1], dateRange[1:])):
    subset = covidAverageDF.filter(
        (pl.col("date") >= dr1[0]) &
        (pl.col("date") < dr1[1])
    )
    if len(subset) == 0:
        jsd_values_array.append(np.array([np.nan,np.nan,np.nan]))
        jsd_values_array2.append(np.array([np.nan,np.nan,np.nan]))
        continue
    jsd_values = JSD(subset.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type: ignore
    jsd_values2 = JSD(subset.select(CODON_ORDER).to_numpy().T, healthyHumanCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type: ignore
    jsd_values_array.append(jsd_values)
    jsd_values_array2.append(jsd_values2)
parts = ax.violinplot(jsd_values_array, positions=range(1, len(dateRange)), showmeans=False, showmedians=False, showextrema=False)
qmqs = [np.percentile(j, [25, 50, 75], axis=0) for j in jsd_values_array]
quartile1, medians, quartile3 = zip(*qmqs)  
whiskers = np.array([
    adjacent_values(sorted_array, q1, q3)
    for sorted_array, q1, q3 in zip(jsd_values_array, quartile1, quartile3)])
whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]
inds = np.arange(1, len(medians) + 1)
ax.vlines(inds, whiskers_min, whiskers_max, color='blue', linestyle='-', lw=2)

parts2 = ax.violinplot(jsd_values_array2, positions=range(1, len(dateRange)), showmeans=False, showmedians=False, showextrema=False)
qmqs2 = [np.percentile(j, [25, 50, 75], axis=0) for j in jsd_values_array2]
quartile1_2, medians_2, quartile3_2 = zip(*qmqs2)
whiskers2 = np.array([
    adjacent_values(sorted_array, q1, q3)
    for sorted_array, q1, q3 in zip(jsd_values_array2, quartile1_2, quartile3_2)])
whiskers_min_2, whiskers_max_2 = whiskers2[:, 0], whiskers2[:, 1]
inds_2 = np.arange(1, len(medians_2) + 1)
ax.vlines(inds_2, whiskers_min_2, whiskers_max_2, color='red', linestyle='-', lw=2)

for pc in parts['bodies']: #type: ignore
        pc.set_facecolor('blue')
        pc.set_edgecolor('blue')
        pc.set_alpha(1)
for pc in parts2['bodies']: #type: ignore
        pc.set_facecolor('red')
        pc.set_edgecolor('red')
        pc.set_alpha(1)
ax.set_xticks(range(1, len(dateRange)))
ax.set_xticklabels([f"{dr[0]}-{dr[1]}" for dr in zip(dateRange[:-1], dateRange[1:])], rotation=45)
ax.set_xlabel("Weeks Since 2020 by 12 Week Intervals")
ax.set_ylabel("JSD from Human Infected with Codon Bias")
ax.set_title("Distribution of JSD from Human Infected with Codon Bias Over Time for Covid Genes")
ax.legend(handles=[mpatches.Patch(color='blue', label='Human with covid codon bias'), mpatches.Patch(color='red', label='Healthy human codon bias')])
fig013.show()


fig023 = plt.figure()
ax = fig023.add_subplot(111)
dateRange = np.arange(0, np.array(covidDates).max()+1, 12) # 4 week intervals
jsd_values_array = []
jsd_values_array2 = []
for i, dr1 in enumerate(zip(dateRange[:-1], dateRange[1:])):
    subset = covidAverageDF.filter(
        (pl.col("date") >= dr1[0]) &
        (pl.col("date") < dr1[1])
    )
    if len(subset) == 0:
        jsd_values_array.append(np.array([np.nan,np.nan,np.nan]))
        jsd_values_array2.append(np.array([np.nan,np.nan,np.nan]))
        continue
    jsd_values = JSD(subset.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type: ignore
    jsd_values2 = JSD(subset.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_expressionOrderedTop100.to_numpy().T) #type: ignore
    jsd_values_array.append(jsd_values)
    jsd_values_array2.append(jsd_values2)
parts = ax.violinplot(jsd_values_array, positions=range(1, len(dateRange)), showmeans=False, showmedians=False, showextrema=False)
qmqs = [np.percentile(j, [25, 50, 75], axis=0) for j in jsd_values_array]
quartile1, medians, quartile3 = zip(*qmqs)  
whiskers = np.array([
    adjacent_values(sorted_array, q1, q3)
    for sorted_array, q1, q3 in zip(jsd_values_array, quartile1, quartile3)])
whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]
inds = np.arange(1, len(medians) + 1)
ax.vlines(inds, whiskers_min, whiskers_max, color='blue', linestyle='-', lw=2)

parts2 = ax.violinplot(jsd_values_array2, positions=range(1, len(dateRange)), showmeans=False, showmedians=False, showextrema=False)
qmqs2 = [np.percentile(j, [25, 50, 75], axis=0) for j in jsd_values_array2]
quartile1_2, medians_2, quartile3_2 = zip(*qmqs2)
whiskers2 = np.array([
    adjacent_values(sorted_array, q1, q3)
    for sorted_array, q1, q3 in zip(jsd_values_array2, quartile1_2, quartile3_2)])
whiskers_min_2, whiskers_max_2 = whiskers2[:, 0], whiskers2[:, 1]
inds_2 = np.arange(1, len(medians_2) + 1)
ax.vlines(inds_2, whiskers_min_2, whiskers_max_2, color='red', linestyle='-', lw=2)

for pc in parts['bodies']: #type: ignore
        pc.set_facecolor('blue')
        pc.set_edgecolor('blue')
        pc.set_alpha(1)
for pc in parts2['bodies']: #type: ignore
        pc.set_facecolor('red')
        pc.set_edgecolor('red')
        pc.set_alpha(1)
ax.set_xticks(range(1, len(dateRange)))
ax.set_xticklabels([f"{dr[0]}-{dr[1]}" for dr in zip(dateRange[:-1], dateRange[1:])], rotation=45)
ax.set_xlabel("Weeks Since 2020 by 12 Week Intervals")
ax.set_ylabel("JSD from Human Infected with Codon Bias")
ax.set_title("Distribution of JSD from Human Infected with Codon Bias Over Time for Covid Genes")
ax.legend(handles=[mpatches.Patch(color='blue', label='Expression Weighted'), mpatches.Patch(color='red', label='Expression Ordered Top 100')])
fig023.show()


fig033 = plt.figure()
ax = fig033.add_subplot(111)
dateRange = np.arange(0, np.array(covidDates).max()+1, 12) # 4 week intervals
jsd_values_array = []
jsd_values_array2 = []
for i, dr1 in enumerate(zip(dateRange[:-1], dateRange[1:])):
    subset = covidAverageDF.filter(
        (pl.col("date") >= dr1[0]) &
        (pl.col("date") < dr1[1])
    )
    if len(subset) == 0:
        jsd_values_array.append(np.array([np.nan,np.nan,np.nan]))
        jsd_values_array2.append(np.array([np.nan,np.nan,np.nan]))
        continue
    jsd_values = JSD(subset.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type: ignore
    jsd_values2 = JSD(subset.select(CODON_ORDER).to_numpy().T, dfs_expressionAdjusted[pigInfluenza][:, np.newaxis]) #type: ignore
    jsd_values_array.append(jsd_values)
    jsd_values_array2.append(jsd_values2)
parts = ax.violinplot(jsd_values_array, positions=range(1, len(dateRange)), showmeans=False, showmedians=False, showextrema=False)
qmqs = [np.percentile(j, [25, 50, 75], axis=0) for j in jsd_values_array]
quartile1, medians, quartile3 = zip(*qmqs)  
whiskers = np.array([
    adjacent_values(sorted_array, q1, q3)
    for sorted_array, q1, q3 in zip(jsd_values_array, quartile1, quartile3)])
whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]
inds = np.arange(1, len(medians) + 1)
ax.vlines(inds, whiskers_min, whiskers_max, color='blue', linestyle='-', lw=2)

parts2 = ax.violinplot(jsd_values_array2, positions=range(1, len(dateRange)), showmeans=False, showmedians=False, showextrema=False)
qmqs2 = [np.percentile(j, [25, 50, 75], axis=0) for j in jsd_values_array2]
quartile1_2, medians_2, quartile3_2 = zip(*qmqs2)
whiskers2 = np.array([
    adjacent_values(sorted_array, q1, q3)
    for sorted_array, q1, q3 in zip(jsd_values_array2, quartile1_2, quartile3_2)])
whiskers_min_2, whiskers_max_2 = whiskers2[:, 0], whiskers2[:, 1]
inds_2 = np.arange(1, len(medians_2) + 1)
ax.vlines(inds_2, whiskers_min_2, whiskers_max_2, color='red', linestyle='-', lw=2)

for pc in parts['bodies']: #type: ignore
        pc.set_facecolor('blue')
        pc.set_edgecolor('blue')
        pc.set_alpha(1)
for pc in parts2['bodies']: #type: ignore
        pc.set_facecolor('red')
        pc.set_edgecolor('red')
        pc.set_alpha(1)
ax.set_xticks(range(1, len(dateRange)))
ax.set_xticklabels([f"{dr[0]}-{dr[1]}" for dr in zip(dateRange[:-1], dateRange[1:])], rotation=45)
ax.set_xlabel("Weeks Since 2020 by 12 Week Intervals")
ax.set_ylabel("JSD from SARS-CoV-2")
ax.set_title("Distribution of JSD from SARS-CoV-2 Codon Bias Over Time for pig and human")
ax.legend(handles=[mpatches.Patch(color='blue', label='Human with covid codon bias'), mpatches.Patch(color='red', label='Pig with covid codon bias')])
fig033.show()

fig04 = plt.figure()
ax = fig04.add_subplot(111) 
colors = tab20b_cmap(np.linspace(0, 1, len(species)))
for i, s in enumerate(species):
    codonBiasHere = dfs_expressionAdjusted[i]
    jsd_values = JSD(covidAverageDF[CODON_ORDER].to_numpy().T, codonBiasHere[:, np.newaxis]) #type: ignore
    ax.scatter(covidDates, jsd_values, s=0.5, color=colors[i])
ax.set_xlabel("Weeks Since Week 1 of 2020")
ax.set_ylabel("JSD from Average Covid Codon Bias")
ax.set_title("JSD from Average Covid Codon Bias Over Time for Different Species")
ax.legend(species)
fig04.show()


## Violinplot version of fig04 as fig05. Subtract the baseline JSD value (first JSD value) from all future samples as well to compare changes of all the species side by side.
fig05 = plt.figure(figsize=(15, 10), layout='constrained')
axes = fig05.subplots(len(species)//2, 2, sharex=True, sharey=True).flatten()
colors = tab20b_cmap(np.linspace(0, 1, len(species)))
dateRange = np.arange(0, np.array(covidDates).max()+1, 12) # 12 week intervals
for i, s in enumerate(species):
    codonBiasHere = dfs_expressionAdjusted[i]
    initialJSD = JSD(codonBiasHere[:, np.newaxis], covidAverageDF.select(CODON_ORDER)[0].to_numpy().T) #type: ignore
    jsd_values_array = []
    for j, dr1 in enumerate(zip(dateRange[:-1], dateRange[1:])):
        subset = covidAverageDF.filter(
            (pl.col("date") >= dr1[0]) &
            (pl.col("date") < dr1[1])
        )
        if len(subset) == 0:
            jsd_values_array.append(np.array([np.nan,np.nan,np.nan]))
            continue
        jsd_values = JSD(subset.select(CODON_ORDER).to_numpy().T, codonBiasHere[:, np.newaxis]) #type: ignore
        jsd_values_array.append(jsd_values - initialJSD[0]) # Subtract baseline JSD value
    parts = axes[i].violinplot(jsd_values_array, positions=range(1, len(dateRange)), showmeans=False, showmedians=False, showextrema=False)
    qmqs = [np.percentile(j, [25, 50, 75], axis=0) for j in jsd_values_array]
    quartile1, medians, quartile3 = zip(*qmqs)
    whiskers = np.array([
        adjacent_values(sorted_array, q1, q3)
        for sorted_array, q1, q3 in zip(jsd_values_array, quartile1, quartile3)])
    whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]
    inds = np.arange(1, len(medians) + 1)
    axes[i].vlines(inds, whiskers_min, whiskers_max, color=colors[i], linestyle='-', lw=2)
    for pc in parts['bodies']: #type: ignore
            pc.set_facecolor(colors[i])
            pc.set_edgecolor(colors[i])
            pc.set_alpha(1)
    axes[i].set_title(s)
    if i > 11:
        axes[i].set_xticks(range(1, len(dateRange)))
        axes[i].set_xticklabels([f"{dr[0]}-{dr[1]}" for dr in zip(dateRange[:-1], dateRange[1:])], rotation=45)
        axes[i].set_xlabel("Weeks Since 2020 by 12 Week Intervals")
    if i == 0:
        axes[i].set_ylabel("Change in JSD from Average Covid Codon Bias")
fig05.suptitle("Change in JSD from Average Covid Codon Bias Over Time for Different Species")
fig05.show()

fig005 = plt.figure(figsize=(15, 10), layout='constrained')
ax = fig005.add_subplot(111)
colors = tab20b_cmap(np.linspace(0, 1, len(species)))
dateRange = np.arange(0, np.array(covidDates).max()+1, 12) # 12 week intervals
for i, s in enumerate(species):
    codonBiasHere = dfs_expressionAdjusted[i]
    initialJSD = JSD(codonBiasHere[:, np.newaxis], covidAverageDF.select(CODON_ORDER)[0].to_numpy().T) #type: ignore
    jsd_values_array = []
    for j, dr1 in enumerate(zip(dateRange[:-1], dateRange[1:])):
        subset = covidAverageDF.filter(
            (pl.col("date") >= dr1[0]) &
            (pl.col("date") < dr1[1])
        )
        if len(subset) == 0:
            jsd_values_array.append(np.array([np.nan,np.nan,np.nan]))
            continue
        jsd_values = JSD(subset.select(CODON_ORDER).to_numpy().T, codonBiasHere[:, np.newaxis]) #type: ignore
        jsd_values_array.append(jsd_values - initialJSD[0]) # Subtract baseline JSD value
    parts = ax.violinplot(jsd_values_array, positions=range(1, len(dateRange)), showmeans=False, showmedians=False, showextrema=False)
    qmqs = [np.percentile(j, [25, 50, 75], axis=0) for j in jsd_values_array]
    quartile1, medians, quartile3 = zip(*qmqs)
    whiskers = np.array([
        adjacent_values(sorted_array, q1, q3)
        for sorted_array, q1, q3 in zip(jsd_values_array, quartile1, quartile3)])
    whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]
    inds = np.arange(1, len(medians) + 1)
    ax.vlines(inds, whiskers_min, whiskers_max, color=colors[i], linestyle='-', lw=2)
    for pc in parts['bodies']: #type: ignore
            pc.set_facecolor(colors[i])
            pc.set_edgecolor(colors[i])
            pc.set_alpha(0.2)
    ax.set_title(s)
    if i > 11:
        ax.set_xticks(range(1, len(dateRange)))
        ax.set_xticklabels([f"{dr[0]}-{dr[1]}" for dr in zip(dateRange[:-1], dateRange[1:])], rotation=45)
        ax.set_xlabel("Weeks Since 2020 by 12 Week Intervals")
    if i == 0:
        ax.set_ylabel("Change in JSD from Average Covid Codon Bias")
fig005.suptitle("Change in JSD from Average Covid Codon Bias Over Time for Different Species")
fig005.legend(species, bbox_to_anchor=(1, 1), loc='upper left')
fig005.show()

## Similar concept to fig05. I want a single plot with each species as a line plot.
## The line should follow the same blocks as fig05 but take the mean at each block and produce a line over time.
fig06 = plt.figure(figsize=(15, 10), layout='constrained')
ax = fig06.add_subplot(111)
colors = tab20b_cmap(np.linspace(0, 1, len(species)))
dateRange = np.arange(0, np.array(covidDates).max()+1, 12) # 12 week intervals
for i, s in enumerate(species):
    codonBiasHere = dfs_expressionAdjusted[i]
    initialJSD = JSD(codonBiasHere[:, np.newaxis], covidAverageDF.select(CODON_ORDER)[0].to_numpy().T) #type: ignore
    mean_jsd_values = []
    for j, dr1 in enumerate(zip(dateRange[:-1], dateRange[1:])):
        subset = covidAverageDF.filter(
            (pl.col("date") >= dr1[0]) &
            (pl.col("date") < dr1[1])
        )
        if len(subset) == 0:
            mean_jsd_values.append(np.nan)
            continue
        jsd_values = JSD(subset.select(CODON_ORDER).to_numpy().T, codonBiasHere[:, np.newaxis]) #type: ignore
        mean_jsd_values.append((jsd_values.mean())) # Subtract baseline JSD value
        # mean_jsd_values.append((jsd_values.mean() - initialJSD)[0]) # Subtract baseline JSD value
    ax.plot(range(1, len(dateRange)), mean_jsd_values, color=colors[i], label=s)
ax.set_xticks(range(1, len(dateRange)))
ax.set_xticklabels([f"{dr[0]}-{dr[1]}" for dr in zip(dateRange[:-1], dateRange[1:])], rotation=45)
ax.set_xlabel("Weeks Since 2020 by 12 Week Intervals")
ax.set_ylabel("Change in Mean JSD from Average Covid Codon Bias")
ax.set_title("Change in Mean JSD from Average Covid Codon Bias Over Time for Different Species")
ax.legend()
fig06.show()


## Codon table of each species using Im show and each column being one species and each row being a codon. Color of the cell is the proportion of that codon in the species. Do hierarchical clustering on the rows and columns to group similar codons and species together. Use the same color scheme as before for the species.
fig07 = plt.figure(figsize=(18,6), layout='constrained')
ax = fig07.add_subplot(111)
codonBiasMatrix = np.array([dfs_expressionAdjusted[i] for i in range(len(species))] + covidAverageDF.select(CODON_ORDER)[0].to_numpy().tolist()) #type: ignore
# Perform hierarchical clustering on rows and columns
from scipy.cluster.hierarchy import dendrogram, linkage
row_linkage = linkage(codonBiasMatrix, method='average')
col_linkage = linkage(codonBiasMatrix.T, method='average')
# Get the order of rows and columns after clustering
row_order = dendrogram(row_linkage, no_plot=True)['leaves']
col_order = dendrogram(col_linkage, no_plot=True)['leaves']
# Reorder the codon bias matrix according to the clustering
ordered_codonBiasMatrix = codonBiasMatrix[row_order][:, col_order].T
im = ax.imshow(ordered_codonBiasMatrix.T, aspect='auto', cmap='viridis')
# Set the ticks and labels
ax.set_yticks(range(len(species) + 1))
ax.set_yticklabels([(species + ["covid"])[i] for i in row_order])
ax.set_xticks(range(len(CODON_ORDER)))
ax.set_xticklabels([CODON_ORDER[i] for i in col_order], rotation=90)
ax.set_title("Codon Bias Matrix of Different Species with Hierarchical Clustering")
fig07.colorbar(im, ax=ax)
fig07.savefig("codon_bias_matrix_clustering.T.png")
fig07.show()

JSD()
JSD(np.array([0,1]), np.array([1,0]))


fig001 = plt.figure()
ax = fig001.add_subplot(111)
ax.scatter(np.array(H1N1["Collection_Date"], dtype="datetime64"), JSD(H1N1.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]), s=0.5) #type: ignore
ax.scatter(np.array(H1N1["Collection_Date"], dtype="datetime64"), JSD(H1N1.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_expressionWeighted[:, np.newaxis]), s=0.5) #type: ignore
ax.legend(["TPM Adjusted", "Expression Weighted"])
ax.set_xlabel("Collection Date")
ax.set_ylabel("JSD from Sick Human Codon Bias")
fig001.show()

fig011 = plt.figure()


fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.scatter(np.array(H1N1["Collection_Date"], dtype="datetime64"), H1N1_div_humanSick, s=0.5)
ax.set_xlabel("Collection Date")
ax.set_ylabel("JSD from Sick Human Codon Bias")
fig1.show()

fig2 = plt.figure()
ax = fig2.add_subplot(111)
for genotype in influenzaDFNewWhole["Genotype"].unique():
    subset2 = influenzaDFNewWhole.filter(pl.col("Genotype") == genotype)
    jsd_values = JSD(subset2.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type:ignore
    ax.scatter(np.array(subset2["Collection_Date"], dtype="datetime64"), jsd_values, s=0.5, label=genotype)
ax.set_xlabel("Collection Date")
ax.set_ylabel("JSD from Sick Human Codon Bias")
ax.legend()
fig2.show()

fig3 = plt.figure()
ax = fig3.add_subplot(111)
for gene in H1N1["geneID"].unique():
    subset2 = H1N1.filter(pl.col("geneID") == gene)
    jsd_values = JSD(subset2.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type:ignore
    ax.scatter(np.array(subset2["Collection_Date"], dtype="datetime64"), jsd_values, s=0.5, label=gene)
ax.set_xlabel("Collection Date")
ax.set_ylabel("JSD from Sick Human Codon Bias")
ax.legend()
fig3.show()

fig4 = plt.figure()
axes = fig4.subplots(5, 2, sharex=True, sharey=True)
for i, gene in enumerate(H1N1["geneID"].unique()):
    subset2 = H1N1.filter(pl.col("geneID") == gene)
    jsd_values = JSD(subset2.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type:ignore
    ax = axes.flatten()[i]
    ax.scatter(np.array(subset2["Collection_Date"], dtype="datetime64"), jsd_values, s=0.5)
    ax.set_title(gene)
    ax.set_xlabel("Collection Date")
    ax.set_ylabel("Codon Bias JSD")
    # ax.set_xlim([np.datetime64("2000-01-01"), np.datetime64("2025-01-01")])
fig4.suptitle("JSD from Sick Human Codon Bias Over Time for H1N1 Genes")
fig4.show()

# violin plot version of fig4 as fig_4
fig_4 = plt.figure(figsize=(20, 5))
ax = fig_4.add_subplot(111)
dateRange = pl.date_range(start=pl.date(2000, 1, 1), end=pl.date(2026, 4, 1), interval="3mo", eager=True)
colors = tab20b_cmap(np.linspace(0, 1, len(H1N1["geneID"].unique())))
for i, gene in enumerate(H1N1["geneID"].unique()):
    jsd_values_array = []
    for j, dr in enumerate(zip(dateRange[:-1], dateRange[1:])):
        subset2 = H1N1.filter(
            (pl.col("geneID") == gene) &
            (pl.col("Collection_Date") >= dr[0]) &
            (pl.col("Collection_Date") < dr[1])
        )
        if len(subset2) == 0:
            jsd_values_array.append(np.array([np.nan,np.nan,np.nan]))
            continue
        jsd_values = JSD(subset2.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type:ignore
        jsd_values_array.append(jsd_values)
    parts = ax.violinplot(jsd_values_array, positions=range(0, len(jsd_values_array)), showmeans=False, showmedians=False, showextrema=False)
    qmqs = [np.percentile(j, [25, 50, 75], axis=0) for j in jsd_values_array]
    quartile1, medians, quartile3 = zip(*qmqs)
    whiskers = np.array([
        adjacent_values(sorted_array, q1, q3)
        for sorted_array, q1, q3 in zip(jsd_values_array, quartile1, quartile3)
    ])
    whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]
    inds = np.arange(1, len(medians) + 1)
    ax.vlines(inds, whiskers_min, whiskers_max, color=colors[i], linestyle='-', lw=2)
    for pc in parts['bodies']: #type: ignore
        pc.set_facecolor(colors[i])
        pc.set_edgecolor(colors[i])
        pc.set_alpha(1)
    ax.set_title(gene)
    ax.set_xlabel("Collection Date")
    ax.set_ylabel("Codon Bias JSD")
    ax.set_xticks(range(0, len(jsd_values_array)))
    # ax.set_xticklabels([f"{dateRange[i][0].year}-{dateRange[i][0].month:02d}" for i in range(len(dateRange)-1)], rotation=45)
fig_4.legend(handles=[mpatches.Patch(color=colors[i], label=H1N1["geneID"].unique()[i]) for i in range(len(H1N1["geneID"].unique()))], bbox_to_anchor=(1, 1), loc='upper left')
fig_4.suptitle("Distribution of JSD from Sick Human Codon Bias Over Time for H1N1 Genes")
fig_4.show()

H1N1[2]["TAC", "TAT"].sum_horizontal()
covidCodonBiasMatrix.mean(axis=2)[0,:][48:50].sum() - covidCodonBiasMatrix.mean(axis=2)[0,:][49]


InfluenzaGenes = ['M2', 'NA', 'PA', 'PB2', 'PB1', 'NEP', 'HA', 'NS1', 'M1', 'NP']
InfluenzaGenes = ['NEP', 'HA', 'NS1', 'M1', 'NP']
CODON_ORDER.index("TAT")
humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[48] + humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[50]

for gene in InfluenzaGenes:
    dates = pl.date_range(start=pl.date(2003, 1, 1), end=pl.date(2013, 4, 1), interval="3mo", eager=True)
    fig = plt.figure(figsize=(20, 5))
    ax = fig.add_subplot(111)
    jsd_values_array = []
    jsd_values_array2 = []
    for i, dateRange in enumerate(zip(dates[:-1], dates[1:])):
        subset2 = H1N1.filter(
            (pl.col("geneID") == gene) &
            (pl.col("Collection_Date") >= dateRange[0]) &
            (pl.col("Collection_Date") < dateRange[1])
        )
        if len(subset2) == 0:
            jsd_values_array.append(np.array([np.nan,np.nan, np.nan]))
            jsd_values_array2.append(np.array([np.nan,np.nan, np.nan]))
            continue
        jsd_values = JSD(subset2.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type:ignore
        jsd_values2 = JSD(subset2.select(CODON_ORDER).to_numpy().T, dfs_expressionAdjusted[pigInfluenza][:, np.newaxis]) #type:ignore
        jsd_values_array.append(jsd_values)
        jsd_values_array2.append(jsd_values2)
    parts = ax.violinplot(jsd_values_array, positions=range(1, len(dates)), showmeans=False, showmedians=False, showextrema=False)
    qmqs = [np.percentile(j, [25, 50, 75], axis=0) for j in jsd_values_array]
    quartile1, medians, quartile3 = zip(*qmqs)  
    whiskers = np.array([
        adjacent_values(sorted_array, q1, q3)
        for sorted_array, q1, q3 in zip(jsd_values_array, quartile1, quartile3)])
    whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]
    inds = np.arange(1, len(medians) + 1)
    ax.vlines(inds, whiskers_min, whiskers_max, color='blue', linestyle='-', lw=2)

    parts2 = ax.violinplot(jsd_values_array2, positions=range(1, len(dates)), showmeans=False, showmedians=False, showextrema=False)
    qmqs2 = [np.percentile(j, [25, 50, 75], axis=0) for j in jsd_values_array2]
    quartile1_2, medians_2, quartile3_2 = zip(*qmqs2)
    whiskers2 = np.array([
        adjacent_values(sorted_array, q1, q3)
        for sorted_array, q1, q3 in zip(jsd_values_array2, quartile1_2, quartile3_2)])
    whiskers_min_2, whiskers_max_2 = whiskers2[:, 0], whiskers2[:, 1]
    inds_2 = np.arange(1, len(medians_2) + 1)
    ax.vlines(inds_2, whiskers_min_2, whiskers_max_2, color='red', linestyle='-', lw=2)

    for pc in parts['bodies']: #type: ignore
            pc.set_facecolor('blue')
            pc.set_edgecolor('blue')
            pc.set_alpha(0.5)
    for pc in parts2['bodies']: #type: ignore
            pc.set_facecolor('red')
            pc.set_edgecolor('red')
            pc.set_alpha(0.5)
    ax.set_xticks(range(1, len(dates)))
    ax.set_xticklabels([f"{dates[i]}-{dates[i+1]}" for i in range(len(dates)-1)], rotation=45)
    ax.set_xlabel("Date")
    ax.set_ylabel("JSD from Influenza Codon Bias")
    ax.set_title(f"Gene: {gene} - Distribution of JSD from SARS-CoV-2 Codon Bias Over Time for pig and human")
    ax.legend(handles=[mpatches.Patch(color='blue', label='Human with covid codon bias'), mpatches.Patch(color='red', label='Pig with covid codon bias')])
    fig.show()


for gene in InfluenzaGenes:
    subset2 = H1N1.filter(pl.col("geneID") == gene).filter(pl.col("Collection_Date") > pl.date(2000, 1, 1))
    jsd_values = JSD(subset2.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type:ignore
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(np.array(subset2["Collection_Date"], dtype="datetime64"), jsd_values, s=0.5)
    ax.set_title(gene)
    ax.set_xlabel("Collection Date")
    ax.set_ylabel("JSD between Influenza Gene and Human Infected with Codon Bias")
    fig.savefig("InfluenzaMedia/H1N1/GeneJSD_2000_Onward/" + gene + "_JSD_over_time.png")


for gene in InfluenzaGenes:
    dates = pl.date_range(start=pl.date(2000, 1, 1), end=pl.date(2026, 4, 1), interval="3mo", eager=True)
    fig = plt.figure(figsize=(40, 10))
    ax = fig.add_subplot(111)
    jsd_values_array = []
    for i, dateRange in enumerate(zip(dates[:-1], dates[1:])):
        subset2 = H1N1.filter(
            (pl.col("geneID") == gene) &
            (pl.col("Collection_Date") >= dateRange[0]) &
            (pl.col("Collection_Date") < dateRange[1])
        )
        if len(subset2) == 0:
            jsd_values_array.append(np.array([0,0,0]))
            continue
        jsd_values = JSD(subset2.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type:ignore
        jsd_values_array.append(jsd_values)
    ax.violinplot(jsd_values_array, positions=range(1, len(dates)))
    ax.set_xticks(range(1, len(dates)))
    ax.set_xticklabels([f"{dateRange[0].year}-{dateRange[0].month:02d}" for dateRange in zip(dates[:-1], dates[1:])], rotation=45)
    ax.set_title(gene)
    ax.set_xlabel("Collection Date")
    ax.set_ylabel("JSD between Influenza Gene and Human Infected with Influenza Codon Bias")
    fig.savefig("InfluenzaMedia/H1N1/GeneJSD_Violin_2000_Onward/" + gene + "_JSD_over_time.png")

for gene in InfluenzaGenes:
    dates = pl.date_range(start=pl.date(2000, 1, 1), end=pl.date(2026, 4, 1), interval="3mo", eager=True)
    fig = plt.figure(figsize=(40, 10))
    ax = fig.add_subplot(111)
    jsd_values_array1 = []
    jsd_values_array2 = []
    for i, dateRange in enumerate(zip(dates[:-1], dates[1:])):
        subset2 = H1N1.filter(
            (pl.col("geneID") == gene) &
            (pl.col("Collection_Date") >= dateRange[0]) &
            (pl.col("Collection_Date") < dateRange[1])
        )
        if len(subset2) == 0:
            jsd_values_array1.append(np.array([0,0,0]))
            jsd_values_array2.append(np.array([0,0,0]))
            continue
        jsd_values1 = JSD(subset2.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type:ignore
        jsd_values2 = JSD(subset2.select(CODON_ORDER).to_numpy().T, healthyHumanCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type:ignore
        jsd_values_array1.append(jsd_values1)
        jsd_values_array2.append(jsd_values2)
    ax.violinplot(jsd_values_array1, positions=range(1, len(dates)))
    ax.violinplot(jsd_values_array2, positions=range(1, len(dates)))

    ax.set_xticks(range(1, len(dates)))
    ax.set_xticklabels([f"{dateRange[0].year}-{dateRange[0].month:02d}" for dateRange in zip(dates[:-1], dates[1:])], rotation=45)
    ax.set_title(gene)
    ax.set_xlabel("Collection Date")
    ax.set_ylabel("JSD between Influenza Gene and Human Infected with Influenza Codon Bias")
    fig.savefig("InfluenzaMedia/H1N1/GeneJSD_Violin_2000_Onward_ExpressionComparison/" + gene + "_JSD_over_time.png")


for gene in InfluenzaGenes:
    dates = pl.date_range(start=pl.date(2008, 1, 1), end=pl.date(2010, 4, 1), interval="1w", eager=True)
    fig = plt.figure(figsize=(40, 10))
    ax = fig.add_subplot(111)
    jsd_values_array = []
    for i, dateRange in enumerate(zip(dates[:-1], dates[1:])):
        subset2 = H1N1.filter(
            (pl.col("geneID") == gene) &
            (pl.col("Collection_Date") >= dateRange[0]) &
            (pl.col("Collection_Date") < dateRange[1])
        )
        if len(subset2) == 0:
            jsd_values_array.append(np.array([0,0,0]))
            continue
        jsd_values = JSD(subset2.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type:ignore
        jsd_values_array.append(jsd_values)
    ax.violinplot(jsd_values_array, positions=range(1, len(dates)))
    ax.set_xticks(range(1, len(dates)))
    ax.set_xticklabels([f"{dateRange[0].year}-{dateRange[0].month:02d}" for dateRange in zip(dates[:-1], dates[1:])], rotation=45)
    ax.set_title(gene)
    ax.set_xlabel("Collection Date")
    ax.set_ylabel("JSD between Influenza Gene and Human Infected with Influenza Codon Bias")
    fig.savefig("InfluenzaMedia/H1N1/GeneJSD_Violin_2008-2010/" + gene + "_JSD_over_time.png")

for gene in InfluenzaGenes:
    dates = pl.date_range(start=pl.date(2008, 1, 1), end=pl.date(2010, 4, 1), interval="1mo", eager=True)
    fig = plt.figure(figsize=(40, 10))
    ax = fig.add_subplot(111)
    jsd_values_array1 = []
    jsd_values_array2 = []
    for i, dateRange in enumerate(zip(dates[:-1], dates[1:])):
        subset2 = H1N1.filter(
            (pl.col("geneID") == gene) &
            (pl.col("Collection_Date") >= dateRange[0]) &
            (pl.col("Collection_Date") < dateRange[1])
        )
        if len(subset2) == 0:
            jsd_values_array1.append(np.array([0,0,0]))
            jsd_values_array2.append(np.array([0,0,0]))
            continue
        jsd_values1 = JSD(subset2.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type:ignore
        jsd_values2 = JSD(subset2.select(CODON_ORDER).to_numpy().T, healthyHumanCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type:ignore
        jsd_values_array1.append(jsd_values1)
        jsd_values_array2.append(jsd_values2)
    ax.violinplot(jsd_values_array1, positions=range(1, len(dates)))
    ax.violinplot(jsd_values_array2, positions=range(1, len(dates)))
    ax.set_xticks(range(1, len(dates)))
    ax.set_xticklabels([f"{dateRange[0].year}-{dateRange[0].month:02d}" for dateRange in zip(dates[:-1], dates[1:])], rotation=45)
    ax.set_title(gene)
    ax.set_xlabel("Collection Date")
    ax.set_ylabel("JSD between Influenza Gene and Human Infected with Influenza Codon Bias")
    fig.savefig("InfluenzaMedia/H1N1/GeneJSD_Violin_2008-2010_ExpressionComparison/" + gene + "_JSD_over_time.png")

for gene in InfluenzaGenes:
    fig = plt.figure(figsize=(36*4, 36))
    axes = [fig.add_subplot(8, 8, i) for i in range(1, 65)]
    dates = pl.date_range(start=pl.date(2000, 1, 1), end=pl.date(2026, 4, 1), interval="3mo", eager=True)
    subset = H1N1.filter(pl.col("geneID") == gene).filter(pl.col("Collection_Date") > pl.date(2000, 1, 1))
    for j, codon in enumerate(CODON_ORDER):
        stuff = []
        for i, dateRange in enumerate(zip(dates[:-1], dates[1:])):
            subset2 = subset.filter(
                (pl.col("Collection_Date") >= dateRange[0]) &
                (pl.col("Collection_Date") < dateRange[1])
            )
            if len(subset2) == 0:
                stuff.append(np.array([0,0,0]))
                continue
            stuff.append(subset2.select(codon).to_numpy().flatten().T)
        # print("hi")
        parts = axes[j].violinplot(stuff, positions=range(1, len(dates)))
        # print("hihi")
        for pc in parts['bodies']: #type: ignore
            pc.set_facecolor('#D43F3A')
            pc.set_edgecolor('black')
            pc.set_alpha(1)
        axes[j].set_title(codon)
        axes[j].legend()
        axes[j].set_xlabel("Collection Date")
        axes[j].set_ylabel("Codon Fraction")
        axes[j].set_xticks(range(1, len(dates)))
        axes[j].set_xticklabels([f"{dateRange[0].year}-{dateRange[0].month:02d}" for dateRange in zip(dates[:-1], dates[1:])], rotation=45)
    fig.savefig("InfluenzaMedia/H1N1/GeneCodonJSD_2000_Onward/" + gene + "_JSD_over_time.png", dpi=300)

for gene in InfluenzaGenes:
    subset2 = H1N1.filter(pl.col("geneID") == gene).filter(pl.col("Collection_Date") > pl.date(2008, 1, 1)).filter(pl.col("Collection_Date") < pl.date(2010, 4, 1))
    fig = plt.figure(figsize=(24, 24))
    axes = [fig.add_subplot(8, 8, i) for i in range(1, 65)]
    for j, codon in enumerate(CODON_ORDER):
        axes[j].scatter(np.array(subset2["Collection_Date"], dtype="datetime64"), subset2[codon], s=0.5, label=f"{codon} fraction")
        axes[j].set_title(codon)
        axes[j].legend()
        axes[j].set_xlabel("Collection Date")
        axes[j].set_ylabel("Codon Fraction")
    fig.savefig("InfluenzaMedia/H1N1/GeneCodonJSD_2008-2010/" + gene + "_JSD_over_time.png", dpi=300)

for gene in InfluenzaGenes:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(range(len(CODON_ORDER)), H1N1.filter(pl.col("geneID") == gene).select(CODON_ORDER).std().to_numpy()[0])
    ax.set_xticks(range(len(CODON_ORDER)))
    ax.set_xticklabels(CODON_ORDER, rotation=45)
    ax.set_title(f"Standard Deviation of Codon Fractions for {gene} Gene")
    ax.set_xlabel("Codon")
    ax.set_ylabel("Standard Deviation")
    fig.savefig(f"InfluenzaMedia/H1N1/GeneCodonJSD_2000_Onward/{gene}_Codon_Fraction_Std.png", dpi=300)

## Need to redo for H1N1
geneCuttoffSTDs0 = {
    "M2": 0.1,
    "NA": 0.08,
    "PA": 0.04,
    "PB2": 0.1,
    "PB1": 0.1,
    "NEP": 0.1,
    "HA": 0.1,
    "NS1": 0.08,
    "M1": 0.1,
    "NP": 0.1
}

for gene in InfluenzaGenes:
    codons = H1N1.filter(pl.col("geneID") == gene).filter(pl.col("Collection_Date") > pl.date(2000, 1, 1)).select(CODON_ORDER).std().to_numpy()[0]
    cutoff = geneCuttoffSTDs0[gene]
    codonsToInclude = [CODON_ORDER[i] for i, std in enumerate(codons) if std > cutoff]
    humanBias = humanWithCovidCodonBias_expressionWeighted[:, np.newaxis] #type:ignore
    # pick out codons from humanBias that are in codonsToInclude
    humanBiasSubset = np.array([humanBias[i] for i, codon in enumerate(CODON_ORDER) if codon in codonsToInclude])
    subset2 = influenzaDFNewWhole.filter(pl.col("geneID") == gene).filter(pl.col("Collection_Date") > pl.date(2000, 1, 1))
    jsd_values = JSD(subset2.select(codonsToInclude).to_numpy().T, humanBiasSubset) 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(np.array(subset2["Collection_Date"], dtype="datetime64"), jsd_values, s=0.5)
    ax.set_title(gene)
    ax.set_xlabel("Collection Date")
    ax.set_ylabel("JSD between Influenza Gene and Human Infected with Codon Bias")
    fig.savefig("InfluenzaMedia/H1N1/GeneJSD_CodonSpecific_highBar_2000_onward/" + gene + "_JSD_over_time.png")


# Need to redo for H1N1
geneCuttoffSTDs1 = {
    "M2": 0.05,
    "NA": 0.05,
    "PA": 0.03,
    "PB2": 0.04,
    "PB1": 0.06,
    "NEP": 0.05,
    "HA": 0.1,
    "NS1": 0.05,
    "M1": 0.06,
    "NP": 0.04
}

for gene in InfluenzaGenes:
    codons = H1N1.filter(pl.col("geneID") == gene).filter(pl.col("Collection_Date") > pl.date(2000, 1, 1)).select(CODON_ORDER).std().to_numpy()[0]
    cutoff = geneCuttoffSTDs1[gene]
    codonsToInclude = [CODON_ORDER[i] for i, std in enumerate(codons) if std > cutoff]
    humanBias = humanWithCovidCodonBias_expressionWeighted[:, np.newaxis] #type:ignore
    # pick out codons from humanBias that are in codonsToInclude
    humanBiasSubset = np.array([humanBias[i] for i, codon in enumerate(CODON_ORDER) if codon in codonsToInclude])
    subset2 = H1N1.filter(pl.col("geneID") == gene).filter(pl.col("Collection_Date") > pl.date(2000, 1, 1))
    jsd_values = JSD(subset2.select(codonsToInclude).to_numpy().T, humanBiasSubset) 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(np.array(subset2["Collection_Date"], dtype="datetime64"), jsd_values, s=0.5)
    ax.set_title(gene)
    ax.set_xlabel("Collection Date")
    ax.set_ylabel("JSD between Influenza Gene and Human Infected with Codon Bias")
    fig.savefig("InfluenzaMedia/H1N1/GeneJSD_CodonSpecific_medBar_2000_onward/" + gene + "_JSD_over_time.png")


# Need to redo for H1N1
geneCuttoffSTDs2 = {
    "M2": 0.02,
    "NA": 0.02,
    "PA": 0.01,
    "PB2": 0.02,
    "PB1": 0.02,
    "NEP": 0.025,
    "HA": 0.02,
    "NS1": 0.02,
    "M1": 0.03,
    "NP": 0.02
}


for gene in InfluenzaGenes:
    codons = H1N1.filter(pl.col("geneID") == gene).filter(pl.col("Collection_Date") > pl.date(2000, 1, 1)).select(CODON_ORDER).std().to_numpy()[0]
    cutoff = geneCuttoffSTDs2[gene]
    codonsToInclude = [CODON_ORDER[i] for i, std in enumerate(codons) if std > cutoff]
    humanBias = humanWithCovidCodonBias_expressionWeighted[:, np.newaxis] #type:ignore
    # pick out codons from humanBias that are in codonsToInclude
    humanBiasSubset = np.array([humanBias[i] for i, codon in enumerate(CODON_ORDER) if codon in codonsToInclude])
    subset2 = H1N1.filter(pl.col("geneID") == gene).filter(pl.col("Collection_Date") > pl.date(2000, 1, 1))
    jsd_values = JSD(subset2.select(codonsToInclude).to_numpy().T, humanBiasSubset) 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(np.array(subset2["Collection_Date"], dtype="datetime64"), jsd_values, s=0.5)
    ax.set_title(gene)
    ax.set_xlabel("Collection Date")
    ax.set_ylabel("JSD between Influenza Gene and Human Infected with Codon Bias")
    fig.savefig("InfluenzaMedia/H1N1/GeneJSD_CodonSpecific_lowBar_2000_onward/" + gene + "_JSD_over_time.png")




codons_of_interest = ["AAT", "AAG"]
codons_of_interest = ["AAT"]
codons_of_interest = ["AAG"]
codons_of_interest = CODON_ORDER
codon_indices = [CODON_ORDER.index(codon) for codon in codons_of_interest]
gene_of_interest = "NA"
subset = H1N1.filter(pl.col("geneID") == gene_of_interest).select(codons_of_interest + ["Collection_Date"]).filter(pl.col("Collection_Date") > pl.date(2000, 1, 1))
dates = pl.date_range(start=pl.date(2000, 1, 1), end=pl.date(2026, 4, 1), interval="3mo", eager=True)

jsd_values_array1 = []
jsd_values_array2 = []
for i, dateRange in enumerate(zip(dates[:-1], dates[1:])):
    subset2 = subset.filter(
        (pl.col("Collection_Date") >= dateRange[0]) &
        (pl.col("Collection_Date") < dateRange[1])
    )
    if len(subset2) == 0:
        jsd_values_array1.append(np.array([0,0,0]))
        jsd_values_array2.append(np.array([0,0,0]))
        continue
    jsd_values1 = JSD(subset2.select(codons_of_interest).to_numpy().T, humanWithInfluenzaCodonBias_expressionWeighted[codon_indices][:, np.newaxis]) #type:ignore
    jsd_values2 = JSD(subset2.select(codons_of_interest).to_numpy().T, pigWithInfluenzaCodonBias_expressionWeighted[codon_indices][:, np.newaxis]) #type:ignore
    jsd_values_array1.append(jsd_values1)
    jsd_values_array2.append(jsd_values2)

fig = plt.figure()
ax = fig.add_subplot(111)
labels = []
def add_label(violin, label, color):
    color = violin["bodies"][0].get_facecolor().flatten()
    labels.append((mpatches.Patch(color=color), label)) #type: ignore
    for pc in violin['bodies']: #type: ignore
        pc.set_facecolor(color)
        # pc.set_edgecolor('black')
        pc.set_alpha(1)
add_label(ax.violinplot(jsd_values_array1, positions=range(1, len(dates)), showextrema=False), "human", "#0008FF")
add_label(ax.violinplot(jsd_values_array2, positions=range(1, len(dates)), showextrema=False), "pig", "#D43F3A")
ax.legend(*zip(*labels), loc=2)

ax.set_xticks(range(1, len(dates)))
ax.set_xticklabels([f"{dateRange[0].year}-{dateRange[0].month:02d}" for dateRange in zip(dates[:-1], dates[1:])], rotation=45)
ax.set_xlabel("Collection Date")
ax.set_ylabel("JSD between Influenza Gene and Human Infected with Influenza Codon Bias")
fig.show()



## All of influenza plots



for gene in InfluenzaGenes:
    subset2 = influenzaDFNewWhole.filter(pl.col("geneID") == gene).filter(pl.col("Collection_Date") > pl.date(2000, 1, 1))
    jsd_values = JSD(subset2.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_expressionWeighted[:, np.newaxis]) #type:ignore
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(np.array(subset2["Collection_Date"], dtype="datetime64"), jsd_values, s=0.5)
    ax.set_title(gene)
    ax.set_xlabel("Collection Date")
    ax.set_ylabel("JSD between Influenza Gene and Human Infected with Codon Bias")
    fig.savefig("InfluenzaMedia/GeneJSD_2000_Onward/" + gene + "_JSD_over_time.png")


for gene in InfluenzaGenes:
    dates = pl.date_range(start=pl.date(2000, 1, 1), end=pl.date(2026, 4, 1), interval="3mo", eager=True)
    fig = plt.figure(figsize=(40, 10))
    ax = fig.add_subplot(111)
    jsd_values_array = []
    for i, dateRange in enumerate(zip(dates[:-1], dates[1:])):
        subset2 = influenzaDFNewWhole.filter(
            (pl.col("geneID") == gene) &
            (pl.col("Collection_Date") >= dateRange[0]) &
            (pl.col("Collection_Date") < dateRange[1])
        )
        if len(subset2) == 0:
            jsd_values_array.append(np.array([0,0,0]))
            continue
        jsd_values = JSD(subset2.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_expressionWeighted[:, np.newaxis]) #type:ignore
        jsd_values_array.append(jsd_values)
    ax.violinplot(jsd_values_array, positions=range(1, len(dates)))
    ax.set_xticks(range(1, len(dates)))
    ax.set_xticklabels([f"{dateRange[0].year}-{dateRange[0].month:02d}" for dateRange in zip(dates[:-1], dates[1:])], rotation=45)
    ax.set_title(gene)
    ax.set_xlabel("Collection Date")
    ax.set_ylabel("JSD between Influenza Gene and Human Infected with Influenza Codon Bias")
    fig.savefig("InfluenzaMedia/GeneJSD_Violin_2000_Onward/" + gene + "_JSD_over_time.png")

for gene in InfluenzaGenes:
    dates = pl.date_range(start=pl.date(2000, 1, 1), end=pl.date(2026, 4, 1), interval="3mo", eager=True)
    fig = plt.figure(figsize=(40, 10))
    ax = fig.add_subplot(111)
    jsd_values_array1 = []
    jsd_values_array2 = []
    for i, dateRange in enumerate(zip(dates[:-1], dates[1:])):
        subset2 = influenzaDFNewWhole.filter(
            (pl.col("geneID") == gene) &
            (pl.col("Collection_Date") >= dateRange[0]) &
            (pl.col("Collection_Date") < dateRange[1])
        )
        if len(subset2) == 0:
            jsd_values_array1.append(np.array([0,0,0]))
            jsd_values_array2.append(np.array([0,0,0]))
            continue
        jsd_values1 = JSD(subset2.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_expressionWeighted[:, np.newaxis]) #type:ignore
        jsd_values2 = JSD(subset2.select(CODON_ORDER).to_numpy().T, healthyHumanCodonBias_expressionWeighted[:, np.newaxis]) #type:ignore
        jsd_values_array1.append(jsd_values1)
        jsd_values_array2.append(jsd_values2)
    ax.violinplot(jsd_values_array1, positions=range(1, len(dates)))
    ax.violinplot(jsd_values_array2, positions=range(1, len(dates)))

    ax.set_xticks(range(1, len(dates)))
    ax.set_xticklabels([f"{dateRange[0].year}-{dateRange[0].month:02d}" for dateRange in zip(dates[:-1], dates[1:])], rotation=45)
    ax.set_title(gene)
    ax.set_xlabel("Collection Date")
    ax.set_ylabel("JSD between Influenza Gene and Human Infected with Influenza Codon Bias")
    fig.savefig("InfluenzaMedia/GeneJSD_Violin_2000_Onward_ExpressionComparison/" + gene + "_JSD_over_time.png")


for gene in InfluenzaGenes:
    dates = pl.date_range(start=pl.date(2008, 1, 1), end=pl.date(2010, 4, 1), interval="1mo", eager=True)
    fig = plt.figure(figsize=(40, 10))
    ax = fig.add_subplot(111)
    jsd_values_array = []
    for i, dateRange in enumerate(zip(dates[:-1], dates[1:])):
        subset2 = influenzaDFNewWhole.filter(
            (pl.col("geneID") == gene) &
            (pl.col("Collection_Date") >= dateRange[0]) &
            (pl.col("Collection_Date") < dateRange[1])
        )
        if len(subset2) == 0:
            jsd_values_array.append(np.array([0,0,0]))
            continue
        jsd_values = JSD(subset2.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_expressionWeighted[:, np.newaxis]) #type:ignore
        jsd_values_array.append(jsd_values)
    ax.violinplot(jsd_values_array, positions=range(1, len(dates)))
    ax.set_xticks(range(1, len(dates)))
    ax.set_xticklabels([f"{dateRange[0].year}-{dateRange[0].month:02d}" for dateRange in zip(dates[:-1], dates[1:])], rotation=45)
    ax.set_title(gene)
    ax.set_xlabel("Collection Date")
    ax.set_ylabel("JSD between Influenza Gene and Human Infected with Influenza Codon Bias")
    fig.savefig("InfluenzaMedia/GeneJSD_Violin_2008-2010/" + gene + "_JSD_over_time.png")

for gene in InfluenzaGenes:
    dates = pl.date_range(start=pl.date(2008, 1, 1), end=pl.date(2010, 4, 1), interval="1mo", eager=True)
    fig = plt.figure(figsize=(40, 10))
    ax = fig.add_subplot(111)
    jsd_values_array1 = []
    jsd_values_array2 = []
    for i, dateRange in enumerate(zip(dates[:-1], dates[1:])):
        subset2 = influenzaDFNewWhole.filter(
            (pl.col("geneID") == gene) &
            (pl.col("Collection_Date") >= dateRange[0]) &
            (pl.col("Collection_Date") < dateRange[1])
        )
        if len(subset2) == 0:
            jsd_values_array1.append(np.array([0,0,0]))
            jsd_values_array2.append(np.array([0,0,0]))
            continue
        jsd_values1 = JSD(subset2.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_expressionWeighted[:, np.newaxis]) #type:ignore
        jsd_values2 = JSD(subset2.select(CODON_ORDER).to_numpy().T, healthyHumanCodonBias_expressionWeighted[:, np.newaxis]) #type:ignore
        jsd_values_array1.append(jsd_values1)
        jsd_values_array2.append(jsd_values2)
    ax.violinplot(jsd_values_array1, positions=range(1, len(dates)))
    ax.violinplot(jsd_values_array2, positions=range(1, len(dates)))
    ax.set_xticks(range(1, len(dates)))
    ax.set_xticklabels([f"{dateRange[0].year}-{dateRange[0].month:02d}" for dateRange in zip(dates[:-1], dates[1:])], rotation=45)
    ax.set_title(gene)
    ax.set_xlabel("Collection Date")
    ax.set_ylabel("JSD between Influenza Gene and Human Infected with Influenza Codon Bias")
    fig.savefig("InfluenzaMedia/GeneJSD_Violin_2008-2010_ExpressionComparison/" + gene + "_JSD_over_time.png")

for gene in InfluenzaGenes:
    fig = plt.figure(figsize=(36*4, 36))
    axes = [fig.add_subplot(8, 8, i) for i in range(1, 65)]
    dates = pl.date_range(start=pl.date(2000, 1, 1), end=pl.date(2026, 4, 1), interval="3mo", eager=True)
    subset = influenzaDFNewWhole.filter(pl.col("geneID") == gene).filter(pl.col("Collection_Date") > pl.date(2000, 1, 1))
    for j, codon in enumerate(CODON_ORDER):
        stuff = []
        for i, dateRange in enumerate(zip(dates[:-1], dates[1:])):
            subset2 = subset.filter(
                (pl.col("Collection_Date") >= dateRange[0]) &
                (pl.col("Collection_Date") < dateRange[1])
            )
            if len(subset2) == 0:
                stuff.append(np.array([0,0,0]))
                continue
            stuff.append(subset2.select(codon).to_numpy().flatten().T)
        # print("hi")
        parts = axes[j].violinplot(stuff, positions=range(1, len(dates)))
        # print("hihi")
        for pc in parts['bodies']: #type: ignore
            pc.set_facecolor('#D43F3A')
            pc.set_edgecolor('black')
            pc.set_alpha(1)
        axes[j].set_title(codon)
        axes[j].legend()
        axes[j].set_xlabel("Collection Date")
        axes[j].set_ylabel("Codon Fraction")
        axes[j].set_xticks(range(1, len(dates)))
        axes[j].set_xticklabels([f"{dateRange[0].year}-{dateRange[0].month:02d}" for dateRange in zip(dates[:-1], dates[1:])], rotation=45)
    fig.savefig("InfluenzaMedia/GeneCodonJSD_2000_Onward/" + gene + "_JSD_over_time.png", dpi=300)

for gene in InfluenzaGenes:
    subset2 = influenzaDFNewWhole.filter(pl.col("geneID") == gene).filter(pl.col("Collection_Date") > pl.date(2008, 1, 1)).filter(pl.col("Collection_Date") < pl.date(2010, 4, 1))
    fig = plt.figure(figsize=(24, 24))
    axes = [fig.add_subplot(8, 8, i) for i in range(1, 65)]
    for j, codon in enumerate(CODON_ORDER):
        axes[j].scatter(np.array(subset2["Collection_Date"], dtype="datetime64"), subset2[codon], s=0.5, label=f"{codon} fraction")
        axes[j].set_title(codon)
        axes[j].legend()
        axes[j].set_xlabel("Collection Date")
        axes[j].set_ylabel("Codon Fraction")
    fig.savefig("InfluenzaMedia/GeneCodonJSD_2008-2010/" + gene + "_JSD_over_time.png", dpi=300)

for gene in InfluenzaGenes:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(range(len(CODON_ORDER)), influenzaDFNewWhole.filter(pl.col("geneID") == gene).select(CODON_ORDER).std().to_numpy()[0])
    ax.set_xticks(range(len(CODON_ORDER)))
    ax.set_xticklabels(CODON_ORDER, rotation=45)
    ax.set_title(f"Standard Deviation of Codon Fractions for {gene} Gene")
    ax.set_xlabel("Codon")
    ax.set_ylabel("Standard Deviation")
    fig.savefig(f"InfluenzaMedia/GeneCodonJSD_2000_Onward/{gene}_Codon_Fraction_Std.png", dpi=300)

geneCuttoffSTDs0 = {
    "M2": 0.1,
    "NA": 0.08,
    "PA": 0.04,
    "PB2": 0.1,
    "PB1": 0.1,
    "NEP": 0.1,
    "HA": 0.1,
    "NS1": 0.08,
    "M1": 0.1,
    "NP": 0.1
}

for gene in InfluenzaGenes:
    codons = influenzaDFNewWhole.filter(pl.col("geneID") == gene).filter(pl.col("Collection_Date") > pl.date(2000, 1, 1)).select(CODON_ORDER).std().to_numpy()[0]
    cutoff = geneCuttoffSTDs0[gene]
    codonsToInclude = [CODON_ORDER[i] for i, std in enumerate(codons) if std > cutoff]
    humanBias = humanWithCovidCodonBias_expressionWeighted[:, np.newaxis] #type:ignore
    # pick out codons from humanBias that are in codonsToInclude
    humanBiasSubset = np.array([humanBias[i] for i, codon in enumerate(CODON_ORDER) if codon in codonsToInclude])
    subset2 = influenzaDFNewWhole.filter(pl.col("geneID") == gene).filter(pl.col("Collection_Date") > pl.date(2000, 1, 1))
    jsd_values = JSD(subset2.select(codonsToInclude).to_numpy().T, humanBiasSubset) 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(np.array(subset2["Collection_Date"], dtype="datetime64"), jsd_values, s=0.5)
    ax.set_title(gene)
    ax.set_xlabel("Collection Date")
    ax.set_ylabel("JSD between Influenza Gene and Human Infected with Codon Bias")
    fig.savefig("InfluenzaMedia/GeneJSD_CodonSpecific_highBar_2000_onward/" + gene + "_JSD_over_time.png")

geneCuttoffSTDs1 = {
    "M2": 0.05,
    "NA": 0.05,
    "PA": 0.03,
    "PB2": 0.04,
    "PB1": 0.06,
    "NEP": 0.05,
    "HA": 0.1,
    "NS1": 0.05,
    "M1": 0.06,
    "NP": 0.04
}

for gene in InfluenzaGenes:
    codons = influenzaDFNewWhole.filter(pl.col("geneID") == gene).filter(pl.col("Collection_Date") > pl.date(2000, 1, 1)).select(CODON_ORDER).std().to_numpy()[0]
    cutoff = geneCuttoffSTDs1[gene]
    codonsToInclude = [CODON_ORDER[i] for i, std in enumerate(codons) if std > cutoff]
    humanBias = humanWithCovidCodonBias_expressionWeighted[:, np.newaxis] #type:ignore
    # pick out codons from humanBias that are in codonsToInclude
    humanBiasSubset = np.array([humanBias[i] for i, codon in enumerate(CODON_ORDER) if codon in codonsToInclude])
    subset2 = influenzaDFNewWhole.filter(pl.col("geneID") == gene).filter(pl.col("Collection_Date") > pl.date(2000, 1, 1))
    jsd_values = JSD(subset2.select(codonsToInclude).to_numpy().T, humanBiasSubset) 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(np.array(subset2["Collection_Date"], dtype="datetime64"), jsd_values, s=0.5)
    ax.set_title(gene)
    ax.set_xlabel("Collection Date")
    ax.set_ylabel("JSD between Influenza Gene and Human Infected with Codon Bias")
    fig.savefig("InfluenzaMedia/GeneJSD_CodonSpecific_medBar_2000_onward/" + gene + "_JSD_over_time.png")


geneCuttoffSTDs2 = {
    "M2": 0.02,
    "NA": 0.02,
    "PA": 0.01,
    "PB2": 0.02,
    "PB1": 0.02,
    "NEP": 0.025,
    "HA": 0.02,
    "NS1": 0.02,
    "M1": 0.03,
    "NP": 0.02
}


for gene in InfluenzaGenes:
    codons = influenzaDFNewWhole.filter(pl.col("geneID") == gene).filter(pl.col("Collection_Date") > pl.date(2000, 1, 1)).select(CODON_ORDER).std().to_numpy()[0]
    cutoff = geneCuttoffSTDs2[gene]
    codonsToInclude = [CODON_ORDER[i] for i, std in enumerate(codons) if std > cutoff]
    humanBias = humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis] #type:ignore
    # pick out codons from humanBias that are in codonsToInclude
    humanBiasSubset = np.array([humanBias[i] for i, codon in enumerate(CODON_ORDER) if codon in codonsToInclude])
    subset2 = influenzaDFNewWhole.filter(pl.col("geneID") == gene).filter(pl.col("Collection_Date") > pl.date(2000, 1, 1))
    jsd_values = JSD(subset2.select(codonsToInclude).to_numpy().T, humanBiasSubset) 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(np.array(subset2["Collection_Date"], dtype="datetime64"), jsd_values, s=0.5)
    ax.set_title(gene)
    ax.set_xlabel("Collection Date")
    ax.set_ylabel("JSD between Influenza Gene and Human Infected with Codon Bias")
    fig.savefig("InfluenzaMedia/GeneJSD_CodonSpecific_lowBar_2000_onward/" + gene + "_JSD_over_time.png")

for gene in InfluenzaGenes:
    dates = pl.date_range(start=pl.date(2000, 1, 1), end=pl.date(2026, 4, 1), interval="3mo", eager=True)
    fig = plt.figure(figsize=(24, 9))
    ax = fig.add_subplot(111)
    jsd_values_array1 = []
    jsd_values_array2 = []
    for i, dateRange in enumerate(zip(dates[:-1], dates[1:])):
        subset2 = H1N1.filter(
            (pl.col("geneID") == gene) &
            (pl.col("Collection_Date") >= dateRange[0]) &
            (pl.col("Collection_Date") < dateRange[1])
        )
        if len(subset2) == 0:
            jsd_values_array1.append(np.array([0,0,0]))
            jsd_values_array2.append(np.array([0,0,0]))
            continue
        jsd_values1 = JSD(subset2.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_expressionWeighted[:, np.newaxis]) #type:ignore
        jsd_values2 = JSD(subset2.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted[:, np.newaxis]) #type:ignore
        jsd_values_array1.append(jsd_values1)
        jsd_values_array2.append(jsd_values2)
        
    p1 = ax.violinplot(jsd_values_array1, positions=range(1, len(dates)))
    p2 = ax.violinplot(jsd_values_array2, positions=range(1, len(dates)))
    for pc in p1['bodies']: #type: ignore
            pc.set_facecolor("#0425FF")
            pc.set_edgecolor('black')
            pc.set_alpha(1)
    for pc in p2['bodies']: #type: ignore
            pc.set_facecolor("#E25700")
            pc.set_edgecolor('black')
            pc.set_alpha(1)
    step = 3
    ax.set_xticks(range(1, len(dates)))
    labels = [f"{dateRange[0].year}-{dateRange[0].month:02d}" for dateRange in zip(dates[:-1], dates[1:])]
    for i in range(0, len(labels)):
        if i % step != 0:
            labels[i] = ""
    ax.set_xticklabels(labels, rotation=45)
    ax.set_title(gene)
    ax.set_xlabel("Collection Date")
    ax.set_ylabel("JSD between Influenza Gene and Human Infected with Influenza Codon Bias")
    ax.legend([mpatches.Patch(color="#0425FF"), mpatches.Patch(color="#E25700")], ["Human with Covid Codon Bias No TPM Adjustment", "Human with Covid Codon Bias TPM Adjustment"], loc=2)
    fig.savefig("InfluenzaMedia/H1N1/TPMVis/" + gene + "_JSD_over_time.png")


for gene in InfluenzaGenes:
    dates = pl.date_range(start=pl.date(2000, 1, 1), end=pl.date(2026, 4, 1), interval="3mo", eager=True)
    fig = plt.figure(figsize=(24, 9))
    ax = fig.add_subplot(111)
    jsd_values_array1 = []
    jsd_values_array2 = []
    for i, dateRange in enumerate(zip(dates[:-1], dates[1:])):
        subset2 = H1N1.filter(
            (pl.col("geneID") == gene) &
            (pl.col("Collection_Date") >= dateRange[0]) &
            (pl.col("Collection_Date") < dateRange[1])
        )
        if len(subset2) == 0:
            jsd_values_array1.append(np.array([0,0,0]))
            jsd_values_array2.append(np.array([0,0,0]))
            continue
        jsd_values1 = JSD(subset2.select(CODON_ORDER).to_numpy().T, humanWithCovidCodonBias_expressionWeighted[:, np.newaxis]) #type:ignore
        jsd_values2 = JSD(subset2.select(CODON_ORDER).to_numpy().T, healthyHumanCodonBias_expressionWeighted[:, np.newaxis]) #type:ignore
        jsd_values_array1.append(jsd_values1)
        jsd_values_array2.append(jsd_values2)
        
    p1 = ax.violinplot(jsd_values_array1, positions=range(1, len(dates)))
    p2 = ax.violinplot(jsd_values_array2, positions=range(1, len(dates)))
    for pc in p1['bodies']: #type: ignore
            pc.set_facecolor("#0425FF")
            pc.set_edgecolor('black')
            pc.set_alpha(1)
    for pc in p2['bodies']: #type: ignore
            pc.set_facecolor("#E25700")
            pc.set_edgecolor('black')
            pc.set_alpha(1)
    step = 3
    ax.set_xticks(range(1, len(dates)))
    labels = [f"{dateRange[0].year}-{dateRange[0].month:02d}" for dateRange in zip(dates[:-1], dates[1:])]
    for i in range(0, len(labels)):
        if i % step != 0:
            labels[i] = ""
    ax.set_xticklabels(labels, rotation=45)
    ax.set_title(gene)
    ax.set_xlabel("Collection Date")
    ax.set_ylabel("JSD between Influenza Gene and Human Infected with Influenza Codon Bias")
    ax.legend([mpatches.Patch(color="#0425FF"), mpatches.Patch(color="#E25700")], ["Human with Covid Codon Bias", "Healthy Human Codon Bias"], loc=2)
    fig.savefig("InfluenzaMedia/H1N1/ExpressionVis/" + gene + "_JSD_over_time.png")