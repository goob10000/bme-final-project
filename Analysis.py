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
from numpy.typing import NDArray

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

covidCodonBiasMatrix = np.load("Results2/Covid/codon_proportions.npy")
covidGenomeDF = pl.read_csv("Genomes/Covid/MT291829_cds.tsv", separator="\t")
covidGenes = covidGenomeDF["gene_id"].to_list()
covidMutations = pl.read_parquet("Genomes/Covid/mutations.parquet")
covidMutations = covidMutations.filter(pl.col("mutations") != covidMutations["mutations"][1]) # Remove null (no mutations?)
covidDates = covidMutations["weeks"]
# (dates, codons, genes) = covidCodonBiasMatrix.shape

# Codon biases for all strains. Averaged across all genes in covid.
covid_small = covidCodonBiasMatrix.mean(axis=2) # Average across genes to get overall codon proportions at each timepoint

dfsWithAdjustmentTPM = adjustTPM(dfs, species, expression)
dfsWithAdjustmentAndMT = adjustTPM(dfs_withMT, species, expression)

species

humanCovid = 4
humanHealthy = 5
humanInfluenza = 6


dfs_expressionAdjusted = [weightByAdjustedExpression(df) for df in dfsWithAdjustmentTPM]
dfs_expressionAdjustedWithMT = [weightByAdjustedExpression(df) for df in dfsWithAdjustmentAndMT]
dfs_expressionUnadjusted = [weightByExpression(df) for df in dfsWithAdjustmentTPM]
dfs_expressionUnadjustedWithMT = [weightByExpression(df) for df in dfsWithAdjustmentAndMT]

## Define TPM adjusted expresison weighted codon biases
humanWithCovidCodonBias_TPMAdjusted_expressionWeighted = dfs_expressionAdjusted[humanCovid]
healthyHumanCodonBias_TPMAdjusted_expressionWeighted = dfs_expressionAdjusted[humanHealthy]
humanWithInfluenzaCodonBias_TPMAdjusted_expressionWeighted = dfs_expressionAdjusted[humanInfluenza]

##Define expression weighted codon biases (no TPM adjustment)
humanWithCovidCodonBias_expressionWeighted = dfs_expressionUnadjusted[humanCovid]
healthyHumanCodonBias_expressionWeighted = dfs_expressionUnadjusted[humanHealthy]
humanWithInfluenzaCodonBias_expressionWeighted = dfs_expressionUnadjusted[humanInfluenza]

n = 200
humanWithCovidCodonBias_expressionOrderedTop100 = dfsWithAdjustmentTPM[humanCovid].sort("adjustedExpression", descending=True).head(n).select(CODON_ORDER).mean()
healthyHumanCodonBias_expressionOrderedTop100 = dfsWithAdjustmentTPM[humanHealthy].sort("adjustedExpression", descending=True).head(n).select(CODON_ORDER).mean()
humanWithInfluenzaCodonBias_expressionOrderedTop100 = dfsWithAdjustmentTPM[humanInfluenza].sort("adjustedExpression", descending=True).head(n).select(CODON_ORDER).mean()

divergencesBetweenCovidAndHumanWithCovid = compare(covid_small, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted, JSD)
divergencesBetweenCovidAndHealthyHuman = compare(covid_small, healthyHumanCodonBias_TPMAdjusted_expressionWeighted, JSD)
divergencesBetweenCovidAndHumanWithInfluenza = compare(covid_small, humanWithInfluenzaCodonBias_TPMAdjusted_expressionWeighted, JSD)

influenza_matrix_path = "Results2/Influenza/codon_proportions.npy"
influenza_meta_path = "Results2/Influenza/codon_proportions_metadata.csv"

if not(os.path.exists(influenza_matrix_path) and os.path.exists(influenza_meta_path)):
    raise ValueError("Influenza data not found. Please ensure the codon proportions matrix and metadata CSV exist at the specified paths.")

influenzaCodonBiasMatrix:NDArray[np.float64] = np.load(influenza_matrix_path)
influenzaMeta:pl.DataFrame = pl.read_csv(influenza_meta_path)
influenzaCodonSmall:NDArray[np.float64] = influenzaCodonBiasMatrix.mean(axis=1)  # (rows, 62)
influenzaGenotypes:list[str] = influenzaMeta["Genotype"].to_list()
influenzaHosts:list[str] = influenzaMeta["Host"].to_list()
influenza_dates:list[datetime] = [datetime.fromordinal(int(x)) for x in influenzaMeta["date_ordinal"].to_list()]
divergencesBetweenInfluenzaAndHumanWithCovid:NDArray[np.float64] = compare(influenzaCodonSmall, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted, JSD)
divergencesBetweenInfluenzaAndHealthyHuman:NDArray[np.float64] = compare(influenzaCodonSmall, healthyHumanCodonBias_TPMAdjusted_expressionWeighted, JSD)
divergencesBetweenInfluenzaAndHumanWithInfluenza:NDArray[np.float64] = compare(influenzaCodonSmall, humanWithInfluenzaCodonBias_TPMAdjusted_expressionWeighted, JSD)

hosts_list: list[str] = []
influenzaColorsHosts = {}
influenzaCalculatedColorsHosts = []
host_groups = [[], [], [], []]

hosts_list = [str(h).strip() if h is not None else "" for h in influenzaHosts]
hosts_list = [h if h and h != "None" else "Unknown" for h in hosts_list]
unique_hosts = sorted(set(hosts_list))
host_palette = tab20b_cmap(np.linspace(0, 1, 20))
influenzaColorsHosts = {host: host_palette[i % len(host_palette)] for i, host in enumerate(unique_hosts)}
influenzaCalculatedColorsHosts = [influenzaColorsHosts[host] for host in hosts_list]

host_counts = {host: hosts_list.count(host) for host in unique_hosts}
hosts_by_frequency = sorted(unique_hosts, key=lambda h: host_counts[h], reverse=True)
host_groups = [[], [], [], []]
for i, host in enumerate(hosts_by_frequency):
    host_groups[i % 4].append(host)

if not(influenza_dates is not None and divergencesBetweenInfluenzaAndHumanWithCovid is not None and influenzaMeta is not None):
    raise ValueError("Influenza data not found. Please ensure the codon proportions matrix and metadata CSV exist at the specified paths.")

dfsWithAdjustmentTPM[humanCovid].sort("adjustedExpression", descending=True).head(5)["gene_id"]
list(zip(humanWithCovidCodonBias_expressionOrderedTop100.columns, humanWithCovidCodonBias_expressionOrderedTop100.to_numpy().tolist()[0]))[0]
d = {x: y for x, y in zip(humanWithCovidCodonBias_expressionOrderedTop100.columns, humanWithCovidCodonBias_expressionOrderedTop100.to_numpy().tolist()[0])}
d

fig1 = plt.figure(layout="constrained", figsize=(10, 8))
ax1 = fig1.add_subplot(211)
ax2 = fig1.add_subplot(212)
ax1.scatter(covidDates, divergencesBetweenCovidAndHumanWithCovid, s=0.5, label="SARS-CoV-2 Codon Bias Divergence compared to Human infected with SARS-CoV-2")
ax1.scatter(covidDates, divergencesBetweenCovidAndHealthyHuman, s=0.5, label="SARS-CoV-2 Codon Bias Divergence compared to Healthy Human")
ax1.scatter(covidDates, divergencesBetweenCovidAndHumanWithInfluenza, s=0.5, label="SARS-CoV-2 Codon Bias Divergence compared to Human infected with Influenza")
ax1.legend()
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.set_xlabel("Weeks since Week 1 of 2020")
ax1.set_ylabel("Jensen-Shannon Divergence")
ax1.set_title("Divergence of SARS-CoV-2 Codon Bias from Human Codon Bias Over Time")
ax2.scatter(covidDates,compare(covid_small, humanWithCovidCodonBias_expressionOrderedTop100.to_numpy()[0,:], JSD), s=0.5, label="SARS-CoV-2 Codon Bias Divergence to Covid Infected Individual (Top 100 Genes)")
ax2.scatter(covidDates,compare(covid_small, healthyHumanCodonBias_expressionOrderedTop100.to_numpy()[0, :], JSD), s=0.5, label="SARS-CoV-2 Codon Bias Divergence to Healthy Individual (Top 100 Genes, TPM adjustment)")
ax2.scatter(covidDates,compare(covid_small, humanWithInfluenzaCodonBias_expressionOrderedTop100.to_numpy()[0, :], JSD), s=0.5, label="SARS-CoV-2 Codon Bias Divergence to Influenza Infected Individual (Top 100 Genes, TPM adjustment)")
ax2.legend()
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.set_xlabel("Weeks since Week 1 of 2020")
ax2.set_ylabel("Jensen-Shannon Divergence")
ax2.set_title("Divergence of SARS-CoV-2 Codon Bias from Human Codon Bias Over Time (Top 100 Genes)")
fig1.show()

fig1_1 = plt.figure(layout="constrained", figsize=(10, 2.5))
ax1 = fig1_1.add_subplot(111)
ax1.scatter(covidDates, divergencesBetweenCovidAndHumanWithCovid, s=0.5, label="SARS-CoV-2 Codon Bias Divergence compared to Human infected with SARS-CoV-2")
ax1.scatter(covidDates, divergencesBetweenCovidAndHealthyHuman, s=0.5, label="SARS-CoV-2 Codon Bias Divergence compared to Healthy Human")
ax1.scatter(covidDates, divergencesBetweenCovidAndHumanWithInfluenza, s=0.5, label="SARS-CoV-2 Codon Bias Divergence compared to Human infected with Influenza")
ax1.legend()
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.set_xlabel("Weeks since Week 1 of 2020")
ax1.set_ylabel("Jensen-Shannon Divergence")
ax1.set_title("Divergence of SARS-CoV-2 Codon Bias from Human Codon Bias Over Time")
fig1_1.savefig(dpi=600, fname="covidDivergenceOverTime.png")



fig2 = plt.figure()
ax3 = fig2.add_subplot(211)
ax3.scatter(influenza_dates, divergencesBetweenInfluenzaAndHumanWithCovid,s=2, label="Influenza H1N1 Codon Bias Divergence in an infected individual") #type: ignore
ax3.set_xlabel("Date")
ax3.set_ylabel("Jensen-Shannon Divergence")
ax3.legend()

ax4 = fig2.add_subplot(212)
ax4.scatter(influenza_dates, divergencesBetweenInfluenzaAndHealthyHuman, s=2, color="orange", label="Influenza H1N1 Codon Bias Divergence in a healthy individual") #type: ignore
ax4.set_xlabel("Date")
ax4.set_ylabel("Jensen-Shannon Divergence")
ax4.legend()
fig2.show()

fig2_1 = plt.figure(layout="constrained", figsize=(10, 2.5))
ax3 = fig2_1.add_subplot(111)
ax3.scatter(influenza_dates, divergencesBetweenInfluenzaAndHumanWithCovid,s=2, label="Influenza H1N1 Codon Bias Divergence in an infected individual") #type: ignore
ax3.set_xlabel("Date")
ax3.set_ylabel("Jensen-Shannon Divergence")
ax3.set_xlim((16000, 20000))
ax3.set_ylim((0.035, 0.07))
ax3.legend()
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
# fig2_1.show()
fig2_1.savefig(dpi=600, fname="influenzaDivergenceInfected.png")

fig3 = plt.figure()
ax5 = fig3.add_subplot(211)
ax5.scatter(covidDates, compare(covid_small, humanWithCovidCodonBias_expressionWeighted, diff), s=0.5, label="SARS-CoV-2 Codon Bias Divergence in an infected individual (no TPM adjustment)")
ax5.set_xlabel("Weeks since Week 1 of 2020")
ax5.set_ylabel("Difference in Codon Proportions")
ax5.legend()
ax6 = fig3.add_subplot(212)
ax6.scatter(covidDates, compare(covid_small, humanWithCovidCodonBias_TPMAdjusted_expressionWeighted, diff), s=0.5, label="SARS-CoV-2 Codon Bias Divergence in an infected individual (WITH TPM adjustment)", color="orange")
ax6.set_xlabel("Weeks since Week 1 of 2020")
ax6.set_ylabel("Difference in Codon Proportions")
ax6.legend()
fig3.show()

### for each species plot the divergence over time
fig4 = plt.figure()
ax7 = fig4.add_subplot(111)
for i, s in enumerate(species):
    divergence = compare(covid_small, weightByAdjustedExpression(dfsWithAdjustmentTPM[i]), diff)
    ax7.scatter(covidDates, divergence, s=0.5, label=s)
ax7.set_xlabel("Weeks since Week 1 of 2020")
ax7.set_ylabel("Difference in Codon Proportions")
ax7.legend()
fig4.show()

### for each species plot the divergence of the top 10 genes and the top 100 genes over time
fig5 = plt.figure()
ax8 = fig5.add_subplot(211)
for i, s in enumerate(species):
    top10Genes = dfsWithAdjustmentTPM[i].sort("adjustedExpression", descending=True).head(10)["gene_id"].to_list()
    top10Divergence = compare(covid_small, weightByAdjustedExpression(dfsWithAdjustmentTPM[i].filter(pl.col("gene_id").is_in(top10Genes))), diff)
    ax8.scatter(covidDates, top10Divergence, s=0.5, label=s)
ax8.set_xlabel("Weeks since Week 1 of 2020")
ax8.set_ylabel("Difference in Codon Proportions (Top 10 Genes)")
ax8.legend()
ax9 = fig5.add_subplot(212)
for i, s in enumerate(species):
    top100Genes = dfsWithAdjustmentTPM[i].sort("adjustedExpression", descending=True).head(100)["gene_id"].to_list()
    top100Divergence = compare(covid_small, weightByAdjustedExpression(dfsWithAdjustmentTPM[i].filter(pl.col("gene_id").is_in(top100Genes))), diff)
    ax9.scatter(covidDates, top100Divergence, s=0.5, label=s)
ax9.set_xlabel("Weeks since Week 1 of 2020")
ax9.set_ylabel("Difference in Codon Proportions (Top 100 Genes)")
ax9.legend()
fig5.show()

fig6 = plt.figure()

if influenza_dates is not None and divergencesBetweenInfluenzaAndHumanWithCovid is not None and divergencesBetweenInfluenzaAndHealthyHuman is not None:
    ax10 = fig6.add_subplot(211)
    ax10.scatter(influenza_dates, divergencesBetweenInfluenzaAndHumanWithCovid, s=0.5, label="Influenza H1N1 Codon Bias Divergence in an infected individual") #type: ignore
    ax10.set_xlabel("Date")
    ax10.set_ylabel("Jensen-Shannon Divergence")
    ax10.legend()

    ax11 = fig6.add_subplot(212)
    ax11.scatter(influenza_dates, divergencesBetweenInfluenzaAndHealthyHuman, s=0.5, label="Influenza H1N1 Codon Bias Divergence in a healthy individual", color="orange") #type: ignore
    ax11.set_xlabel("Date")
    ax11.set_ylabel("Jensen-Shannon Divergence")
    ax11.legend()

    fig6.show()

covid1 = covid_small[0,:]

arrs = []
divs = []

for i, sp in enumerate(species):
    print(f"{dfsWithAdjustmentTPM[i].shape=}")
    arr = np.zeros((500, 62), dtype=np.float64)
    for j in range(0, 4999, 10):
        arr[j//10] = dfsWithAdjustmentTPM[i].sort("adjustedExpression", descending=True).head(j+1).select(CODON_ORDER).mean() # Average codon bias of top j genes
    arrs.append(arr)

for i, sp in enumerate(species):
    divergence = compare(arr, covid1, JSD)
    divs.append(divergence)
    print(divergence.shape)
    print("inf") if np.isinf(divergence).any() else "no inf"

arrsSmall = []

for i, sp in enumerate(species):
    arr = np.zeros((10, 62), dtype=np.float64)
    for j in range(0, 20, 2):
        arr[j//2] = dfsWithAdjustmentTPM[i].sort("adjustedExpression", descending=True).head(j+1).select(CODON_ORDER).mean() # Average codon bias of top j genes
    arrsSmall.append(arr)

## Plot showing divergence by top-n genes compared to single covid sequence
fig7 = plt.figure()
ax12 = fig7.add_subplot(111)
for i in range(len(species)):
    ax12.plot(np.arange(0, 5000, 10), divs[i], label=species[i], c=species_colors[i])
ax12.set_yscale("log")
ax12.set_xlabel("Number of Top Genes Included")
ax12.set_ylabel("Jensen-Shannon Divergence to Single SARS-CoV-2 Sequence")
ax12.legend()
ax12.set_ylim(divergence.min()/2, divergence.max()*3)

fig7.show()

gly = [arrs[i][:,39:43] for i in range(len(species))]

fig8 = plt.figure()
ax16 = fig8.add_subplot(221)
ax17 = fig8.add_subplot(222)
ax18 = fig8.add_subplot(223)
ax19 = fig8.add_subplot(224)

for i in range(1, len(species)):
    ax16.plot(np.arange(0, 5000, 10), gly[i][:,0], label=species[i], c=species_colors[i])
    ax17.plot(np.arange(0, 5000, 10), gly[i][:,1], label=species[i], c=species_colors[i])
    ax18.plot(np.arange(0, 5000, 10), gly[i][:,2], label=species[i], c=species_colors[i])
    ax19.plot(np.arange(0, 5000, 10), gly[i][:,3], label=species[i], c=species_colors[i])
ax16.set_xlabel("Number of Top Genes Included")
ax16.set_ylabel("Proportion of GGT Codon (Glycine)")
ax17.set_xlabel("Number of Top Genes Included")
ax17.set_ylabel("Proportion of GGC Codon (Glycine)")
ax18.set_xlabel("Number of Top Genes Included")
ax18.set_ylabel("Proportion of GGA Codon (Glycine)")
ax19.set_xlabel("Number of Top Genes Included")
ax19.set_ylabel("Proportion of GGG Codon (Glycine)")

ax16.set_xscale("log")
ax17.set_xscale("log")
ax18.set_xscale("log")
ax19.set_xscale("log")
# ax16.legend()
# ax17.legend()
# ax18.legend()
# ax19.legend()
fig8.show()

divs[1]

ax13 = fig7.add_subplot(222)
ax14 = fig7.add_subplot(223)
ax15 = fig7.add_subplot(224)



fig9 = plt.figure()

if influenza_dates is not None and divergencesBetweenInfluenzaAndHumanWithCovid is not None and divergencesBetweenInfluenzaAndHealthyHuman is not None:
    influenza_dates_arr = np.array(influenza_dates)
    influenza_hCC_arr = np.array(divergencesBetweenInfluenzaAndHumanWithCovid)
    hosts_arr = np.array(hosts_list)

    axes = [
        fig9.add_subplot(221),
        fig9.add_subplot(222),
        fig9.add_subplot(223),
        fig9.add_subplot(224),
    ]

    for idx, ax in enumerate(axes):
        subset_hosts = host_groups[idx]
        for host in subset_hosts:
            mask = hosts_arr == host
            if not np.any(mask):
                continue
            ax.scatter(
                influenza_dates_arr[mask],
                influenza_hCC_arr[mask],
                s=0.5,
                color=influenzaColorsHosts[host],
                label=host,
            )
        ax.set_title(f"Host subset {idx + 1} ({len(subset_hosts)} hosts)")
        ax.set_xlabel("Date")
        ax.set_ylabel("Jensen-Shannon Divergence")
        ax.legend(fontsize=6, markerscale=4, frameon=False)

    fig9.tight_layout()

    fig9.show()


fig10 = plt.figure()


if influenza_dates is not None and divergencesBetweenInfluenzaAndHumanWithCovid is not None and influenzaMeta is not None:
    influenza_dates_arr = np.array(influenza_dates)
    influenza_hCC_arr = np.array(divergencesBetweenInfluenzaAndHumanWithCovid)
    hosts_arr = np.array(hosts_list)
    
    countries_list = [str(c).strip() if c is not None else "" for c in influenzaMeta["Country"].to_list()]
    countries_list = [c if c and c != "None" else "Unknown" for c in countries_list]
    
    # Filter for humans only
    human_mask = hosts_arr == "Homo sapiens"
    human_dates = influenza_dates_arr[human_mask]
    human_hCC = influenza_hCC_arr[human_mask]
    human_countries = [countries_list[i] for i in range(len(countries_list)) if human_mask[i]]
    
    if len(human_countries) > 0:
        unique_countries = sorted(set(human_countries))
        country_palette = tab20b_cmap(np.linspace(0, 1, 20))
        influenzaColorsCountries = {country: country_palette[i % len(country_palette)] for i, country in enumerate(unique_countries)}
        
        country_counts = {country: human_countries.count(country) for country in unique_countries}
        countries_by_frequency = sorted(unique_countries, key=lambda c: country_counts[c], reverse=True)
        
        country_groups = [[], [], [], []]
        for i, country in enumerate(countries_by_frequency):
            country_groups[i % 4].append(country)
        
        axes = [
            fig10.add_subplot(221),
            fig10.add_subplot(222),
            fig10.add_subplot(223),
            fig10.add_subplot(224),
        ]
        
        for idx, ax in enumerate(axes):
            subset_countries = country_groups[idx]
            for country in subset_countries:
                mask = np.array([c == country for c in human_countries])
                if not np.any(mask):
                    continue
                ax.scatter(
                    human_dates[mask],
                    human_hCC[mask],
                    s=0.5,
                    color=influenzaColorsCountries[country],
                    label=country,
                )
            ax.set_title(f"Country subset {idx + 1} ({len(subset_countries)} countries)")
            ax.set_xlabel("Date")
            ax.set_ylabel("Jensen-Shannon Divergence")
            ax.legend(fontsize=6, markerscale=4, frameon=False)
        
        fig10.tight_layout()
    
    fig10.show()



gly = arrs[6][:,39:43]
glySmall = arrsSmall[6][:,39:43]
gly_humanInfluenza = dfs_expressionAdjusted[humanInfluenza][39:43]
fig11 = plt.figure()
ax20 = fig11.add_subplot(221)
ax21 = fig11.add_subplot(222)
ax22 = fig11.add_subplot(223)
ax23 = fig11.add_subplot(224)

ax20.plot(np.arange(20, 5000, 10), gly[2:,0], label="Expression Sorted Bias")
ax20.plot(np.arange(0,20,2), glySmall[:,0], label="Expression Sorted Bias (Top 20 Genes)")
ax20.plot((0,5000), (gly_humanInfluenza[0], gly_humanInfluenza[0]), label="Expression weighted bias", color="black")
ax20.set_xlabel("Number of Top Genes Included")
ax20.set_ylabel("Proportion of GGT Codon (Glycine)")
ax20.set_xscale("log")
ax20.legend()

ax21.plot(np.arange(20, 5000, 10), gly[2:,1], label="Expression Sorted Bias")
ax21.plot(np.arange(0,20,2), glySmall[:,1], label="Expression Sorted Bias (Top 20 Genes)")
ax21.plot((0,5000), (gly_humanInfluenza[1], gly_humanInfluenza[1]), label="Expression weighted bias", color="black")
ax21.set_xlabel("Number of Top Genes Included")
ax21.set_ylabel("Proportion of GGC Codon (Glycine)")
ax21.set_xscale("log")
ax21.legend()

ax22.plot(np.arange(20, 5000, 10), gly[2:,2], label="Expression Sorted Bias")
ax22.plot(np.arange(0,20,2), glySmall[:,2], label="Expression Sorted Bias (Top 20 Genes)")
ax22.plot((0,5000), (gly_humanInfluenza[2], gly_humanInfluenza[2]), label="Expression weighted bias", color="black")
ax22.set_xlabel("Number of Top Genes Included")
ax22.set_ylabel("Proportion of GGA Codon (Glycine)")
ax22.set_xscale("log")
ax22.legend()

ax23.plot(np.arange(20, 5000, 10), gly[2:,3], label="Expression Sorted Bias")
ax23.plot(np.arange(0,20,2), glySmall[:,3], label="Expression Sorted Bias (Top 20 Genes)")
ax23.plot((0,5000), (gly_humanInfluenza[3], gly_humanInfluenza[3]), label="Expression weighted bias", color="black")
ax23.set_xlabel("Number of Top Genes Included")
ax23.set_ylabel("Proportion of GGG Codon (Glycine)")
ax23.set_xscale("log")
ax23.legend()

fig11.show()

fig12 = plt.figure(figsize=(24, 24), layout="constrained")

axes = [
    fig12.add_subplot(8,8,i+1) for i in range(62)
]

for i, ax in enumerate(axes):
    ax.plot(np.arange(20,5000,10), arrs[6][2:,i], label="Expression Sorted Bias")
    ax.plot(np.arange(0,20,2), arrsSmall[6][:,i], label="Expression Sorted Bias (Top 20 Genes)")
    ax.plot((0,5000), (dfs_expressionAdjusted[humanInfluenza][i], dfs_expressionAdjusted[humanInfluenza][i]), label="Expression weighted bias", color="black")
    ax.set_xscale("log")
    ax.set_title(CODON_ORDER[i])

fig12.savefig("codonBiasExpressionSortedWeightedComparison.png", dpi=300)

fig13 = plt.figure()

if influenza_dates is not None and divergencesBetweenInfluenzaAndHumanWithCovid is not None and influenzaMeta is not None:
    influenzaMetaMichigan = influenzaMeta.filter(pl.col("Geo_Location") == "USA: Michigan")
    influenzaMetaCalifornia = influenzaMeta.filter(pl.col("Geo_Location") == "USA: California")
    iMJapan = influenzaMeta.filter(pl.col("Geo_Location") == "Japan")
    iMUSA = influenzaMeta.filter(pl.col("Geo_Location") == "USA")
    influenza_dates_arr = np.array(influenza_dates)
    influenza_hCC_arr = np.array(divergencesBetweenInfluenzaAndHumanWithCovid)
    hosts_arr = np.array(hosts_list)
    
    influenzaRowsMichigan = influenzaMetaMichigan["row_index"].to_numpy()
    influenzaRowsCalifornia = influenzaMetaCalifornia["row_index"].to_numpy()
    iRJapan = iMJapan["row_index"].to_numpy()
    iRUSA = iMUSA["row_index"].to_numpy()

    a = np.arange(len(influenza_hCC_arr)) 
    mich = [False if i not in influenzaRowsMichigan else True for i in a]
    cal = [False if i not in influenzaRowsCalifornia else True for i in a]
    japan = [False if i not in iRJapan else True for i in a]
    usa = [False if i not in iRUSA else True for i in a]

    michNumpy = influenza_hCC_arr[mich]
    calNumpy = influenza_hCC_arr[cal]
    japanNumpy = influenza_hCC_arr[japan]
    usaNumpy = influenza_hCC_arr[usa]
    
    ax = fig13.add_subplot(111)
    ax.scatter(influenza_dates_arr[mich], influenza_hCC_arr[mich], label="Michigan", s=0.2)
    ax.scatter(influenza_dates_arr[cal], influenza_hCC_arr[cal], label="California", s=0.2)
    ax.scatter(influenza_dates_arr[japan], influenza_hCC_arr[japan], label="Japan", s=0.2)
    ax.scatter(influenza_dates_arr[usa], influenza_hCC_arr[usa], label="USA (unspecified location)", s=0.2)
    ax.set_xlabel("Date")
    ax.legend()
        
    fig13.tight_layout()
    fig13.show()

fig14 = plt.figure()
if influenza_dates is not None and divergencesBetweenInfluenzaAndHumanWithCovid is not None and influenzaMeta is not None:
    influenzaMetaMichigan = influenzaMeta.filter(pl.col("Geo_Location") == "USA: Michigan")
    influenzaMetaCalifornia = influenzaMeta.filter(pl.col("Geo_Location") == "USA: California")
    iMJapan = influenzaMeta.filter(pl.col("Geo_Location") == "Japan")
    iMUSA = influenzaMeta.filter(pl.col("Geo_Location") == "USA")
    influenza_dates_arr = np.array(influenza_dates)
    influenza_hCC_arr = np.array(divergencesBetweenInfluenzaAndHumanWithCovid)
    hosts_arr = np.array(hosts_list)
    
    influenzaRowsMichigan = influenzaMetaMichigan["row_index"].to_numpy()
    influenzaRowsCalifornia = influenzaMetaCalifornia["row_index"].to_numpy()
    iRJapan = iMJapan["row_index"].to_numpy()
    iRUSA = iMUSA["row_index"].to_numpy()

    a = np.arange(len(influenza_hCC_arr)) 
    mich = [False if i not in influenzaRowsMichigan else True for i in a]
    cal = [False if i not in influenzaRowsCalifornia else True for i in a]
    japan = [False if i not in iRJapan else True for i in a]
    usa = [False if i not in iRUSA else True for i in a]

    michNumpy = influenza_hCC_arr[mich]
    calNumpy = influenza_hCC_arr[cal]
    japanNumpy = influenza_hCC_arr[japan]
    usaNumpy = influenza_hCC_arr[usa]
    
    ax = fig14.add_subplot(111)
    d = {"Michigan": influenza_hCC_arr[mich], "California": influenza_hCC_arr[cal], "Japan": influenza_hCC_arr[japan], "USA (unspecified location)": influenza_hCC_arr[usa]}

    ax.violinplot(d.values()) #type: ignore
    ax.set_xlabel("Date")
    ax.legend(["Michigan", "California", "Japan", "USA (unspecified location)"])
        
    fig14.tight_layout()
    fig14.show()

fig16 = plt.figure()
ax1 = fig16.add_subplot(111)
ax1.scatter(np.array(influenza_dates), divergencesBetweenInfluenzaAndHumanWithCovid, s=0.5, label="Influenza A/H1N1 Codon Bias Divergence in an infected individual")
ax1.scatter(np.array(influenza_dates), divergencesBetweenInfluenzaAndHealthyHuman, s=0.5, label="Influenza A/H1N1 Codon Bias Divergence in a healthy individual")
ax1.scatter(np.array(influenza_dates), divergencesBetweenInfluenzaAndHumanWithInfluenza, s=0.5, label="Influenza A/H1N1 Codon Bias Divergence in an individual with Influenza")
ax1.legend()
ax1.set_xlabel("Year")
ax1.set_ylabel("Jensen-Shannon Divergence")
ax1.set_title("Divergence of Influenza A/H1N1 Codon Bias from Human Codon Bias Over Time")
fig16.show()

fig17 = plt.figure()
ax1 = fig17.add_subplot(111)
l = []
l.append(covid_small[0,:])
l.append(influenzaCodonSmall[0,:])
for i in controlVirusCBs:
    l.append(i)
for i, v in enumerate(viruses):
    print(f"{v=}")
    print(f"{l[i].shape=}")

d = [JSD(humanWithCovidCodonBias_TPMAdjusted_expressionWeighted, virusCB) for virusCB in l]
d2 = [JSD(healthyHumanCodonBias_TPMAdjusted_expressionWeighted, virusCB) for virusCB in l]
x = np.arange(len(viruses))
width = 0.4
ax1.bar(x - width/2, d2, width=width, alpha=0.5, label="Divergence from Healthy Human Codon Bias")
ax1.bar(x + width/2, d, width=width, alpha=0.9, label="Divergence from Covid Infected Human Codon Bias")
ax1.set_xticks(x)
ax1.set_xticklabels(viruses)
ax1.set_ylabel("Jensen-Shannon Divergence")
ax1.set_title("Divergence of Virus Codon Bias from Covid Infected Human Codon Bias")
ax1.legend()
fig17.show()


fig18 = plt.figure(layout="constrained")
ax1 = fig18.add_subplot(111)
arr = np.zeros((len(species), len(species)))

for i in range(len(species)):
    for j in range(len(species)):
        arr[i, j] = JSD(dfs_expressionAdjusted[i], dfs_expressionAdjusted[j])

im = ax1.imshow(arr, cmap="viridis")
ax1.set_xticks(np.arange(len(species)))
ax1.set_yticks(np.arange(len(species)))
ax1.set_xticklabels(species, rotation=90)
ax1.set_yticklabels(species)
ax1.set_title("Jensen-Shannon Divergence Between Species Based on Expression-Weighted Codon Bias")
fig18.colorbar(im, ax=ax1, label="Jensen-Shannon Divergence")
fig18.show()

fig19 = plt.figure(layout="constrained")
ax1 = fig19.add_subplot(111)
arrNoMicrobes = np.zeros((len(species)-2, len(species)-2))
arrNoMicrobes[0:3, 0:3] = arr[0:3, 0:3]
arrNoMicrobes[3:, 3:] = arr[4:-1, 4:-1]
arrNoMicrobes[3:, 0:3] = arr[4:-1, 0:3]
arrNoMicrobes[0:3, 3:] = arr[0:3, 4:-1]
speciesNoMicrobes = species[0:3] + species[4:-1]
im = ax1.imshow(arrNoMicrobes, cmap="viridis")
ax1.set_xticks(np.arange(len(speciesNoMicrobes)))
ax1.set_yticks(np.arange(len(speciesNoMicrobes)))
ax1.set_xticklabels(speciesNoMicrobes, rotation=90)
ax1.set_yticklabels(speciesNoMicrobes)
ax1.set_title("Jensen-Shannon Divergence Between Species Based on Expression-Weighted Codon Bias")
fig19.colorbar(im, ax=ax1, label="Jensen-Shannon Divergence")
fig19.show()

fig20 = plt.figure(layout="constrained")
ax1 = fig20.add_subplot(111)
species
speciesToRemove = ["Yeast_S288C", 'Pangolin_ManisPentadactyla', "Ecoli_K12_MG1655", "HumanHealthy1"]
indicesToRemove = [species.index(s) for s in speciesToRemove]
arrZeroedMicrobes = arr.copy()
for idx in indicesToRemove:
    arrZeroedMicrobes[idx, :] = 0
    arrZeroedMicrobes[:, idx] = 0
im = ax1.imshow(arrZeroedMicrobes, cmap="viridis")
ax1.set_xticks(np.arange(len(species)))
ax1.set_yticks(np.arange(len(species)))
ax1.set_xticklabels(species, rotation=90)
ax1.set_yticklabels(species)
ax1.set_title("Jensen-Shannon Divergence Between Species Based on Expression-Weighted Codon Bias")
fig20.colorbar(im, ax=ax1, label="Jensen-Shannon Divergence")
fig20.show()

from scipy.cluster.hierarchy import dendrogram, linkage

linkage_matrix = linkage(arr, method='ward')
fig21 = plt.figure(figsize=(10, 7))
dendrogram(linkage_matrix, labels=species, leaf_rotation=90)
plt.title("Hierarchical Clustering of Species Based on Expression-Weighted Codon Bias")
plt.xlabel("Species")
plt.ylabel("Distance")
plt.tight_layout()
fig21.show()


all_species = species + viruses
all_arr = np.zeros((len(all_species), len(all_species)))
for i in range(len(all_species)):
    for j in range(len(all_species)):
        if i < len(species) and j < len(species):
            all_arr[i, j] = JSD(dfs_expressionAdjusted[i], dfs_expressionAdjusted[j])
        elif i >= len(species) and j >= len(species):
            all_arr[i, j] = JSD(l[i - len(species)], l[j - len(species)])
        elif i < len(species) and j >= len(species):
            all_arr[i, j] = JSD(dfs_expressionAdjusted[i], l[j - len(species)])
        else:
            all_arr[i, j] = JSD(l[i - len(species)], dfs_expressionAdjusted[j])
linkage_matrix_all = linkage(all_arr, method='ward')
fig22 = plt.figure(figsize=(12, 8))
dendrogram(linkage_matrix_all, labels=all_species, leaf_rotation=90)
plt.title("Hierarchical Clustering of Species and Viruses Based on Expression-Weighted Codon Bias")
plt.xlabel("Species/Virus")
plt.ylabel("Distance")
# plt.yscale("log")
plt.tight_layout()
fig22.show()

fig23 = plt.figure(figsize=(12, 8))
im = plt.imshow(all_arr, cmap="viridis")
plt.xticks(np.arange(len(all_species)), all_species, rotation=90)
plt.yticks(np.arange(len(all_species)), all_species)
plt.title("Jensen-Shannon Divergence Between Species and Viruses Based on Expression-Weighted Codon Bias")
plt.colorbar(im, label="Jensen-Shannon Divergence")
plt.tight_layout()
fig23.show()


## Heat graph of divergence between all species but ordered by the hierarchical clustering of all the species and viruses together
fig24 = plt.figure(figsize=(12, 8))
# Get the order of species based on hierarchical clustering
dendro = dendrogram(linkage_matrix_all, labels=all_species, no_plot=True)
ordered_species = [all_species[i] for i in dendro['leaves']]
ordered_arr = np.zeros((len(all_species), len(all_species)))
for i in range(len(all_species)):
    for j in range(len(all_species)):
        idx_i = all_species.index(ordered_species[i])
        idx_j = all_species.index(ordered_species[j])
        ordered_arr[i, j] = all_arr[idx_i, idx_j]

im = plt.imshow(ordered_arr, cmap="viridis")
plt.xticks(np.arange(len(all_species)), ordered_species, rotation=90)
plt.yticks(np.arange(len(all_species)), ordered_species)
plt.title("Jensen-Shannon Divergence Between Species and Viruses Based on Expression-Weighted Codon Bias (Ordered by Clustering)")
plt.colorbar(im, label="Jensen-Shannon Divergence")
fig24.show()


# fig15 = plt.figure()

# biases = []
# biases.append(arrsSmall[6][-1])
# for i in range(24):
#     biases.append(arrs[6][i*2+2])
# biases = np.array(biases)
# divsTwentyInfluenza = compare(biases, , JSD)




## Surface plot of divergence over time and average codon bias of top n genes for each species
fig25 = plt.figure()
ax10 = fig25.add_subplot(111, projection='3d')
for i, s in enumerate(species):
    print(s)
    arr = np.zeros((500, 62), dtype=np.float64)
    for j in range(0, 4999, 10):
        arr[j//10] = dfsWithAdjustmentTPM[i].select(CODON_ORDER)[:j+1].mean() # Average codon bias of top j genes
    ## I want a numpy way to calculate Jensen-Shannon Divergence for each codon proportions in covidDates to every row in arr
    X, Y = np.meshgrid(covidDates, arr)
    print(X)
    divergence = compare(covid_small, arr, JSD) 


## Check slopes

a = np.polyfit(covidDates, divergencesBetweenCovidAndHumanWithCovid, deg=1)
b = np.polyfit(covidDates, divergencesBetweenCovidAndHealthyHuman, deg=1)
c = np.polyfit(covidDates, divergencesBetweenCovidAndHumanWithInfluenza, deg=1)
d = np.polyfit(covidDates, compare(covid_small, humanWithCovidCodonBias_expressionOrderedTop100.to_numpy()[0,:], JSD), deg=1)
e = np.polyfit(covidDates, compare(covid_small, healthyHumanCodonBias_expressionOrderedTop100.to_numpy()[0, :], JSD), deg=1)
f = np.polyfit(covidDates, compare(covid_small, humanWithInfluenzaCodonBias_expressionOrderedTop100.to_numpy()[0, :], JSD), deg=1)

a
b
c
d
e
f

### Tried to recreate CBI vs expression plot but didn't really work.
fig26 = plt.figure()
ax1 = fig26.add_subplot(111)
ax1.scatter(1 - JSD(humanWithCovidCodonBias_expressionOrderedTop100.to_numpy().T, dfsWithAdjustmentTPM[humanCovid].select(CODON_ORDER).to_numpy().T), dfsWithAdjustmentTPM[humanCovid]["expression"], s=0.5, label = "Similarity to SARS-CoV-2 Codon Bias (Top 100 Genes)")
ax1.scatter(1 - JSD(humanWithCovidCodonBias_expressionWeighted[:, np.newaxis], dfsWithAdjustmentTPM[humanCovid].select(CODON_ORDER).to_numpy().T), dfsWithAdjustmentTPM[humanCovid]["expression"], s=0.5, label = "Similarity to SARS-CoV-2 Codon Bias (Expression Weighted)")
ax1.set_xlabel("Jensen-Shannon Divergence to SARS-CoV-2 Codon Bias")
ax1.set_ylabel("Gene Expression (TPM)")
ax1.set_title("Relationship Between Gene Expression and Codon Bias in Infected Individual")
# ax1.set_xscale("log")
# ax1.set_yscale("log")
fig26.show()

### CAI vs expression
fig27 = plt.figure()

GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
}

human_covid_df = dfsWithAdjustmentTPM[humanCovid]
human_covid_matrix = human_covid_df.select(CODON_ORDER).to_numpy()
human_covid_expression = human_covid_df["expression"].to_numpy()

synonymous_codons_by_aa: dict[str, list[str]] = {}
for codon in CODON_ORDER:
    aa = GENETIC_CODE[codon]
    synonymous_codons_by_aa.setdefault(aa, []).append(codon)

synonymous_codons_by_aa = {
    aa: codons
    for aa, codons in synonymous_codons_by_aa.items()
    if aa != "*" and len(codons) > 1
}

reference_count = max(1, int(np.ceil(len(human_covid_df) * 0.2)))
reference_df = human_covid_df.sort("adjustedExpression", descending=True).head(reference_count)
reference_matrix = np.nan_to_num(reference_df.select(CODON_ORDER).to_numpy(), nan=0.0, posinf=0.0, neginf=0.0)

codon_weights: dict[str, float] = {}
preferred_codons: dict[str, str] = {}
for aa, codons in synonymous_codons_by_aa.items():
    codon_indices = [CODON_ORDER.index(codon) for codon in codons]
    aa_present = reference_matrix[:, codon_indices].sum(axis=1) > 0
    if not np.any(aa_present):
        continue

    ref_means = reference_matrix[aa_present][:, codon_indices].mean(axis=0)
    ref_means = np.clip(ref_means, 1e-12, None)
    max_ref = float(ref_means.max())
    if max_ref <= 0:
        continue

    preferred_codons[aa] = codons[int(np.argmax(ref_means))]
    for codon, mean_value in zip(codons, ref_means):
        codon_weights[codon] = float(mean_value / max_ref)

print(codon_weights)
codon_weights = {codon: weight for codon, weight in zip(CODON_ORDER, humanWithCovidCodonBias_expressionOrderedTop100.to_numpy()[0, :]) if weight > 0}
print(codon_weights)


def _gene_cai_cbi_fop_scores(row: np.ndarray) -> tuple[float, float, float]:
    row = np.nan_to_num(row.astype(np.float64, copy=False), nan=0.0, posinf=0.0, neginf=0.0)
    cai_scores: list[float] = []
    cbi_scores: list[float] = []
    fop_scores: list[float] = []

    for aa, codons in synonymous_codons_by_aa.items():
        codon_indices = [CODON_ORDER.index(codon) for codon in codons]
        codon_values = row[codon_indices].astype(np.float64)
        aa_total = float(codon_values.sum())
        if aa_total <= 0:
            continue

        codon_proportions = codon_values / aa_total
        cai_log_sum = 0.0
        for codon, proportion in zip(codons, codon_proportions):
            weight = max(codon_weights.get(codon, 1e-12), 1e-12)
            cai_log_sum += float(proportion) * np.log(weight)
        cai_scores.append(float(np.exp(cai_log_sum)))

        preferred_codon = preferred_codons.get(aa)
        if preferred_codon is None:
            continue
        preferred_index = CODON_ORDER.index(preferred_codon)
        preferred_proportion = float(row[preferred_index] / aa_total)
        fop_scores.append(preferred_proportion)

        codon_count = float(len(codons))
        if codon_count > 1:
            cbi_scores.append((preferred_proportion - (1.0 / codon_count)) / (1.0 - (1.0 / codon_count)))

    cai = float(np.mean(cai_scores)) if cai_scores else 0.0
    cbi = float(np.mean(cbi_scores)) if cbi_scores else 0.0
    fop = float(np.mean(fop_scores)) if fop_scores else 0.0
    return cai, cbi, fop


score_rows = np.array([_gene_cai_cbi_fop_scores(row) for row in human_covid_matrix], dtype=np.float64)
score_rows = np.nan_to_num(score_rows, nan=0.0, posinf=0.0, neginf=0.0)


score_labels = ["CAI", "CBI", "FOP"]
score_colors = ["tab:blue", "tab:green", "tab:orange"]
score_limits = []
for i, label in enumerate(score_labels):
    valid_mask = np.isfinite(score_rows[:, i]) & np.isfinite(human_covid_expression) & (human_covid_expression > 0)
    if np.any(valid_mask):
        score_limits.append((float(np.nanmin(score_rows[valid_mask, i])), float(np.nanmax(score_rows[valid_mask, i]))))
    else:
        score_limits.append((0.0, 1.0))

axes = [fig27.add_subplot(131), fig27.add_subplot(132), fig27.add_subplot(133)]
for idx, ax in enumerate(axes):
    valid_mask = np.isfinite(score_rows[:, idx]) & np.isfinite(human_covid_expression) & (human_covid_expression > 0)
    ax.scatter(score_rows[valid_mask, idx], human_covid_expression[valid_mask], s=0.5, alpha=0.6, color=score_colors[idx])
    ax.set_xlabel(score_labels[idx])
    if idx == 0:
        ax.set_ylabel("Gene Expression (TPM)")
    ax.set_title(f"{score_labels[idx]} vs Expression")
    ax.set_xlim(score_limits[idx])
    ax.set_yscale("log")

fig27.suptitle("Codon Bias Indices vs Gene Expression in the SARS-CoV-2 Infected Human Sample")
fig27.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
fig27.show()

fig28 = plt.figure()
ax1 = fig28.add_subplot(111)
arr = np.array((np.array(covidDates), divergencesBetweenCovidAndHumanWithCovid))
bins = [np.where(arr[0] < 50), 
        np.where((arr[0] >= 50) & (arr[0] < 100)), 
        np.where((arr[0] >= 100) & (arr[0] < 150)), 
        np.where((arr[0] >= 150) & (arr[0] < 200)), 
        np.where((arr[0] >= 200) & (arr[0] < 250)),
        np.where((arr[0] >= 250))]
bin_labels = ["0-50 days", "50-100 days", "100-150 days", "150-200 days", "200-250 days", "250+ days"]
for i in range(len(bins)):
    if len(bins[i][0]) > 0:
        ax1.violinplot(arr[1][bins[i]], positions=[i], widths=0.8)
        ax1.set_xticks(np.arange(len(bin_labels)))
        ax1.set_xticklabels(bin_labels, rotation=45)
        ax1.set_ylabel("Jensen-Shannon Divergence")
        ax1.set_title("Distribution of Divergence from Covid Infected Human Codon Bias Over Time")
fig28.tight_layout()
fig28.show()


# ax10.set_xlabel("Weeks since Week 1 of 2020")
# ax10.set_ylabel("Average Codon Bias of Top 100 Genes")
# ax10.set_zlabel("Difference in Codon Proportions (Top 100 Genes)")
# ax10.legend()
# fig6.show()

# dfs_expressionAdjusted[0]

# arr = np.zeros((500, 62), dtype=np.float64)
# for j in range(0, 4999, 10):
#     arr[j//10] = dfsWithAdjustmentTPM[i].select(CODON_ORDER)[:j+1].mean() # Average codon bias of top j genes
# ## I want a numpy way to calculate Jensen-Shannon Divergence for each codon proportions in covidDates to every row in arr
# arr.shape
# X, Y = np.meshgrid(covidDates, arr)
# X.shape
# Y.shape

# X.reshape(500, 62, 3646)
# JSD(X, Y)

# divergence = compare(covid_small, arr, JSD) 
