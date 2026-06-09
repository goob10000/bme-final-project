from xgboost import XGBRegressor, plot_tree, XGBClassifier
import numpy as np
import polars as pl
import pandas as pd
import matplotlib.pyplot as plt
from AnalysisFunctions import CODON_ORDER, adjustTPM, joinMatrixGenesAndCounts, joinMatrixGenesAndCountsWithMT, weightByExpression, weightByAdjustedExpression, diff, compare, JSD
import os
from datetime import datetime
from explainerdashboard import RegressionExplainer, ExplainerDashboard, ClassifierExplainer
mtGenesList = pl.read_csv("mtGenes.tsv", separator="\t").filter(pl.col("gene_biotype") == "protein_coding")

species = [d for d in os.listdir("Results2/") if (d != "Covid" and d != "Influenza")]
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

covid_small = covidCodonBiasMatrix.mean(axis=2) # Average across genes to get overall codon proportions at each timepoint

dfsWithAdjustmentTPM = adjustTPM(dfs, species, expression)
dfsWithAdjustmentAndMT = adjustTPM(dfs_withMT, species, expression)

humanCovid = 4
humanHealthy = 5
humanInfluenza = 6

dfs_expressionAdjusted = [weightByAdjustedExpression(df) for df in dfsWithAdjustmentTPM]
dfs_expressionAdjustedWithMT = [weightByAdjustedExpression(df) for df in dfsWithAdjustmentAndMT]
dfs_expressionUnadjusted = [weightByExpression(df) for df in dfsWithAdjustmentTPM]
dfs_expressionUnadjustedWithMT = [weightByExpression(df) for df in dfsWithAdjustmentAndMT]

hCN = dfs_expressionAdjusted[humanCovid]
hHN = dfs_expressionAdjusted[humanHealthy]
hIN = dfs_expressionAdjusted[humanInfluenza]
hCN_arr = dfs_expressionAdjusted[humanCovid]
hHN_arr = dfs_expressionAdjusted[humanHealthy]
hIN_arr = dfs_expressionAdjusted[humanInfluenza]
hCN_arrNoTPM = dfs_expressionUnadjusted[humanCovid]
hHN_arrNoTPM = dfs_expressionUnadjusted[humanHealthy]
hIN_arrNoTPM = dfs_expressionUnadjusted[humanInfluenza]

hCC = compare(covid_small, hCN_arr, JSD)
hHC = compare(covid_small, hHN_arr, JSD)
hIC = compare(covid_small, hIN_arr, JSD)

influenza_matrix_path = "Results2/Influenza/codon_proportions.npy"
influenza_meta_path = "Results2/Influenza/codon_proportions_metadata.csv"

influenza_dates = None
influenza_hCC = None
influenza_hHC = None
influenzaGenotypes = None
influenzaHosts = None
influenzaMeta = None

if os.path.exists(influenza_matrix_path) and os.path.exists(influenza_meta_path):
    influenzaCodonBiasMatrix = np.load(influenza_matrix_path)
    influenzaMeta = pl.read_csv(influenza_meta_path)
    influenzaCodonSmall = influenzaCodonBiasMatrix.mean(axis=1)  # (rows, 62)
    influenzaGenotypes = influenzaMeta["Genotype"]
    influenzaHosts = influenzaMeta["Host"]
    influenza_dates = [datetime.fromordinal(int(x)) for x in influenzaMeta["date_ordinal"].to_list()]
    influenza_hCI = compare(influenzaCodonSmall, hCN_arr, JSD)
    influenza_hHI = compare(influenzaCodonSmall, hHN_arr, JSD)
    influenza_hII = compare(influenzaCodonSmall, hIN_arr, JSD)
influenzaColorsGenotypes = {}
hosts_list: list[str] = []
influenzaColorsHosts = {}
influenzaCalculatedColorsHosts = []
host_groups = [[], [], [], []]

if influenzaGenotypes is not None:
    influenzaColorsGenotypes = {genotype: color for genotype, color in zip(set(influenzaGenotypes.unique()), species_colors)}

if influenzaHosts is not None:
    hosts_list = [str(h).strip() if h is not None else "" for h in influenzaHosts.to_list()]
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


X = influenzaMeta.drop(["strain_key", "Assembly", "date_ordinal", "Organism_Name", "Genotype", "Nuc_Completeness", "segments_present", "segment_list", "row_index", "Length", "Accession", "Segment", "Isolate", "parsed_date", "Molecule_type", "Collection_Date"]).to_pandas()
# X["strain_key"]


plt.scatter(np.asarray(influenza_dates), influenza_hII, s=0.2)
plt.show()

Y = (influenza_hII < 0.045).astype(int)

# Convert the metadata-heavy feature table into a numeric matrix so XGBoost,
# SHAP, and explainerdashboard all see the same stable dtype layout.
X_d = pd.get_dummies(X.fillna("Unknown"), dummy_na=False)
X_d.shape

bst = XGBClassifier(objective="binary:logistic", eval_metric="logloss", n_estimators=50, max_depth=5, random_state=2, device="cpu")
bst.fit(X_d, Y)

# Use a real pandas sample so the one-hot feature matrix stays aligned for SHAP.
background = X_d.sample(n=min(100, len(X)), random_state=2)

# Explainerdashboard will build the tree SHAP explainer internally from these args.
explainer = ClassifierExplainer(bst, X_d, Y, shap="tree", X_background=background, shap_kwargs={"check_additivity": False})
db = ExplainerDashboard(explainer, simple=True, hide_pdp=True)
db.run()