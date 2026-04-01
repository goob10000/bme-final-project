import polars as pl
import numpy as np
import matplotlib.pyplot as plt

species = ["Bat_MyotisLucifugus", 
           "Bat_RhinolophusFerrumequinum", 
           "Human1",
           "Human2",
           "Mouse",
           "Pig",
           "Covid"]
nonCovidSpecies = ["Bat_MyotisLucifugus", 
           "Bat_RhinolophusFerrumequinum", 
           "Human1",
           "Human2",
           "Mouse",
           "Pig"]
paths = ["Results/Bat_MyotisLucifugus/topn_codon_proportions_Control.csv", 
         "Results/Bat_RhinolophusFerrumequinum/topn_codon_proportions_Control.csv", 
         "Results/Human1/topn_codon_proportions_COVID_1.csv",
         "Results/Human2/topn_codon_proportions_H1.csv",
         "Results/Mouse/topn_codon_proportions_GSM8741172_WT_1.csv",
         "Results/Pig/topn_codon_proportions_GSM7286742_Control_1.csv",
         "Results/Covid/topn_codon_proportions_Count.csv"]

dfs = {species[i]: pl.read_csv(paths[i]) for i in range(len(species))}

vectors5 = 0
vectors10 = 1
vectors50 = 2
vectors100 = 3
vectors500 = 4

vec = [5,10,15,20,25,30,35,40,45,50,100,200,300,400,500]

codons = [x+y+z for x in "ACTG" for y in "ACTG" for z in "ACTG"]

codonIndex = {codon: i for i, codon in enumerate(codons)}

arr = np.zeros((len(species), len(vec), 64))

dfs[species[5]].filter(pl.col("top_n") == 5).shape[0]

for s in nonCovidSpecies:
    print(f"{s=}")
    for i, val in enumerate(vec):
        for row in dfs[s].filter(pl.col("top_n") == val).sort("codon").iter_rows():
            codon = row[3]
            codon_proportion = row[6]
            arr[species.index(s), i, codonIndex[codon]] = codon_proportion

for row in dfs["Covid"].filter(pl.col("top_n") == 5).sort("codon").iter_rows():
    codon = row[3]
    codon_proportion = row[6]
    for i in range(len(vec)):
        arr[species.index("Covid"), i, codonIndex[codon]] = codon_proportion

# for i in range(arr.shape[1]):
#     plt.figure(figsize=(10, 6))
#     for j in range(arr.shape[0]):
#         plt.plot(arr[j, i], label=species[j])
#     plt.title(f"Top {vec[i]} Codon Proportions")
#     plt.xlabel("Codon Index")
#     plt.ylabel("Proportion")
#     plt.legend()
#     plt.show()

div = np.zeros((len(species)-1, len(vec)))

for i in range(arr.shape[1]):
    dist = arr[:, i, :]
    for s in nonCovidSpecies:
        covid_dist = arr[-1, i, :]
        species_dist = arr[species.index(s), i, :]
        # print(f"covid_dist={covid_dist.sum()}")
        # print(f"{species} dist={species_dist.sum()}")
        # print(f"{s + '-' + str(vec[i]) + ': ' + str(abs((np.log1p(covid_dist)-np.log1p(species_dist)).sum()))}")
        kl_divergence = np.sum(np.where(covid_dist != 0, covid_dist * np.log(covid_dist / (species_dist + 1e-10)), 0))
        print(f"KL Divergence between Covid and {s} for top {vec[i]} genes: {kl_divergence}")
        div[species.index(s), i] = kl_divergence

plt.plot(vec, div.T)
plt.xlabel("Top N Genes")
plt.ylabel("KL Divergence from Covid")
plt.title("KL Divergence of Codon Proportions from Covid")
plt.legend(nonCovidSpecies)
plt.show()


l = [x.filter(pl.col("aa") == "A") for x in dfs.values()]

fig = plt.figure(figsize=(10, 6), layout="constrained")
for i, c in enumerate(l[0]["codon"]):
    print(c)
    a = fig.add_subplot(2, 2, i+1)
    for j, s in enumerate(species):
        a.plot(l[j].filter(pl.col("codon") == c)["top_n"], l[j].filter(pl.col("codon") == c)["codon_proportion"], label=s, marker="o")
        a.set_title(f"Codon {c}")
        # a.set_ylim(0, 0.7)
        a.set_xlabel("log(Top N Genes)")
        a.set_ylabel("Codon Proportion")
        # a.set_xticks(np.exp([0, 1, 2, 3, 4, 5]))
        a.set_xscale("log")
        # a.set_xticks(range(len(species)))
        # a.set_xticklabels(species, rotation=45)

plt.legend()
plt.show()
