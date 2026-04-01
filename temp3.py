import polars as pl
import matplotlib.pyplot as plt
import numpy as np

dfPrec = pl.read_csv("Expression/Covid/all_runs/human_geneid_count_pairs.csv")
df2 = pl.DataFrame(np.cast[int](dfPrec.drop("GeneID").to_numpy().T), schema=dfPrec["GeneID"].to_list())
keep = df2.to_numpy().sum(axis=1) > 1
df3 = pl.DataFrame(df2.to_numpy()[keep], schema=dfPrec["GeneID"].to_list())


normalized = df3  / df3.sum_horizontal()

count = np.cast[int]((normalized.mean()*10000).to_numpy()).flatten()
genes = normalized.columns
dfOut = pl.DataFrame({"GeneID": genes, "Count": count})
dfOut.write_csv("Expression/Covid/human_geneid_count_pairs_normalized.csv")

plt.plot(normalized)
plt.legend(labels=normalized.columns)
plt.show()