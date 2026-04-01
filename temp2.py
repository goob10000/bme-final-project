import polars as pl
import matplotlib.pyplot as plt
import numpy as np

df = pl.read_csv("Expression/Human1/GSE318035_raw_counts.csv")

#"Gene_ID,COVID_1,COVID_2,COVID_3,COVID_4,COVID_5,Flu_1,Flu_2,Flu_3,Flu_4,Flu_5,Flu_6,Flu_7,Flu_8,Flu_9,Flu_10"

meanCOVID = df.select(pl.col("COVID_1", "COVID_2", "COVID_3", "COVID_4", "COVID_5")).mean_horizontal()
meanFlu = df.select(pl.col("Flu_1", "Flu_2", "Flu_3", "Flu_4", "Flu_5", "Flu_6", "Flu_7", "Flu_8", "Flu_9", "Flu_10")).mean_horizontal()
df = df.with_columns(
    meanCOVID.alias("meanCOVID"),
    meanFlu.alias("meanFlu")
)

plt.scatter(np.log1p(df["meanCOVID"]), np.log1p(df["meanFlu"]), s=0.5)
plt.xlabel("Mean COVID Counts")
plt.ylabel("Mean Flu Counts")
plt.title("Mean COVID vs Mean Flu Counts")
plt.show()


