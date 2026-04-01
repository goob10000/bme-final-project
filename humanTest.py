import numpy as np
import polars as pl



df = pl.read_csv("HumanGenome/GSE318035_raw_counts.csv")

df["COVID_1"].max()
