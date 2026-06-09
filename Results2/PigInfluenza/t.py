import polars as pl
df = pl.read_csv("Results2/PigInfluenza/ExpressionData.csv")
df.columns
df = df.select("Gene_ID", "Total_Counts", "Gene_Name", "Mean_Per_Cell")
df.write_csv("Results2/PigInfluenza/ExpressionData.csv")