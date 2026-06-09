import os
import polars as pl
import numpy as np
import json
import matplotlib.pyplot as plt

CODON_ORDER = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']

def joinMatrixGenesAndCounts (species, stuff, mtGenesList):
    '''
    Convert codon matrix into dataframe with codon as column name and then add gene_id and length columns.
    Filters out mitochondrial genes in the process.
    '''
    dfs = []
    for i in range(len(species)):
        df = pl.DataFrame({codon:col for codon, col in zip(CODON_ORDER, stuff[i]["matrix"].T)})
        df = df.with_columns(
            pl.Series(stuff[i]["gene_ids"]).alias("gene_id"),
            pl.Series(stuff[i]["counts"]).alias("length")
        )
        df = df.filter(pl.col("gene_id").is_in(list(mtGenesList["gene_id"])).not_()) # Filter out mitochondiral genes
        dfs.append(df)
    return dfs

def joinMatrixGenesAndCountsWithMT (species, stuff):
    '''
    Same as join MatrixGenesAndCounts but without filtering out mitochondrial genes.
    '''
    dfs = []
    for i in range(len(species)):
        df = pl.DataFrame({codon:col for codon, col in zip(CODON_ORDER, stuff[i]["matrix"].T)})
        df = df.with_columns(
            pl.Series(stuff[i]["gene_ids"]).alias("gene_id"),
            pl.Series(stuff[i]["counts"]).alias("length")
        )
        dfs.append(df)
    return dfs

def adjustTPM (dfs: list[pl.DataFrame], species: list[str], expression: list) -> list[pl.DataFrame]:
    '''
    Join on expression (retaining only genes with expression data), then add a new column "adjustedExpression" = expression / length, which is proportional to TPM. Then drop nulls (genes without expression data).
    '''
    dfsNew = []
    for i, _ in enumerate(species):
        dfsNew.append(dfs[i].join(expression[i], on="gene_id").with_columns(
            (pl.col("expression") / pl.col("length")).alias("adjustedExpression")).drop_nulls())
            # (pl.col("expression") / pl.col("length")).alias("adjustedExpression")))
    return dfsNew

def weightByAdjustedExpression(df: pl.DataFrame) -> pl.DataFrame:
    return (df.select(CODON_ORDER).to_numpy() * df["adjustedExpression"].to_numpy()[:, None]).sum(0)/ df["adjustedExpression"].sum()

def weightByExpression(df: pl.DataFrame) -> pl.DataFrame:
    return (df.select(CODON_ORDER).to_numpy() * df["expression"].to_numpy()[:, None]).sum(0)/ df["expression"].sum()

def _normalize_distribution(x, axis=0, eps=1e-12):
    x = np.asarray(x, dtype=np.float64)
    # Guard against tiny negative numerical noise before normalization.
    x = np.clip(x, 0.0, None)
    x = x + eps
    total = np.sum(x, axis=axis, keepdims=True)
    return x / total

def _kl_one_way(p, q, axis=0):
    return np.sum(p * (np.log(p) - np.log(q)), axis=axis)

def klDivergence (p, q, axis=0):
    p = _normalize_distribution(p, axis=axis)
    q = _normalize_distribution(q, axis=axis)
    kl_div = _kl_one_way(p, q, axis=axis) + _kl_one_way(q, p, axis=axis)
    return kl_div

def JSD (p, q, axis=0):
    p = _normalize_distribution(p, axis=axis)
    q = _normalize_distribution(q, axis=axis)
    m = 0.5 * (p + q)
    # Jensen-Shannon divergence uses one-way KL terms: KL(p||m) and KL(q||m).
    kl_pm = _kl_one_way(p, m, axis=axis)
    kl_qm = _kl_one_way(q, m, axis=axis)
    jsd = 0.5 * (kl_pm + kl_qm)
    return jsd

def kl_divergence(b, a, axis=0):
    return np.sum(np.where(a != 0, a * np.log(a / (b + 1e-10)), 0), axis=0) + np.sum(np.where(b != 0, b * np.log(b / (a + 1e-10)), 0), axis=axis)


def diff (a,b, axis=0):
    return np.sum(np.abs(a-b), axis=axis)

def compare (a, b, fun):
    return np.array([fun(a[i], b) for i in range(a.shape[0])])
