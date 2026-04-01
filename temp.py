import json, numpy as np
from codon_frequency import build_codon_frequency_matrix, CODON_ORDER, GENETIC_CODE
import polars as pl
import matplotlib.pyplot as plt

gene_dict = json.load(open("Genomes/Human/gene_id_to_sequence_transcripts.json"))
gene_ids, matrix = build_codon_frequency_matrix(gene_dict)
# matrix[i, j] = count of codon CODON_ORDER[j] in gene gene_ids[i]

gene_codon_freqs = dict(zip(gene_ids, matrix))  # dict mapping gene_id to codon frequency array

gene_codon_freqs["ENSG00000139618"]  # codon frequency array for BRCA2 gene

matrix


expression = pl.read_csv("Expression/Human1/GSE318035_raw_counts.csv")
fig = plt.figure()
ax1 = fig.add_subplot(4, 1, 1)
ax1.hist(expression["COVID_1"].log10(), bins=500, log=True)
ax2 = fig.add_subplot(4, 1, 2)
ax2.hist(expression["COVID_2"].log10(), bins=500, log=True)
ax3 = fig.add_subplot(4, 1, 3)
ax3.hist(expression["COVID_3"].log10(), bins=500, log=True)
ax4 = fig.add_subplot(4, 1, 4)
ax4.hist(expression["COVID_4"].log10(), bins=500, log=True)
plt.show()

fig = plt.figure()
ax1 = fig.add_subplot(2, 1, 1)
ax1.hist(expression["COVID_1"].log10(), bins=500, log=True)
ax2 = fig.add_subplot(2, 1, 2)
ax2.hist(expression["Flu_1"].log10(), bins=500, log=True)
plt.show()