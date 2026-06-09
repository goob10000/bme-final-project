import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import json
from codon_frequency import GENETIC_CODE, CODON_INDEX
from collections import defaultdict

CODON_ORDER = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']


mutations = pl.read_parquet("Genomes/Covid/mutations2.parquet")
mutations = mutations.filter(pl.col("mutations") != mutations["mutations"][1]) # Remove null (no mutations?)
arr = np.zeros((mutations.height,62,11), dtype=np.float64) # 62 codons (no tyr or met cause they have one codon), 11 genes
SARS_genome = pl.read_csv("Genomes/Covid/MT291829_cds.tsv", separator="\t")
genes = SARS_genome["gene_id"].to_list()

covidDict = json.load(open(r"Genomes\Covid\genes.json"))
seqs = []
# covidDict["ORF8"][list(covidDict["ORF8"].keys())[0]]["sequence"]

for i, gene in enumerate(SARS_genome["gene_id"]):
    d = covidDict[gene]
    seq = d[list(d.keys())[0]]["sequence"]
    seqs.append(seq)

df = SARS_genome.with_columns(pl.Series(seqs).alias("sequence"))

info = {
        gene: (df.filter(pl.col("gene_id") == gene)["start"][0],df.filter(pl.col("gene_id") == gene)["stop"][0],df.filter(pl.col("gene_id") == gene)["sequence"].to_list()[0]) for gene in df["gene_id"]
    }

def codon_proportions(seq: str) -> np.ndarray:
    counts = {codon: 0.0 for codon in [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"]}
    total_codons = 0
    sequence = seq.upper()
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]
        if all(base in "ACGT" for base in codon):
            counts[codon] += 1.0
            total_codons += 1
    if total_codons == 0:
        return np.zeros(len(CODON_ORDER), dtype=np.float64)
    aa_totals: dict[str, float] = defaultdict(float)
    for codon, count in counts.items():
        aa_totals[GENETIC_CODE[codon]] += float(count)
    arr = np.zeros(len(CODON_ORDER), dtype=np.float64)
    for codon, count in counts.items():
        idx = CODON_INDEX.get(codon)
        aa_total = aa_totals[GENETIC_CODE[codon]]
        if idx is not None and aa_total > 0:
            arr[idx] = float(count) / aa_total

    return arr

for i in range(mutations.height):
    # dim 1: time
    print(f"Processing mutation {i+1} of {mutations.height}")
    infoCopy = info.copy() # dict of [gene → sequence] for the current timepoint
    mutationsNow = mutations["mutations"][i]["nuc"] ## list of all mutations NucPosNuc, e.g. "C241T"
    for mutationNow in mutationsNow:
        pos = int(mutationNow[1:-1])
        for gene in genes:
            if pos >= infoCopy[gene][0] and pos <= infoCopy[gene][1]: # If the mutation is within the gene
                seq = infoCopy[gene][2]
                posInGene = pos - infoCopy[gene][0] # Position of the mutation within the gene
                codonStart = (posInGene // 3) * 3 # Start of the codon containing the mutation
                codonEnd = codonStart + 3 # End of the codon containing the mutation
                mutatedCodon = seq[codonStart:codonEnd] # The original codon sequence
                mutatedCodonList = list(mutatedCodon)
                mutatedCodonList[posInGene % 3] = mutationNow[-1] # Apply the mutation to the codon
                mutatedCodon = "".join(mutatedCodonList)
                infoCopy[gene] = (infoCopy[gene][0], infoCopy[gene][1], seq[:codonStart] + mutatedCodon + seq[codonEnd:]) # Update the sequence for this gene in the copy of the info dict
    for gene in genes:
        arr[i,:,genes.index(gene)] = codon_proportions(infoCopy[gene][2]) # Calculate codon proportions for this gene at this timepoint and store in the array
    # dim 2: codon proportions
    # dim 3: gene
np.save("Results2/Covid/codon_proportions.npy", arr)