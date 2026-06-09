import re
import polars as pl
from Bio import SeqIO
import numpy as np
from functools import reduce
from collections import Counter
from AnalysisFunctions import CODON_ORDER
from codon_frequency import GENETIC_CODE
import json


reference_influenza_genome_path = "Genomes/Influenza/genbank reference sequence/ncbi_dataset/data/GCA_046377675.1/GCA_046377675.1_ASM4637767v1_genomic.fna"
reference_influenza_gtf_path = "Genomes/Influenza/genbank reference sequence/ncbi_dataset/data/GCA_046377675.1/genomic.gtf"

gtf = pl.read_csv(reference_influenza_gtf_path, separator="\t", skip_lines=3, new_columns=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
gtf_cds = gtf.filter(pl.col("feature") == "CDS")
gtf_cds = gtf_cds.with_columns(
    gene_id=pl.col("attribute").str.extract(r'gene_id "([^"]+)"')
)


## Example fasta header:
# >LC660656.1 Influenza A virus (H1N1) A/Iwate/1130/2009 PB2 gene for polymerase PB2, complete cds
# ATGGAGAGAATAAAAGAACTGAGAGATCTAATGTCGCAGTCCCGCACTCGCGAGATACTCACTAAGACCACTGTGGACCA
# TATGGCCATAATCAAAAAGTACACATCAGGAAGGCAAGAGAAGAACCCCGCACTCAGAATGAAGTGGATGATGGCAATGA
# GATACCCAATTACAGCAGACAAGAGAATAATGGACATGATTCCAGAGAGGAATGAACAAGGACAAACCCTCTGGAGCAAA

# Take cds from gtf_cds. If gene doesn't exist in the dictionary, add it with the sequence.
# If it does exist, append the sequence.
# Key for the dictionary is the gene_id which can be found as the first element in the attribute column, which is a string of key-value pairs separated by ";". The gene_id is in the format "gene_id "GENE_NAME"".
gene_ids = gtf_cds["gene_id"].unique()

gene_sequences = {}
for gene_id in gene_ids:
    gene_sequences[gene_id] = ""

# Iterate through gtf_cds and append sequences to the corresponding gene_id in gene_sequences
for row in gtf_cds.iter_rows():
    gene_id = row[9]  # gene_id is the 10th column (index 9)
    seqname = row[0]  # seqname is the 1st column (index 0)
    start = row[3]  # start is the 4th column (index 3)
    end = row[4]  # end is the 5th column (index 4)
    strand = row[6]  # strand is the 7th column (index 6)

    # Extract the sequence from the reference genome
    for record in SeqIO.parse(reference_influenza_genome_path, "fasta"):
        if record.id == seqname:
            if strand == "+":
                gene_sequences[gene_id] += str(record.seq[start-1:end+3])  # -1 because of 0-based indexing
            else:
                gene_sequences[gene_id] += str(record.seq[start-1:end+3].reverse_complement())

gene_lengths = {gene_id: len(seq) for gene_id, seq in gene_sequences.items()}
length_to_gene = {length: gene_id for gene_id, length in gene_lengths.items()}

gene_lengths

non_stop_codons = []
nucleotides = ['A', 'T', 'C', 'G']
for i in nucleotides:
    for j in nucleotides:
        for k in nucleotides:
            codon = i + j + k
            if codon not in ['TAA', 'TAG', 'TGA']:
                non_stop_codons.append(codon)

influenza_sequences_to_analyze_path = "Genomes/Influenza/sequences.fasta"
influenza_sequences_to_analyze_meta = pl.read_csv("Genomes/Influenza/sequences.csv", schema_overrides={"Accession": pl.String, "GenBank_RefSeq": pl.String, "Assembly": pl.String, "Organism_Name": pl.String, "Species": pl.String, "Genus": pl.String, "Genotype": pl.String, "Isolate": pl.String, "Segment": pl.String, "Length": pl.Int64, "Nuc_Completeness": pl.String, "Geo_Location": pl.String, "Country": pl.String, "Host": pl.String, "Collection_Date": pl.String, "Molecule_type": pl.String}, separator=",")

influenza_sequences_to_analyze_meta.columns

# First few lines of sequences.csv:
# Accession,GenBank_RefSeq,Assembly,Organism_Name,Species,Genus,Genotype,Isolate,Segment,Length,Nuc_Completeness,Geo_Location,Country,Host,Collection_Date,Molecule_type
# NC_026431,RefSeq,GCF_001343785.1,Influenza A virus (A/California/07/2009(H1N1)),Alphainfluenzavirus influenzae,Alphainfluenzavirus,H1N1,,7,982,complete,USA: California state,USA,Homo sapiens,2009-04-09,ssRNA(-)
# NC_026432,RefSeq,GCF_001343785.1,Influenza A virus (A/California/07/2009(H1N1)),Alphainfluenzavirus influenzae,Alphainfluenzavirus,H1N1,,8,863,complete,USA: California state,USA,Homo sapiens,2009-04-09,ssRNA(-)
# NC_026433,RefSeq,GCF_001343785.1,Influenza A virus (A/California/07/2009(H1N1)),Alphainfluenzavirus influenzae,Alphainfluenzavirus,H1N1,,4,1701,complete,USA: California state,USA,Homo sapiens,2009-04-09,ssRNA(-)

codonPossibilities = json.load(open("codon_possibilities.json", "r"))

def codon_proportion(sequence):
    codonCounter = Counter()
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        codonCounter[codon] += 1
    codonCounts = np.array([codonCounter[codon] for codon in CODON_ORDER])
    codonProportions = np.zeros_like(codonCounts, dtype=np.float64)
    c = Counter()
    for j, codon in enumerate(CODON_ORDER):
        c[GENETIC_CODE[codon]] += codonCounts[j]
    for j, codon in enumerate(CODON_ORDER):
        if c[GENETIC_CODE[codon]] > 0:
            codonProportions[j] = codonCounts[j] / c[GENETIC_CODE[codon]] / codonPossibilities[codon]
    return codonProportions


cols = influenza_sequences_to_analyze_meta.columns + CODON_ORDER + ["geneID", "orf"]
outputPath = "influenzaAlignmentM2.csv"
with open(outputPath, "w") as f:
    f.write(",".join(cols) + "\n")
    for i, sequence in enumerate(SeqIO.parse(influenza_sequences_to_analyze_path, "fasta")):
        seq = str(sequence.seq)
        if i%10000 == 0:
            print(f"Processing sequence {i} / {len(influenza_sequences_to_analyze_meta)}")
        meta = influenza_sequences_to_analyze_meta[i]
        if meta["Segment"][0] != "7":
            continue
        if len(seq) != 982:
            continue
        M2Orf = seq[:26] + seq[715:]
        gene_id = "M2"
        proportion = codon_proportion(M2Orf)
        if gene_id == None:
            continue
        dfAdding = meta
        dfAdding2 = pl.DataFrame({
            "geneID": [gene_id],
            "orf" : [M2Orf]
        })
        dfAdding3 = pl.DataFrame({CODON_ORDER[i]: [proportion[i]] for i in range(len(CODON_ORDER))})
        dfAdding = dfAdding.join(dfAdding2, how="cross").join(dfAdding3, how="cross")
        x = f.write(",".join([str(dfAdding[x][0]).replace(",", "") if dfAdding[x][0] is not None else "" for x in cols]) + "\n")
