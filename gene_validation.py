import polars as pl
import matplotlib.pyplot as plt
import numpy as np
import json


covid = "C:/Projects/UCSC/BME205/Final Project/Expression/Covid/human_geneid_count_pairs_normalized.csv"
human1 = "C:/Projects/UCSC/BME205/Final Project/Expression/Human1/GSE318035_raw_counts.csv"
human2 = "C:/Projects/UCSC/BME205/Final Project/Expression/Human2/GSE275264_raw_healthy.csv"
mouse = "C:/Projects/UCSC/BME205/Final Project/Expression/Mouse/ExpressionData.csv"
pig = "C:/Projects/UCSC/BME205/Final Project/Expression/Pig/ExpressionData.csv"
bat_rhinolophus = "C:/Projects/UCSC/BME205/Final Project/Expression/Bat_RhinolophusFerrumequinum/ExpressionData.csv"
bat_myotis = "C:/Projects/UCSC/BME205/Final Project/Expression/Bat_MyotisLucifugus/ExpressionData.csv"

covid_cds = "C:/Projects/UCSC/BME205/Final Project/Genomes/Covid/gene_id_to_sequence_transcripts.json"
human_cds = "C:/Projects/UCSC/BME205/Final Project/Genomes/Human/gene_id_to_sequence_transcripts.json"
mouse_cds = "C:/Projects/UCSC/BME205/Final Project/Genomes/Mouse/gene_id_to_sequence_transcripts.json"
pig_cds = "C:/Projects/UCSC/BME205/Final Project/Genomes/Pig/gene_id_to_sequence_transcripts.json"
myotis_cds = "C:/Projects/UCSC/BME205/Final Project/Genomes/Bat_MyotisLucifugus/gene_id_to_transcripts.json"
rhino_cds = "C:/Projects/UCSC/BME205/Final Project/Genomes/Bat_RhinolophusFerrumequinum/gene_id_to_transcripts.json"

organisms = [covid, human1, human2, mouse, pig, bat_rhinolophus, bat_myotis]
organismNames = ["Covid", "Human1", "Human2", "Mouse", "Pig", "Bat_RhinolophusFerrumequinum", "Bat_MyotisLucifugus"]
cds_list = [covid_cds, human_cds, human_cds, mouse_cds, pig_cds, rhino_cds, myotis_cds]
cdss = [json.load(open(file)) for file in cds_list]


dfs = [pl.read_csv(file) for file in organisms]
def top_5_genes(df, cds):
    col_name = df.columns[1]
    gene_id = df.columns[0]
    df = df.filter(pl.col(gene_id).is_in(list(cds.keys())))
    top_5 = df.sort(col_name, descending=True).head(5)
    return top_5.select(gene_id).to_series().to_list()

top_genes = [top_5_genes(df, cds) for df, cds in zip(dfs, cdss)]

for organism, genes in zip(organismNames, top_genes):
    print(f"Organism: {organism}")
    for gene in genes:
        print(f"Gene: {gene}")
