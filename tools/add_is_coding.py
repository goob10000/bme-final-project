#!/usr/bin/env python3
"""Add `is_coding` boolean column to node_mutations Parquet files.

For Covid: use `Genomes/Covid/MT291829_cds.gtf` to map nucleotide positions to CDS.
For Influenza: mark `is_coding=True` only when `gene` column is non-null (no CDS map available).

Backs up original Parquet files before writing.
"""
import os
import re
from Bio import SeqIO
import polars as pl


def parse_gtf(gtf_path):
    genes = {}
    with open(gtf_path) as fh:
        for line in fh:
            line=line.strip()
            if not line or line.startswith('#'):
                continue
            cols=line.split('\t')
            if len(cols)<9 or cols[2] != 'CDS':
                continue
            start=int(cols[3]); end=int(cols[4])
            m=re.search(r'gene_id "([^\"]+)"', cols[8])
            if not m:
                continue
            gene=m.group(1)
            genes.setdefault(gene, []).append((start,end))
    for g in genes:
        genes[g]=sorted(genes[g])
    return genes


def build_coding_pos_set(genes_ranges):
    s=set()
    for ranges in genes_ranges.values():
        for start,end in ranges:
            s.update(range(start,end+1))
    return s


def process_covid(parquet_path, gtf_path):
    print('Processing Covid parquet', parquet_path)
    df=pl.read_parquet(parquet_path)
    genes=parse_gtf(gtf_path)
    coding_pos=build_coding_pos_set(genes)
    rows=df.to_dicts()
    for r in rows:
        is_coding=False
        if r.get('gene'):
            is_coding=True
        else:
            if r.get('mut_type')=='nuc' and r.get('pos'):
                try:
                    if int(r['pos']) in coding_pos:
                        is_coding=True
                except Exception:
                    pass
        r['is_coding']=is_coding
    bak=parquet_path+'.pre_is_coding.bak'
    if not os.path.exists(bak):
        os.rename(parquet_path, bak)
    pl.DataFrame(rows).write_parquet(parquet_path)
    print('Wrote', parquet_path)


def process_influenza(parquet_path):
    print('Processing Influenza parquet', parquet_path)
    df=pl.read_parquet(parquet_path)
    rows=df.to_dicts()
    for r in rows:
        r['is_coding'] = bool(r.get('gene'))
    bak=parquet_path+'.pre_is_coding.bak'
    if not os.path.exists(bak):
        os.rename(parquet_path, bak)
    pl.DataFrame(rows).write_parquet(parquet_path)
    print('Wrote', parquet_path)


def main():
    covid_parquet='Results/Covid/node_mutations.parquet'
    covid_gtf='Genomes/Covid/MT291829_cds.gtf'
    influenza_parquet='Results/Influenza/node_mutations.parquet'
    if os.path.exists(covid_parquet) and os.path.exists(covid_gtf):
        process_covid(covid_parquet, covid_gtf)
    else:
        print('Skipping Covid - missing files')
    if os.path.exists(influenza_parquet):
        process_influenza(influenza_parquet)
    else:
        print('Skipping Influenza - missing parquet')


if __name__=='__main__':
    main()
