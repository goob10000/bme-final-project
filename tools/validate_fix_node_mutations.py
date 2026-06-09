#!/usr/bin/env python3
"""Validate and fix Results/*/node_mutations.parquet files.

Recomputes missing trinucleotide contexts and codon/AA fields for Covid output
using the reference and GTF. Writes fixed parquet (backups original).
"""
import os
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
import polars as pl
import re


def load_reference(path):
    for r in SeqIO.parse(path, 'fasta'):
        return str(r.seq)
    raise FileNotFoundError(path)


def parse_gtf(gtf_path):
    genes = {}
    with open(gtf_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 9:
                continue
            if cols[2] != 'CDS':
                continue
            start = int(cols[3]); end = int(cols[4])
            m = re.search(r'gene_id "([^\"]+)"', cols[8])
            if not m:
                continue
            gene = m.group(1)
            genes.setdefault(gene, []).append((start, end))
    for g in genes:
        genes[g] = sorted(genes[g])
    return genes


def build_gene_position_map(genes_ranges):
    gene_pos = {}
    for g, ranges in genes_ranges.items():
        positions = []
        for start, end in ranges:
            positions.extend(list(range(start, end + 1)))
        gene_pos[g] = positions
    return gene_pos


def trinuc_context_safe(ref_seq, pos):
    L = len(ref_seq)
    if pos is None or pos < 1 or pos > L:
        return None
    i = pos - 1
    left = ref_seq[i - 1] if i - 1 >= 0 else 'N'
    mid = ref_seq[i] if 0 <= i < L else 'N'
    right = ref_seq[i + 1] if i + 1 < L else 'N'
    tri = (left + mid + right).upper()
    base = tri[1]
    if base in ('A', 'G'):
        comp = str(Seq(tri).complement())
        tri = comp[::-1]
    return tri


def codon_and_aa_for_nuc(ref_seq, gene_map, pos):
    for gene, positions in gene_map.items():
        if pos in positions:
            idx = positions.index(pos)
            codon_index = idx // 3
            codon_positions = positions[codon_index * 3: codon_index * 3 + 3]
            if len(codon_positions) < 3:
                return gene, codon_positions, None, None
            codon = ''.join(ref_seq[p - 1] for p in codon_positions)
            aa = str(Seq(codon).translate())
            return gene, codon_positions, codon.upper(), aa
    return None, None, None, None


def recompute_covid(covid_parquet, ref_fasta, gtf_path):
    print('Loading', covid_parquet)
    df = pl.read_parquet(covid_parquet)
    before_nulls = {c: df.lazy().filter(pl.col(c).is_null()).select(pl.count()).collect().to_series()[0] for c in df.columns}
    print('Null counts before:', {k: int(v) for k, v in before_nulls.items()})

    ref_seq = load_reference(ref_fasta)
    genes = parse_gtf(gtf_path)
    gene_map = build_gene_position_map(genes)

    rows = df.to_dicts()
    fixed = 0
    for r in rows:
        if r.get('mut_type') == 'nuc':
            pos = r.get('pos')
            # fix ctx96
            if not r.get('ctx96') or r.get('ctx96') in ('NNN', None):
                tri = trinuc_context_safe(ref_seq, pos)
                if tri:
                    r['ctx96'] = tri
                    fixed += 1
            # fix codon_from / aa_from if missing
            if not r.get('codon_from'):
                gene, codon_positions, codon, aa = codon_and_aa_for_nuc(ref_seq, gene_map, pos)
                if codon:
                    r['codon_from'] = codon
                    r['aa_from'] = aa
                    r['gene'] = gene if r.get('gene') is None else r.get('gene')
            # fix codon_to / aa_to when alt available and codon_positions present
            if r.get('codon_from') and not r.get('codon_to') and r.get('alt'):
                # reconstruct codon_positions
                gene = r.get('gene')
                positions = gene_map.get(gene) if gene else None
                if positions and pos in positions:
                    idx = positions.index(pos)
                    codon_index = idx // 3
                    codon_positions = positions[codon_index * 3: codon_index * 3 + 3]
                    if len(codon_positions) == 3:
                        codon_list = list(r['codon_from'])
                        try:
                            idx_in_codon = codon_positions.index(pos)
                            codon_list[idx_in_codon] = r.get('alt') or codon_list[idx_in_codon]
                            codon_to = ''.join(codon_list).upper()
                            r['codon_to'] = codon_to
                            r['aa_to'] = str(Seq(codon_to).translate())
                            fixed += 1
                        except Exception:
                            pass
        # attempt to fill missing time_year/time_iso by propagating from tips to internal nodes
        # build node_id -> times mapping
        node_times = {}
        for r in rows:
            nid = r.get('node_id')
            ty = r.get('time_year')
            if ty is not None:
                node_times.setdefault(nid, []).append(ty)
        # average where multiple
        node_time_avg = {nid: sum(v)/len(v) for nid, v in node_times.items()}
        # fill rows lacking time_year using node averages if available
        for r in rows:
            if r.get('time_year') is None:
                nid = r.get('node_id')
                if nid in node_time_avg:
                    r['time_year'] = node_time_avg[nid]
                    fixed += 1
            # derive time_iso from time_year if missing
            if (not r.get('time_iso')) and r.get('time_year') is not None:
                y = float(r['time_year'])
                yr = int(y)
                mon = int((y - yr) * 12) + 1
                if mon < 1: mon = 1
                if mon > 12: mon = 12
                r['time_iso'] = f"{yr:04d}-{mon:02d}"

    # write backup and overwrite parquet
    bak = covid_parquet + '.bak'
    if not os.path.exists(bak):
        os.rename(covid_parquet, bak)
    new_df = pl.DataFrame(rows)
    new_df.write_parquet(covid_parquet)
    after_nulls = {c: new_df.lazy().filter(pl.col(c).is_null()).select(pl.count()).collect().to_series()[0] for c in new_df.columns}
    print('Fixed approx rows:', fixed)
    print('Null counts after:', {k: int(v) for k, v in after_nulls.items()})


def main():
    covid_parquet = 'Results/Covid/node_mutations.parquet'
    covid_ref = 'Genomes/Covid/covidGenomeFasta.txt'
    covid_gtf = 'Genomes/Covid/MT291829_cds.gtf'
    if os.path.exists(covid_parquet):
        recompute_covid(covid_parquet, covid_ref, covid_gtf)
    else:
        print('Covid parquet not found:', covid_parquet)


if __name__ == '__main__':
    main()
