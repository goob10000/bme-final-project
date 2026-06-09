#!/usr/bin/env python3
"""Generate per-branch mutation Parquet files for Covid and Influenza.

Writes:
- Results/Covid/node_mutations.parquet
- Results/Influenza/node_mutations.parquet

Requirements: polars, biopython
"""
import json
import os
import re
from collections import Counter, defaultdict
from datetime import datetime
import math

import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq


def load_reference(path):
    rec = None
    for r in SeqIO.parse(path, "fasta"):
        if rec is None:
            rec = r
            break
    if rec is None:
        raise FileNotFoundError(path)
    return str(rec.seq)


def parse_gtf(gtf_path):
    genes = defaultdict(list)
    with open(gtf_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 9:
                continue
            feature = cols[2]
            if feature != 'CDS':
                continue
            start = int(cols[3])
            end = int(cols[4])
            attrs = cols[8]
            m = re.search(r'gene_id "([^"]+)"', attrs)
            gene = m.group(1) if m else None
            if gene:
                genes[gene].append((start, end))
    # sort ranges
    for g in genes:
        genes[g] = sorted(genes[g])
    return genes


def build_gene_position_map(genes_ranges):
    # for each gene produce a list mapping CDS index (1-based nucleotide) -> genomic pos
    gene_pos = {}
    for g, ranges in genes_ranges.items():
        positions = []
        for start, end in ranges:
            positions.extend(list(range(start, end + 1)))
        gene_pos[g] = positions
    return gene_pos


def traverse_nextstrain_json(json_path):
    with open(json_path) as fh:
        data = json.load(fh)
    root = data.get('tree')
    nodes = []
    counter = 0

    def recurse(node, parent_id=None):
        nonlocal counter
        node_id = f'NODE_{counter:06d}'
        counter += 1
        node_name = node.get('name')
        node_attrs = node.get('node_attrs', {})
        branch_attrs = node.get('branch_attrs', {})
        nodes.append({'node_id': node_id, 'name': node_name,
                      'node_attrs': node_attrs, 'branch_attrs': branch_attrs,
                      'parent_id': parent_id})
        for c in node.get('children', []):
            recurse(c, node_id)

    recurse(root, None)
    return nodes


def parse_nuc_change(s):
    # examples: C241T, 241T, T104A
    m = re.match(r'([ACGT])?(\d+)([ACGT])?', s)
    if not m:
        return None
    ref, pos, alt = m.group(1), int(m.group(2)), m.group(3)
    return {'ref': ref, 'pos': pos, 'alt': alt}


def trinuc_context(ref_seq, pos):
    # pos is 1-based
    L = len(ref_seq)
    # guard against out-of-range positions
    if pos is None or pos < 1 or pos > L:
        return 'NNN'
    i = pos - 1
    left = ref_seq[i - 1] if i - 1 >= 0 else 'N'
    mid = ref_seq[i] if 0 <= i < L else 'N'
    right = ref_seq[i + 1] if i + 1 < L else 'N'
    tri = (left + mid + right).upper()
    # canonicalize to pyrimidine-centered (C/T) convention for 96 context
    base = tri[1]
    if base in ('A', 'G'):
        # complement
        comp = str(Seq(tri).complement())
        tri = comp[::-1]  # reverse to keep 5'->3'
    return tri


def codon_and_aa_for_nuc(ref_seq, gene_map, pos):
    # pos genomic 1-based
    for gene, positions in gene_map.items():
        if pos in positions:
            idx = positions.index(pos)  # 0-based nucleotide index within CDS
            codon_index = idx // 3
            codon_positions = positions[codon_index * 3: codon_index * 3 + 3]
            if len(codon_positions) < 3:
                return gene, None, None, None
            codon = ''.join(ref_seq[p - 1] for p in codon_positions)
            aa = str(Seq(codon).translate())
            return gene, codon_positions, codon.upper(), aa
    return None, None, None, None


def build_influenza_consensus(fasta_path, n=5, max_n=10):
    seqs = list(SeqIO.parse(fasta_path, 'fasta'))
    if len(seqs) == 0:
        raise FileNotFoundError(fasta_path)
    take = min(n, len(seqs))
    for attempt in range(take, max_n + 1):
        selection = seqs[:attempt]
        lengths = {len(s.seq) for s in selection}
        if len(lengths) == 1:
            L = lengths.pop()
            cons = []
            for i in range(L):
                bases = [str(s.seq[i]).upper() for s in selection]
                c = Counter(bases).most_common(1)[0][0]
                cons.append(c)
            return ''.join(cons)
    # fallback: use first sequence
    return str(seqs[0].seq)


def float_year_from_date_string(date_str):
    # try to parse YYYY-MM-DD or YYYY-MM or decimal
    try:
        if isinstance(date_str, (int, float)):
            return float(date_str)
        if re.match(r'^\d{4}\.\d+$', str(date_str)):
            return float(date_str)
        if re.match(r'^\d{4}-\d{2}-\d{2}$', date_str):
            dt = datetime.strptime(date_str, '%Y-%m-%d')
            return dt.year + (dt.timetuple().tm_yday - 1) / (366 if dt.year % 4 == 0 else 365)
        if re.match(r'^\d{4}-\d{2}$', date_str):
            dt = datetime.strptime(date_str, '%Y-%m')
            return dt.year + (dt.timetuple().tm_yday - 1) / (366 if dt.year % 4 == 0 else 365)
    except Exception:
        return None
    return None


def process_covid(json_path, ref_fasta, gtf_path, out_path):
    print('Processing Covid tree...')
    ref_seq = load_reference(ref_fasta)
    genes = parse_gtf(gtf_path)
    gene_map = build_gene_position_map(genes)
    nodes = traverse_nextstrain_json(json_path)
    rows = []
    for n in nodes:
        node_id = n['node_id']
        name = n.get('name')
        node_attrs = n.get('node_attrs', {})
        branch = n.get('branch_attrs', {})
        # time extraction
        time_iso = None
        time_year = None
        na = node_attrs
        # try date, num_date, epiweek
        if 'num_date' in na:
            time_year = float(na['num_date']['value'] if isinstance(na['num_date'], dict) else na['num_date'])
        if 'date' in na:
            try:
                d = na['date']
                if len(d) >= 7:
                    time_iso = d[:7]
                else:
                    time_iso = d
            except Exception:
                pass
        if time_year is None and 'date' in na:
            time_year = float_year_from_date_string(na['date'])

        muts = branch.get('mutations', {})
        # process nuc
        for s in muts.get('nuc', []) if isinstance(muts.get('nuc', []), list) else []:
            parsed = parse_nuc_change(s)
            if not parsed:
                continue
            pos = parsed['pos']
            ref = parsed['ref']
            alt = parsed['alt']
            ctx = trinuc_context(ref_seq, pos)
            gene, codon_positions, codon, aa = codon_and_aa_for_nuc(ref_seq, gene_map, pos)
            codon_from = codon
            if codon_positions:
                # produce alt codon
                alt_seq = list(codon)
                # find index of pos in codon_positions
                idx = codon_positions.index(pos)
                alt_seq[idx] = alt if alt else alt_seq[idx]
                codon_to = ''.join(alt_seq).upper()
                aa_from = aa
                aa_to = str(Seq(codon_to).translate())
            else:
                codon_to = None
                aa_from = None
                aa_to = None
            rows.append({'time_iso': time_iso, 'time_year': time_year, 'ctx96': ctx,
                         'codon_from': codon_from, 'codon_to': codon_to,
                         'aa_from': aa_from, 'aa_to': aa_to,
                         'gene': gene, 'node_id': node_id, 'node_name': name,
                         'mut_type': 'nuc', 'pos': pos, 'ref': ref, 'alt': alt, 'raw': s})
        # process per-gene aa entries
        for gene, aa_list in muts.items():
            if gene == 'nuc':
                continue
            if not isinstance(aa_list, list):
                continue
            for s in aa_list:
                m = re.match(r'([A-Z\*])?(\d+)([A-Z\*])?', s)
                if not m:
                    continue
                aa_from = m.group(1)
                aa_pos = int(m.group(2))
                aa_to = m.group(3)
                # map aa_pos to genomic positions via gene_map
                positions = gene_map.get(gene)
                codon_from = None
                codon_to = None
                if positions:
                    start_idx = (aa_pos - 1) * 3
                    codon_positions = positions[start_idx:start_idx + 3]
                    if len(codon_positions) == 3:
                        codon_from = ''.join(ref_seq[p - 1] for p in codon_positions).upper()
                        # cannot know alt nt without nucleotide info; leave codon_to None
                rows.append({'time_iso': time_iso, 'time_year': time_year, 'ctx96': None,
                             'codon_from': codon_from, 'codon_to': codon_to,
                             'aa_from': aa_from, 'aa_to': aa_to,
                             'gene': gene, 'node_id': node_id, 'node_name': name,
                             'mut_type': 'aa', 'pos': aa_pos, 'ref': None, 'alt': None, 'raw': s})

    df = pl.DataFrame(rows)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df.write_parquet(out_path)
    print('Wrote', out_path)


def process_influenza(nwk_path, fasta_path, out_path):
    print('Processing Influenza consensus and tree (mutations if present)...')
    # build simple consensus
    consensus = build_influenza_consensus(fasta_path)
    # try to parse Newick for node labels/comments using Biopython
    try:
        from Bio import Phylo
        tree = Phylo.read(nwk_path, 'newick')
    except Exception:
        tree = None
    rows = []
    node_counter = 0
    if tree is not None:
        for clade in tree.find_clades(order='preorder'):
            node_id = f'NODE_{node_counter:06d}'
            node_counter += 1
            name = clade.name
            # attempt to extract mutations from comments
            raw = getattr(clade, 'comment', None)
            rows.append({'time_iso': None, 'time_year': None, 'ctx96': None,
                         'codon_from': None, 'codon_to': None,
                         'aa_from': None, 'aa_to': None,
                         'gene': None, 'node_id': node_id, 'node_name': name,
                         'mut_type': None, 'pos': None, 'ref': None, 'alt': None, 'raw': raw})
    else:
        # fallback: produce one consensus-only row
        rows.append({'time_iso': None, 'time_year': None, 'ctx96': None,
                     'codon_from': None, 'codon_to': None,
                     'aa_from': None, 'aa_to': None,
                     'gene': None, 'node_id': 'NODE_000000', 'node_name': None,
                     'mut_type': None, 'pos': None, 'ref': None, 'alt': None, 'raw': None})

    df = pl.DataFrame(rows)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df.write_parquet(out_path)
    print('Wrote', out_path)


def main():
    # Paths (workspace-relative)
    covid_json = 'Genomes/Covid/ncov_open_global_6m.json'
    covid_ref = 'Genomes/Covid/covidGenomeFasta.txt'
    covid_gtf = 'Genomes/Covid/MT291829_cds.gtf'
    covid_out = 'Results/Covid/node_mutations.parquet'

    influenza_nwk = 'Genomes/Influenza/nextstrain_seasonal-flu_h1n1pdm_ha_2y_timetree.nwk'
    influenza_fasta = 'Genomes/Influenza/sequences.fasta'
    influenza_out = 'Results/Influenza/node_mutations.parquet'

    process_covid(covid_json, covid_ref, covid_gtf, covid_out)
    process_influenza(influenza_nwk, influenza_fasta, influenza_out)


if __name__ == '__main__':
    main()
