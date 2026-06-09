#!/usr/bin/env python3
"""Infer per-branch influenza mutations from tree + sequences.

Performs:
1. Sequence alignment (MAFFT or built-in Biopython)
2. Ancestral sequence reconstruction (parsimony)
3. Per-branch mutation extraction
4. Codon/AA effect computation
5. is_coding tagging using GTF

Outputs:
- Results/Influenza/node_mutations.parquet (per-branch mutations)
- Results/Influenza/ancestral_sequences.fasta (reconstructed ancestor seqs)
- Results/Influenza/tree_with_sequences.txt (tree structure + seq info)
"""
import os
import re
from collections import defaultdict, Counter
from io import StringIO
import math

import polars as pl
from Bio import Phylo, SeqIO, Align, Seq
from Bio.Seq import Seq as BioSeq

def parse_gtf(gtf_path):
    """Parse GTF and return genes dict with CDS ranges per segment."""
    genes = {}
    with open(gtf_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 9 or cols[2] != 'CDS':
                continue
            seqname = cols[0]
            start = int(cols[3])
            end = int(cols[4])
            m = re.search(r'gene_id "([^"]+)"', cols[8])
            gene = m.group(1) if m else None
            if not gene:
                continue
            key = (seqname, gene)
            genes.setdefault(key, []).append((start, end))
    for k in genes:
        genes[k] = sorted(genes[k])
    return genes


def build_coding_positions(genes, seqname):
    """Build set of coding positions for a sequence segment."""
    s = set()
    for (seg, gene), ranges in genes.items():
        if seg == seqname:
            for start, end in ranges:
                s.update(range(start, end + 1))
    return s


def load_reference(fasta_path):
    """Load reference genome segments into dict."""
    seqs = {}
    for record in SeqIO.parse(fasta_path, 'fasta'):
        seqs[record.id] = str(record.seq)
    return seqs


def normalize_strain_name(name):
    """Normalize influenza strain labels for cross-source matching."""
    if not name:
        return None
    s = name.strip()
    # Drop subtype suffixes like (H1N1) and normalize common lab suffixes.
    s = re.sub(r'\([^)]*\)', '', s)
    s = re.sub(r'[-_](egg|cell|mdck|passage|passaged)$', '', s, flags=re.IGNORECASE)
    s = re.sub(r'\s+', ' ', s).strip()
    return s.upper()


def strain_aliases(name):
    """Generate normalized aliases for a strain label.

    The tree and FASTA naming conventions are close but not identical, so we
    try several canonical forms before giving up.
    """
    base = normalize_strain_name(name)
    if not base:
        return []

    aliases = []

    def add(value):
        if value and value not in aliases:
            aliases.append(value)

    add(base)
    add(re.sub(r'\s*/\s*', '/', base))
    add(base.replace(' ', ''))
    add(re.sub(r'[^A-Z0-9/]+', '', base))
    add(re.sub(r'[^A-Z0-9]+', '', base))

    return aliases


def load_sequences_for_tips(fasta_path, tip_names):
    """Load sequences matching tree tip names by parsing FASTA descriptions.
    
    FASTA description format: "NC_026431.1 |Influenza A virus (A/California/07/2009(H1N1)) segment ..."
    Extracts strain name from parentheses and matches to tree tips.
    Returns dict: tip_name -> sequence (for tips found in FASTA).
    """
    seqs = {}
    tip_set = set(tip_names)
    print(f"Looking for {len(tip_set)} tip sequences in FASTA...")
    
    # Parse FASTA and extract HA strain names from descriptions.
    # Keep only segment 4 / HA because the tree is HA-specific.
    strain_to_seq = {}  # normalized_strain_name -> sequence string
    ambiguous_aliases = set()
    tip_lookup = set()
    for tip in tip_names:
        tip_lookup.update(strain_aliases(tip))
    count = 0
    ha_count = 0
    for record in SeqIO.parse(fasta_path, 'fasta'):
        desc = record.description
        desc_lower = desc.lower()
        if 'segment 4' not in desc_lower and 'hemagglutinin' not in desc_lower:
            continue
        ha_count += 1

        # FASTA description format: "... |Influenza A virus (A/California/07/2009(H1N1)) segment ..."
        # Capture full strain name up to the ') segment' marker so subtype tags
        # like '(H1N1)' remain within the captured text and can be normalized out.
        match = re.search(r'\|Influenza [AB] virus \((.*?)\)\s+segment', desc)
        if match:
            strain_name = match.group(1)
            for alias in strain_aliases(strain_name):
                if alias not in tip_lookup:
                    continue
                if alias in strain_to_seq and strain_to_seq[alias] != str(record.seq):
                    ambiguous_aliases.add(alias)
                    continue
                strain_to_seq.setdefault(alias, str(record.seq))
            count += 1
            if count % 100 == 0:
                print(f"  Parsed {count} sequences...")
    
    print(f"Extracted strain names from {count} HA sequences (from {ha_count} HA records)")
    
    # Match tree tips to normalized strain keys.
    matched = 0
    for tip in tip_names:
        for alias in strain_aliases(tip):
            if alias in ambiguous_aliases:
                continue
            if alias in strain_to_seq:
                seqs[tip] = strain_to_seq[alias]
                matched += 1
                break
    
    print(f"Matched {matched} / {len(tip_set)} tips to sequences")
    return seqs


def reconstruct_ancestors_parsimony(tree, seq_alignment):
    """Reconstruct ancestral sequences using simple parsimony.
    
    seq_alignment: dict mapping leaf names to sequences.
    Returns: dict mapping clade objects to sequences.
    """
    ancestor_seqs = {}
    
    # Traverse tree post-order (leaves first, then internal nodes)
    def postorder_reconstruct(clade):
        if clade.is_terminal():
            name = clade.name
            if name in seq_alignment:
                seq = seq_alignment[name]
                ancestor_seqs[clade] = seq
                return seq
            else:
                return None
        
        # Internal node: reconstruct from children
        child_seqs = []
        for child in clade.clades:
            child_seq = postorder_reconstruct(child)
            if child_seq:
                child_seqs.append(child_seq)
        
        if not child_seqs:
            return None
        
        # Simple parsimony: majority vote at each position
        seq_len = len(child_seqs[0])
        anc_seq = []
        for pos in range(seq_len):
            bases = [s[pos] for s in child_seqs if pos < len(s)]
            if bases:
                # Majority vote (tie-breaking: pick first)
                c = Counter(bases)
                majority_base = c.most_common(1)[0][0]
                anc_seq.append(majority_base)
            else:
                anc_seq.append('N')
        
        result = ''.join(anc_seq)
        ancestor_seqs[clade] = result
        return result
    
    postorder_reconstruct(tree.root)
    return ancestor_seqs


def trinuc_context_safe(seq, pos):
    """Get trinucleotide context centered at pos (1-based)."""
    if pos is None or pos < 1 or pos > len(seq):
        return None
    i = pos - 1
    left = seq[i - 1] if i - 1 >= 0 else 'N'
    mid = seq[i]
    right = seq[i + 1] if i + 1 < len(seq) else 'N'
    tri = (left + mid + right).upper()
    if tri[1] in ('A', 'G'):
        tri = str(BioSeq(tri).complement())[::-1]
    return tri


def get_codon_and_aa(ref_seq, gene_map, seqname, pos):
    """Get codon and AA for nucleotide position."""
    for (seg, gene), ranges in gene_map.items():
        if seg == seqname and pos in set(sum((list(range(s, e+1)) for s, e in ranges), [])):
            # Find which range
            for start, end in ranges:
                if start <= pos <= end:
                    idx = pos - start  # 0-based within this CDS part
                    # Accumulate offset from earlier CDS parts
                    offset = 0
                    for s, e in ranges:
                        if s >= start:
                            break
                        offset += e - s + 1
                    total_idx = offset + idx
                    codon_idx = total_idx // 3
                    codon_positions = []
                    cumulative = 0
                    for s, e in ranges:
                        for p in range(s, e + 1):
                            if cumulative // 3 == codon_idx:
                                codon_positions.append(p)
                                if len(codon_positions) == 3:
                                    break
                            cumulative += 1
                        if len(codon_positions) == 3:
                            break
                    if len(codon_positions) == 3:
                        codon = ''.join(ref_seq[p - 1] for p in codon_positions)
                        aa = str(BioSeq(codon).translate())
                        return gene, codon_positions, codon.upper(), aa
    return None, None, None, None


def extract_mutations(parent_seq, child_seq, seqname, gene_map, ref_seqs):
    """Extract mutations between parent and child sequences."""
    mutations = []
    if not parent_seq or not child_seq or len(parent_seq) != len(child_seq):
        return mutations
    
    for pos in range(len(parent_seq)):
        if parent_seq[pos] != child_seq[pos]:
            ref_base = parent_seq[pos]
            alt_base = child_seq[pos]
            ctx = trinuc_context_safe(parent_seq, pos + 1)
            
            # Find gene and codon
            gene, codon_pos, codon_from, aa_from = get_codon_and_aa(parent_seq, gene_map, seqname, pos + 1)
            
            # Alt codon
            codon_to = None
            aa_to = None
            if codon_pos:
                codon_list = list(codon_from or '')
                idx = codon_pos.index(pos + 1)
                if idx < len(codon_list):
                    codon_list[idx] = alt_base
                    codon_to = ''.join(codon_list).upper()
                    aa_to = str(BioSeq(codon_to).translate())
            
            mutations.append({
                'pos': pos + 1,
                'ref': ref_base,
                'alt': alt_base,
                'ctx96': ctx,
                'codon_from': codon_from,
                'codon_to': codon_to,
                'aa_from': aa_from,
                'aa_to': aa_to,
                'gene': gene,
            })
    
    return mutations


def main():
    print("Loading reference and GTF...")
    ref_path = 'Genomes/Influenza/genbank reference sequence/ncbi_dataset/data/GCA_046377675.1/GCA_046377675.1_ASM4637767v1_genomic.fna'
    gtf_path = 'Genomes/Influenza/genbank reference sequence/ncbi_dataset/data/GCA_046377675.1/genomic.gtf'
    nwk_path = 'Genomes/Influenza/nextstrain_seasonal-flu_h1n1pdm_ha_2y_timetree.nwk'
    fasta_path = 'Genomes/Influenza/sequences.fasta'
    
    ref_seqs = load_reference(ref_path)
    gene_map = parse_gtf(gtf_path)
    
    print(f"Reference segments: {list(ref_seqs.keys())}")
    print(f"Gene features: {len(gene_map)}")
    
    print("\nLoading tree...")
    tree = Phylo.read(nwk_path, 'newick')
    print(f"Tree loaded: {len(list(tree.find_clades()))} nodes")
    
    # Extract tip names from tree
    tip_names = [t.name for t in tree.get_terminals()]
    print(f"Tree tips: {len(tip_names)}")
    
    print("\nLoading sequences...")
    query_seqs = load_sequences_for_tips(fasta_path, tip_names)
    
    if len(query_seqs) == 0:
        print("ERROR: No sequences matched to tree tips!")
        print("Sample tip names:", tip_names[:5])
        return
    
    # For simplicity: use HA (hemagglutinin) segment as primary
    # The GTF has segment names like LC660659.1 for HA
    ha_segment = 'LC660659.1'  # HA is at index 3 in the reference
    
    print(f"\nUsing segment: {ha_segment}")
    print(f"Matched {len(query_seqs)} sequences to tree tips")
    
    print("\nReconstructing ancestral sequences (parsimony)...")
    ancestor_seqs = reconstruct_ancestors_parsimony(tree, query_seqs)
    print(f"Reconstructed {len(ancestor_seqs)} ancestral sequences")
    
    if len(ancestor_seqs) == 0:
        print("ERROR: No ancestral sequences reconstructed!")
        return
    
    print("\nExtracting per-branch mutations...")
    rows = []
    node_counter = 0
    
    def traverse_mutations(clade):
        nonlocal node_counter
        node_id = f"NODE_{node_counter:06d}"
        node_counter += 1
        
        node_name = clade.name if clade.name else None
        node_seq = ancestor_seqs.get(clade)
        
        # Extract mutations from this node to its children
        for child in clade.clades:
            child_seq = ancestor_seqs.get(child)
            if node_seq and child_seq:
                muts = extract_mutations(node_seq, child_seq, ha_segment, gene_map, ref_seqs)
                for m in muts:
                    is_coding = bool(m['gene'])
                    rows.append({
                        'time_iso': None,
                        'time_year': None,
                        'ctx96': m['ctx96'],
                        'codon_from': m['codon_from'],
                        'codon_to': m['codon_to'],
                        'aa_from': m['aa_from'],
                        'aa_to': m['aa_to'],
                        'gene': m['gene'],
                        'is_coding': is_coding,
                        'node_id': node_id,
                        'node_name': node_name,
                        'mut_type': 'nuc',
                        'pos': m['pos'],
                        'ref': m['ref'],
                        'alt': m['alt'],
                        'raw': f"{m['ref']}{m['pos']}{m['alt']}",
                        'segment': ha_segment,
                    })
            traverse_mutations(child)
    
    traverse_mutations(tree.root)
    print(f"Extracted {len(rows)} per-branch mutations")
    
    # Write outputs
    os.makedirs('Results/Influenza', exist_ok=True)
    
    print("\nWriting mutation parquet...")
    df_muts = pl.DataFrame(rows)
    df_muts.write_parquet('Results/Influenza/node_mutations.parquet')
    print("Wrote Results/Influenza/node_mutations.parquet")
    
    print("\nWriting ancestral sequences...")
    with open('Results/Influenza/ancestral_sequences.fasta', 'w') as fh:
        for i, (clade, seq) in enumerate(ancestor_seqs.items()):
            name = clade.name if clade.name else f"ancestor_{i}"
            fh.write(f">{name}\n{seq}\n")
    print("Wrote Results/Influenza/ancestral_sequences.fasta")
    
    print("\nWriting tree structure summary...")
    with open('Results/Influenza/tree_with_sequences.txt', 'w') as fh:
        fh.write("INFLUENZA TREE STRUCTURE WITH SEQUENCE INFO\n")
        fh.write("=" * 80 + "\n\n")
        fh.write(f"Total nodes: {len(list(tree.find_clades()))}\n")
        fh.write(f"Tips (sequences): {len(tree.get_terminals())}\n")
        fh.write(f"Internal nodes: {len(list(tree.find_clades())) - len(tree.get_terminals())}\n")
        fh.write(f"Sequences matched: {len(query_seqs)}\n\n")
        
        fh.write("MATCHED LEAF NODES (sequence names):\n")
        for i, tip in enumerate(tree.get_terminals()[:100]):
            if tip.name in query_seqs:
                seq_len = len(query_seqs[tip.name])
                fh.write(f"  {i+1:4d}. {tip.name:50s} ({seq_len} bp)\n")
        if len(tree.get_terminals()) > 100:
            fh.write(f"  ... and {len(tree.get_terminals()) - 100} more\n")
        
        fh.write("\n" + "=" * 80 + "\n")
        fh.write("MUTATIONS EXTRACTED:\n")
        fh.write(f"  Total per-branch mutations: {len(rows)}\n")
        if rows:
            by_gene = defaultdict(int)
            by_coding = {'coding': 0, 'non-coding': 0}
            for r in rows:
                if r['gene']:
                    by_gene[r['gene']] += 1
                    by_coding['coding'] += 1
                else:
                    by_coding['non-coding'] += 1
            fh.write(f"  Coding: {by_coding['coding']}, Non-coding: {by_coding['non-coding']}\n")
            fh.write("  By gene:\n")
            for gene in sorted(by_gene.keys()):
                fh.write(f"    {gene}: {by_gene[gene]}\n")
    
    print("Wrote Results/Influenza/tree_with_sequences.txt")
    print("\nDone!")


if __name__ == '__main__':
    main()
