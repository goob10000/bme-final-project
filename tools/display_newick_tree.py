#!/usr/bin/env python3
"""Show Newick tree structure - compact ASCII format."""
from Bio import Phylo
from io import StringIO

nwk_path = 'Genomes/Influenza/nextstrain_seasonal-flu_h1n1pdm_ha_2y_timetree.nwk'

with open(nwk_path) as fh:
    newick_str = fh.read()

tree = Phylo.read(StringIO(newick_str), 'newick')

print("INFLUENZA H1N1 TREE STRUCTURE")
print("=" * 80)
print()

# Summary stats
total_clades = len(list(tree.find_clades()))
tips = len(tree.get_terminals())
internal = total_clades - tips
print(f"Total nodes: {total_clades}")
print(f"Tips (leaf nodes): {tips}")
print(f"Internal nodes: {internal}")
print()

# Show tree structure with ASCII art (first few levels)
def show_clade(clade, depth=0, max_depth=6, count=[0]):
    if count[0] >= 100:
        return
    if depth > max_depth:
        return
    
    count[0] += 1
    indent = "  " * depth
    node_type = "LEAF" if clade.is_terminal() else "NODE"
    name_str = clade.name if clade.name else f"internal_{count[0]}"
    branch_len = f" ({clade.branch_length:.4f})" if clade.branch_length else ""
    
    # Limit long names
    if len(name_str) > 50:
        name_str = name_str[:47] + "..."
    
    print(f"{indent}├─ {node_type}: {name_str}{branch_len}")
    
    for child in clade.clades:
        show_clade(child, depth + 1, max_depth, count)

print("TREE HIERARCHY (first 100 nodes, up to 6 levels deep):")
print()
show_clade(tree.root)

print()
print("=" * 80)
print("SAMPLE LEAF NAMES (first 20):")
print()
for i, tip in enumerate(tree.get_terminals()[:20]):
    print(f"  {i+1:2d}. {tip.name}")

print()
print("Tree branches have comments/annotations:", any(c.comment for c in tree.find_clades()))
