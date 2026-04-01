"""codon_frequency.py – Convert CDS sequences to codon-count features.

For each gene_id, per-transcript codon counts are summed into a codon-count
dictionary. Codons that are the *only* codon for their amino acid are excluded
(currently `ATG` for M and `TGG` for W).

NumPy arrays are still available, with one element per retained codon, ordered
alphabetically by codon.

Typical use
-----------
    python codon_frequency.py HumanGenome/gene_id_to_sequence_transcripts.json

Outputs
-------
* Builds per-gene codon-count dictionaries.
* Optionally returns/saves a matrix where columns follow `CODON_ORDER`.
"""

import argparse
import json
import numpy as np
from collections import Counter, defaultdict
from pathlib import Path

# ---------------------------------------------------------------------------
# Standard genetic code (uppercase DNA)
# ---------------------------------------------------------------------------
GENETIC_CODE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# ---------------------------------------------------------------------------
# Keep only codons for amino acids with synonymous alternatives.
AA_TO_CODONS: dict[str, list[str]] = defaultdict(list)
for _codon, _aa in GENETIC_CODE.items():
    AA_TO_CODONS[_aa].append(_codon)

CODON_ORDER: list[str] = sorted(
    [
        codon
        for codon, aa in GENETIC_CODE.items()
        if len(AA_TO_CODONS[aa]) > 1
    ]
)
"""Filtered codon column order (excludes singleton-AA codons: ATG, TGG)."""

CODON_INDEX: dict[str, int] = {codon: i for i, codon in enumerate(CODON_ORDER)}
"""Reverse lookup: codon → column index in the frequency array."""

_VALID_BASES = frozenset("ACGT")


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def transcript_codon_counts(sequence: str, offset: int) -> Counter:
    """Count codons in *sequence* starting at *offset* (frame-corrected 0-based).

    Codons containing ambiguous bases (N, etc.) are skipped.
    """
    counter: Counter = Counter()
    seq = sequence.upper()
    for i in range(offset, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        if len(codon) == 3 and _VALID_BASES.issuperset(codon):
            counter[codon] += 1
    return counter


def gene_codon_count_dict(transcripts: dict) -> dict[str, int]:
    """Return filtered codon counts for one gene.

    Counts are summed across all transcripts. Each transcript's
    ``first_complete_codon_offset`` is used so reading begins on a full codon.

    Only codons from amino acids with multiple synonymous codons are kept.

    Parameters
    ----------
    transcripts:
        The inner dict for one gene_id from the JSON file, i.e.
        ``{transcript_id: {sequence, first_complete_codon_offset, ...}}``.
    """
    total: Counter = Counter()
    for record in transcripts.values():
        offset = record.get("first_complete_codon_offset", 0)
        total += transcript_codon_counts(record["sequence"], offset)

    return {codon: int(total[codon]) for codon in CODON_ORDER if total[codon] > 0}


def gene_codon_frequency_array(transcripts: dict) -> np.ndarray:
    """Return a filtered codon-count array for one gene.

    Array columns follow `CODON_ORDER`.
    """
    counts = gene_codon_count_dict(transcripts)

    arr = np.zeros(len(CODON_ORDER), dtype=np.int64)
    for codon, count in counts.items():
        idx = CODON_INDEX.get(codon)
        if idx is not None:
            arr[idx] = count
    return arr


def build_gene_codon_count_dict(gene_dict: dict) -> dict[str, dict[str, int]]:
    """Return `{gene_id: {codon: count}}` for all genes.

    Only codons in `CODON_ORDER` are included.
    """
    return {
        gene_id: gene_codon_count_dict(transcripts)
        for gene_id, transcripts in gene_dict.items()
    }


def build_codon_frequency_matrix(
    gene_dict: dict,
) -> tuple[list[str], np.ndarray]:
    """Compute codon frequency arrays for every gene in *gene_dict*.

    Returns
    -------
    gene_ids : list[str]
        Gene IDs in sorted order — the row labels for *matrix*.
    matrix : np.ndarray, shape (n_genes, n_codons)
        Row i corresponds to gene_ids[i]; column j to CODON_ORDER[j].
    """
    gene_ids = sorted(gene_dict.keys())
    matrix = np.zeros((len(gene_ids), len(CODON_ORDER)), dtype=np.int64)
    for i, gene_id in enumerate(gene_ids):
        matrix[i] = gene_codon_frequency_array(gene_dict[gene_id])
    return gene_ids, matrix


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert CDS JSON to per-gene codon counts (excluding singleton-AA codons)."
    )
    parser.add_argument(
        "json_file",
        nargs="?",
        default="HumanGenome/gene_id_to_sequence_transcripts.json",
        help="Path to the gene_id_to_transcripts JSON (default: %(default)s)",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Optional output path (.npz for matrix or .json for dict).",
    )
    parser.add_argument(
        "--show-codons",
        action="store_true",
        help="Print retained codon order and exit.",
    )
    args = parser.parse_args()

    if args.show_codons:
        for i, codon in enumerate(CODON_ORDER):
            aa = GENETIC_CODE[codon]
            print(f"{i:2d}  {codon}  {aa}")
        return

    print(f"Loading {args.json_file} …")
    with open(args.json_file, encoding="utf-8") as fh:
        gene_dict = json.load(fh)

    print(f"Computing codon counts for {len(gene_dict):,} genes …")
    gene_count_dict = build_gene_codon_count_dict(gene_dict)
    gene_ids, matrix = build_codon_frequency_matrix(gene_dict)

    print(f"Matrix shape: {matrix.shape}  (genes × {len(CODON_ORDER)} codons)")
    print(f"Total codons counted: {matrix.sum():,}")

    if args.output:
        out_path = Path(args.output)
        if out_path.suffix.lower() == ".json":
            with open(out_path, "w", encoding="utf-8") as fh:
                json.dump(gene_count_dict, fh, indent=2)
            print(f"Saved codon-count dict to {out_path}")
        else:
            np.savez_compressed(
                out_path,
                matrix=matrix,
                gene_ids=np.array(gene_ids),
                codon_order=np.array(CODON_ORDER),
            )
            print(f"Saved matrix to {out_path}")
    else:
        # Print a small summary for the first 5 genes
        print("\nSample output (first 5 genes, count of retained codons per gene):")
        for i in range(min(5, len(gene_ids))):
            print(f"  {gene_ids[i]}  total_codons={matrix[i].sum()}")


if __name__ == "__main__":
    main()
