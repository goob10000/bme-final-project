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


def gene_codon_count_dict(transcripts: dict) -> tuple[dict[str, int], float]:
    """Return filtered codon counts for one gene.
    Each transcript's ``first_complete_codon_offset`` is used so reading begins on a full codon.

    Only codons from amino acids with multiple synonymous codons are kept.

    Parameters
    ----------
    transcripts:
        The inner dict for one gene_id from the JSON file, i.e.
        ``{transcript_id: {sequence, first_complete_codon_offset, ...}}``.
    """
    total: Counter = Counter()
    numTranscripts = len(transcripts.values()) ## Average transcripts if there are more than one.
    cdsLength = 0
    for record in transcripts.values():
        offset = record.get("first_complete_codon_offset", 0)
        total += transcript_codon_counts(record["sequence"], offset)
        cdsLength += record["cds_length"]

    return ({codon: int(total[codon]/numTranscripts) for codon in CODON_ORDER if total[codon] > 0}, cdsLength/numTranscripts)


def gene_codon_frequency_array(transcripts: dict) -> tuple[np.ndarray, float]:
    """Return a filtered codon-proportion array for one gene.

    For each amino acid, codon values are normalized to proportions that sum to 1
    across that amino acid's synonymous codons present in the gene.
    Array columns follow `CODON_ORDER`.
    """
    counts, length = gene_codon_count_dict(transcripts)

    arr = np.zeros(len(CODON_ORDER), dtype=np.float64)
    aa_totals: dict[str, float] = defaultdict(float)
    for codon, count in counts.items():
        aa_totals[GENETIC_CODE[codon]] += float(count)

    for codon, count in counts.items():
        idx = CODON_INDEX.get(codon)
        aa_total = aa_totals[GENETIC_CODE[codon]]
        if idx is not None and aa_total > 0:
            arr[idx] = float(count) / aa_total
    return arr, length


def build_gene_codon_count_dict(gene_dict: dict) -> dict[str, dict[str, int]]:
    """Return `{gene_id: {codon: count}}` for all genes.

    Only codons in `CODON_ORDER` are included.

    Returns a dict that takes the gene ID and returns a dict [codon → count] for that gene.
    """
    return {
        gene_id: gene_codon_count_dict(transcripts)[0] for gene_id, transcripts in gene_dict.items()
    }


def build_codon_frequency_matrix(
    gene_dict: dict,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Compute codon frequency arrays for every gene in *gene_dict*.

    Returns
    -------
    gene_ids : list[str]
        Gene IDs in sorted order — the row labels for *matrix*.
    matrix : np.ndarray, shape (n_genes, n_codons)
        Row i corresponds to gene_ids[i]; column j to CODON_ORDER[j].
    """
    gene_ids = sorted(gene_dict.keys())
    matrix = np.zeros((len(gene_ids), len(CODON_ORDER)), dtype=np.float64)
    counts = np.zeros(len(gene_ids), dtype=np.float64)
    for i, gene_id in enumerate(gene_ids):
        matrix[i], counts[i] = gene_codon_frequency_array(gene_dict[gene_id])
    return gene_ids, matrix, counts

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert CDS JSON to per-gene codon counts (excluding singleton-AA codons)."
    )
    parser.add_argument(
        "--json_file",
        nargs="?",
        help="Path to the gene_id_to_transcripts JSON (default: %(default)s)",
    )
    parser.add_argument(
        "-o", "--output",
        help="Required output path (.npz for matrix or .json for dict).",
        required=True,
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
    gene_count_dict = build_gene_codon_count_dict(gene_dict) # dict[gene_id → dict[codon → count]]
    gene_ids, matrix, counts = build_codon_frequency_matrix(gene_dict)

    print(f"Matrix shape: {matrix.shape}  (genes × {len(CODON_ORDER)} codons)")
    print(f"Total codons counted: {matrix.sum():,}")
    print(f"Average codons per gene: {counts.mean():.1f}")

    out_path = Path(args.output)
    np.savez_compressed(
        out_path,
        matrix=matrix,
        gene_ids=np.array(gene_ids),
        codon_order=np.array(CODON_ORDER),
        counts=counts,
    )
    print(f"Saved matrix to {out_path}")

if __name__ == "__main__":
    main()
