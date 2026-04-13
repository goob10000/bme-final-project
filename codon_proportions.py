"""topn_codon_proportions.py - Expression-weighted codon usage for top-N genes.

Builds a tidy table where each row is one codon contribution for one amino acid
at a specific top-N cutoff (and sample). This is designed for plotting curves of
codon proportion vs number of top-expressed genes included.

Core idea
---------
1) Build per-gene codon counts from transcript CDS JSON using codon_frequency.py.
2) Rank genes by expression for a sample.
3) For top N genes, weight each gene's codon counts by its expression.
4) Within each amino acid, convert codon weighted counts to proportions.

Output columns
--------------
sample, top_n, aa, codon, weighted_codon_count, aa_total_weighted, codon_proportion
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

from codon_frequency import GENETIC_CODE, build_gene_codon_count_dict


def normalize_gene_id(gene_id: str) -> str:
    """Return stable ENSG id without version suffix."""
    return gene_id.split(".", 1)[0]


def load_expression_table(csv_path: str) -> tuple[list[dict[str, str]], list[str]]:
    """Load raw counts table.

    Returns
    -------
    rows:
        Raw rows as dicts.
    sample_columns:
        All expression sample columns (everything except Gene_ID).
    """
    with open(csv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None or "Gene_ID" not in reader.fieldnames:
            raise ValueError("CSV must contain a Gene_ID column")
        rows = list(reader)
        sample_columns = [c for c in reader.fieldnames if c != "Gene_ID"]
    return rows, sample_columns


def sorted_genes_for_sample(
    expression_rows: list[dict[str, str]],
    sample: str,
    available_genes: set[str],
) -> list[tuple[str, float]]:
    """Return (gene_id, expression) sorted descending for a sample.

    Genes missing in codon dictionary are dropped.
    """
    ranked: list[tuple[str, float]] = []
    for row in expression_rows:
        gid = normalize_gene_id(row["Gene_ID"])
        if gid not in available_genes:
            continue
        raw = row.get(sample, "0")
        try:
            expr = float(raw)
        except (TypeError, ValueError):
            expr = 0.0
        ranked.append((gid, expr))

    ranked.sort(key=lambda x: x[1], reverse=True)
    return ranked


def codon_proportions_for_top_n(
    ranked_genes: list[tuple[str, float]],
    gene_codon_counts: dict[str, dict[str, int]],
    top_n: int,
    sample: str,
) -> list[dict[str, Any]]:
    """Build a plotting-ready long table for one sample and one top-N cutoff."""
    selected = ranked_genes[:top_n]

    # Weighted codon totals across selected genes.
    aa_codon_weighted: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))

    for gene_id, expr in selected:
        if expr <= 0:
            continue
        codon_counts = gene_codon_counts.get(gene_id)
        if not codon_counts:
            continue

        for codon, count in codon_counts.items():
            aa = GENETIC_CODE[codon]
            aa_codon_weighted[aa][codon] += expr * float(count)

    rows: list[dict[str, Any]] = []
    for aa in sorted(aa_codon_weighted.keys()):
        codon_map = aa_codon_weighted[aa]
        aa_total = sum(codon_map.values())
        if aa_total <= 0:
            continue

        for codon in sorted(codon_map.keys()):
            weighted = codon_map[codon]
            rows.append(
                {
                    "sample": sample,
                    "top_n": int(top_n),
                    "aa": aa,
                    "codon": codon,
                    "weighted_codon_count": weighted,
                    "aa_total_weighted": aa_total,
                    "codon_proportion": weighted / aa_total,
                }
            )

    return rows


def build_top_n_series(
    expression_rows: list[dict[str, str]],
    gene_codon_counts: dict[str, dict[str, int]],
    sample: str,
    n_values: list[int],
) -> list[dict[str, Any]]:
    """Compute codon proportion table for multiple top-N cutoffs."""
    available = set(gene_codon_counts.keys())
    ranked = sorted_genes_for_sample(expression_rows, sample, available)

    all_rows: list[dict[str, Any]] = []
    for n in n_values:
        if n <= 0:
            continue
        all_rows.extend(
            codon_proportions_for_top_n(
                ranked_genes=ranked,
                gene_codon_counts=gene_codon_counts,
                top_n=min(n, len(ranked)),
                sample=sample,
            )
        )
    return all_rows


def try_to_polars(rows: list[dict[str, Any]]):
    """Return a Polars DataFrame if Polars is available, else None."""
    try:
        import polars as pl  # type: ignore
    except ModuleNotFoundError:
        return None

    if not rows:
        return pl.DataFrame(
            {
                "sample": [],
                "top_n": [],
                "aa": [],
                "codon": [],
                "weighted_codon_count": [],
                "aa_total_weighted": [],
                "codon_proportion": [],
            }
        )

    df = pl.DataFrame(rows)
    return df.sort(["sample", "aa", "codon", "top_n"])


def write_table(rows: list[dict[str, Any]], output_path: str) -> None:
    """Write long table to JSON or CSV, and Parquet if Polars is available."""
    out = Path(output_path)
    suffix = out.suffix.lower()

    if suffix == ".json":
        with open(out, "w", encoding="utf-8") as fh:
            json.dump(rows, fh, indent=2)
        return

    if suffix == ".csv":
        fieldnames = [
            "sample",
            "top_n",
            "aa",
            "codon",
            "weighted_codon_count",
            "aa_total_weighted",
            "codon_proportion",
        ]
        with open(out, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        return

    if suffix == ".parquet":
        pl_df = try_to_polars(rows)
        if pl_df is None:
            raise RuntimeError(
                "Polars is not installed. Install it or use .csv/.json output."
            )
        pl_df.write_parquet(out)
        return

    raise ValueError("Output must end with .csv, .json, or .parquet")


def aa_to_filename_token(aa: str) -> str:
    """Map amino-acid code to a safe filename token."""
    if aa == "*":
        return "STOP"
    return aa

def parse_n_values(top_n: int, top_n_series: str | None) -> list[int]:
    """Return sorted unique N values."""
    values: set[int] = {int(top_n)}
    if top_n_series:
        for token in top_n_series.split(","):
            token = token.strip()
            if not token:
                continue
            values.add(int(token))
    return sorted(v for v in values if v > 0)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Expression-weighted codon proportion table for top-N genes."
    )
    parser.add_argument(
        "--transcript-json",
        help="Path to transcript-level CDS JSON.",
    )
    parser.add_argument(
        "--counts-csv",
        help="Path to expression raw counts CSV.",
    )
    parser.add_argument(
        "--sample",
        help="Sample column to use.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output table path (.csv, .json, or .parquet).",
    )
    args = parser.parse_args()

    print(f"Loading transcript CDS JSON: {args.transcript_json}")
    with open(args.transcript_json, encoding="utf-8") as fh:
        transcript_gene_dict = json.load(fh)

    print("Building per-gene codon counts from CDS data...")
    gene_codon_counts = build_gene_codon_count_dict(transcript_gene_dict)

    print(f"Loading expression counts: {args.counts_csv}")
    expression_rows, sample_cols = load_expression_table(args.counts_csv)

    sample = args.sample or sample_cols[0]
    if sample not in sample_cols:
        raise ValueError(
            f"Sample '{sample}' not in CSV columns. Available: {', '.join(sample_cols)}"
        )

    top_n = args.top_n
    print(f"Computing codon proportions for sample={sample}, top_n={top_n}")

    rows = build_top_n_series(
        expression_rows=expression_rows,
        gene_codon_counts=gene_codon_counts,
        sample=sample,
        n_values=n_values,
    )

    print(f"Writing table with {len(rows):,} rows to {args.output}")
    write_table(rows, args.output)

    print("Done.")


if __name__ == "__main__":
    main()
