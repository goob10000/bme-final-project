import argparse
import csv
import gzip
from pathlib import Path


def parse_geo_tsv_row(line: str):
    reader = csv.reader([line.rstrip("\n")], delimiter="\t", quotechar='"')
    return next(reader)


def parse_series_metadata(series_matrix_path: Path):
    geo_accessions = []
    titles = []
    disease_states = []

    with series_matrix_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("!Sample_geo_accession"):
                geo_accessions = parse_geo_tsv_row(line)[1:]
            elif line.startswith("!Sample_title"):
                titles = parse_geo_tsv_row(line)[1:]
            elif line.startswith("!Sample_characteristics_ch1") and "disease state:" in line.lower():
                vals = parse_geo_tsv_row(line)[1:]
                disease_states = [v.split(":", 1)[1].strip() if ":" in v else v.strip() for v in vals]

    if not geo_accessions:
        raise ValueError(f"Could not find !Sample_geo_accession in {series_matrix_path}")

    if len(titles) != len(geo_accessions):
        titles = [""] * len(geo_accessions)
    if len(disease_states) != len(geo_accessions):
        disease_states = [""] * len(geo_accessions)

    out = {}
    for gsm, title, disease in zip(geo_accessions, titles, disease_states):
        out[gsm] = {"title": title, "disease": disease}
    return out


def sample_label_from_title(title: str, fallback: str):
    # GEO titles look like: "Lung airway epithelium, H1"
    if "," in title:
        tail = title.split(",")[-1].strip()
        if tail:
            return tail.replace(" ", "_")
    return fallback


def open_text_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def resolve_sample_file(human2_dir: Path, prefix: str, stem: str):
    # Prefer extracted plain-text files; fall back to gzipped files if needed.
    # Support both naming styles:
    #   <prefix>.matrix.mtx / <prefix>.features.tsv
    #   <prefix>_matrix.mtx / <prefix>_features.tsv
    for suffix in ["", ".gz"]:
        candidates = [
            human2_dir / f"{prefix}.{stem}{suffix}",
            human2_dir / f"{prefix}_{stem}{suffix}",
        ]
        for candidate in candidates:
            if candidate.exists():
                return candidate
    return None


def read_features(features_path: Path):
    gene_ids = []
    with open_text_maybe_gzip(features_path) as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            gene_ids.append(fields[0])
    return gene_ids


def sum_gene_counts_from_mtx(matrix_path: Path, gene_count: int):
    sums = [0] * gene_count

    with open_text_maybe_gzip(matrix_path) as handle:
        dims_seen = False
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("%"):
                continue

            if not dims_seen:
                # Matrix Market dimensions line: n_genes n_barcodes nnz
                dims = line.split()
                if len(dims) != 3:
                    raise ValueError(f"Invalid dimensions line in {matrix_path}: {line}")
                n_rows = int(dims[0])
                if n_rows != gene_count:
                    raise ValueError(
                        f"Gene count mismatch in {matrix_path}: features has {gene_count}, matrix has {n_rows}"
                    )
                dims_seen = True
                continue

            parts = line.split()
            if len(parts) != 3:
                continue
            gene_idx = int(parts[0]) - 1
            val = int(float(parts[2]))
            sums[gene_idx] += val

    return sums


def find_sample_prefixes(human2_dir: Path):
    prefixes = set()

    for p in human2_dir.glob("GSM*.matrix.mtx"):
        prefixes.add(p.name.replace(".matrix.mtx", ""))
    for p in human2_dir.glob("GSM*.matrix.mtx.gz"):
        prefixes.add(p.name.replace(".matrix.mtx.gz", ""))
    for p in human2_dir.glob("GSM*_matrix.mtx"):
        prefixes.add(p.name.replace("_matrix.mtx", ""))
    for p in human2_dir.glob("GSM*_matrix.mtx.gz"):
        prefixes.add(p.name.replace("_matrix.mtx.gz", ""))

    return sorted(prefixes)


def build_sample_columns(human2_dir: Path, sample_prefixes, metadata_by_gsm):
    ordered_gene_ids = None
    columns = []

    for prefix in sample_prefixes:
        gsm = prefix.split("_")[0]
        meta = metadata_by_gsm.get(gsm, {})
        title = meta.get("title", "")
        label = sample_label_from_title(title, prefix)

        features_path = resolve_sample_file(human2_dir, prefix, "features.tsv")
        matrix_path = resolve_sample_file(human2_dir, prefix, "matrix.mtx")

        if features_path is None or matrix_path is None:
            continue

        gene_ids = read_features(features_path)
        sums = sum_gene_counts_from_mtx(matrix_path, len(gene_ids))

        if ordered_gene_ids is None:
            ordered_gene_ids = gene_ids
            gene_index = {gid: i for i, gid in enumerate(ordered_gene_ids)}
        else:
            if gene_ids != ordered_gene_ids:
                # Align if gene ordering differs between samples.
                aligned = [0] * len(ordered_gene_ids)
                for gid, val in zip(gene_ids, sums):
                    idx = gene_index.get(gid)
                    if idx is not None:
                        aligned[idx] = val
                sums = aligned

        columns.append((label, sums))

    if ordered_gene_ids is None or not columns:
        raise ValueError("No valid samples found to process.")

    return ordered_gene_ids, columns


def write_counts_csv(output_path: Path, gene_ids, columns):
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Gene_ID"] + [name for name, _ in columns])

        for row_idx, gid in enumerate(gene_ids):
            row = [gid]
            for _, col in columns:
                row.append(col[row_idx])
            writer.writerow(row)


def read_gene_lengths(gene_lengths_path: Path):
    gene_lengths = {}
    with gene_lengths_path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or len(row) < 2:
                continue
            gene_id = row[0].strip()
            if not gene_id or gene_id.lower() == "gene_id":
                continue
            try:
                gene_lengths[gene_id] = float(row[1])
            except ValueError:
                continue
    return gene_lengths


def normalize_columns(gene_ids, columns, norm, gene_lengths=None):
    if norm == "raw":
        return columns

    normalized = []

    for sample_name, counts in columns:
        total_counts = sum(counts)

        if norm == "cpm":
            if total_counts == 0:
                values = [0.0] * len(counts)
            else:
                scale = 1_000_000.0 / float(total_counts)
                values = [c * scale for c in counts]
            normalized.append((sample_name, values))
            continue

        if norm == "tpm":
            if gene_lengths is None:
                raise ValueError("TPM requested but no gene lengths were supplied.")

            rpk = []
            for gid, c in zip(gene_ids, counts):
                glen = gene_lengths.get(gid)
                if glen is None or glen <= 0:
                    rpk.append(0.0)
                else:
                    rpk.append(float(c) / (glen / 1000.0))

            rpk_sum = sum(rpk)
            if rpk_sum == 0:
                values = [0.0] * len(rpk)
            else:
                scale = 1_000_000.0 / rpk_sum
                values = [v * scale for v in rpk]

            normalized.append((sample_name, values))
            continue

        raise ValueError(f"Unsupported normalization: {norm}")

    return normalized


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct a Human1-like raw counts CSV from GEO GSE275264 10x files "
            "(gene rows, donor columns)."
        )
    )
    parser.add_argument(
        "--human2-dir",
        default="Expression/Human2",
        help="Directory with GSM*.features.tsv/.gz and GSM*.matrix.mtx/.gz files.",
    )
    parser.add_argument(
        "--needs-series-matrix",
        action="store_true",
        help="Whether to require the GEO series matrix for sample metadata. If not set, will attempt to process any sample files found in --human2-dir regardless of whether they have metadata entries.",
    )
    parser.add_argument(
        "--series-matrix",
        default="Expression/Human2/GSE275264_series_matrix.txt",
        help="Path to unzipped GEO series matrix TXT.",
    )
    parser.add_argument(
        "--group",
        choices=["healthy", "ipf", "all"],
        default="healthy",
        help=(
            "Which samples to include in output columns when --needs-series-matrix is set "
            "(healthy, ipf, or all). Ignored otherwise."
        ),
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output CSV path. Defaults to Expression/Human2/GSE275264_raw_counts_<group>.csv",
    )
    parser.add_argument(
        "--norm",
        choices=["raw", "cpm", "tpm"],
        default="raw",
        help="Normalization to apply per sample column before writing output.",
    )
    parser.add_argument(
        "--gene-lengths",
        default=None,
        help=(
            "Tab-delimited file with gene length in bases: gene_id<TAB>length_bp. "
            "Required when --norm tpm is used."
        ),
    )
    args = parser.parse_args()

    human2_dir = Path(args.human2_dir)
    metadata_by_gsm = {}
    if args.needs_series_matrix and not Path(args.series_matrix).exists():
        raise ValueError(f"--needs-series-matrix is set but {args.series_matrix} does not exist.")
    if args.needs_series_matrix:
        series_matrix_path = Path(args.series_matrix)
        metadata_by_gsm = parse_series_metadata(series_matrix_path)

    effective_group = args.group if args.needs_series_matrix else "all"

    out_path = (
        Path(args.output)
        if args.output
        else human2_dir / f"GSE275264_{args.norm}_{effective_group}.csv"
    )

    all_prefixes = find_sample_prefixes(human2_dir)

    if not args.needs_series_matrix:
        # Without series metadata, include every sample detected in the directory.
        selected = list(all_prefixes)
        if args.group != "all":
            print(
                "Note: --group is ignored unless --needs-series-matrix is set; processing all samples found in --human2-dir."
            )
    else:
        selected = []
        for prefix in all_prefixes:
            gsm = prefix.split("_")[0]
            disease = metadata_by_gsm.get(gsm, {}).get("disease", "").strip().lower()

            if args.group == "all":
                selected.append(prefix)
            elif args.group == "healthy" and disease == "healthy":
                selected.append(prefix)
            elif args.group == "ipf" and disease == "ipf":
                selected.append(prefix)

    gene_ids, columns = build_sample_columns(human2_dir, selected, metadata_by_gsm)

    gene_lengths = None
    if args.norm == "tpm":
        if not args.gene_lengths:
            raise ValueError("--norm tpm requires --gene-lengths.")
        gene_lengths = read_gene_lengths(Path(args.gene_lengths))

    out_columns = normalize_columns(gene_ids, columns, args.norm, gene_lengths)
    write_counts_csv(out_path, gene_ids, out_columns)

    print(f"Wrote {out_path}")
    print(f"Genes: {len(gene_ids)}")
    print(f"Samples: {len(out_columns)}")
    print(f"Normalization: {args.norm}")
    print("Sample columns:", ", ".join(name for name, _ in out_columns))


if __name__ == "__main__":
    main()
