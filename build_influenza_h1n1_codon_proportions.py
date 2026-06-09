import csv
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np

from codon_frequency import CODON_ORDER


SEGMENT_ORDER = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]
SEGMENT_FROM_NUMBER = {
    "1": "PB2",
    "2": "PB1",
    "3": "PA",
    "4": "HA",
    "5": "NP",
    "6": "NA",
    "7": "M",
    "8": "NS",
}
SEGMENT_INDEX = {segment: i for i, segment in enumerate(SEGMENT_ORDER)}


def parse_accession_from_header(header: str) -> str:
    token = header[1:].split(" ", 1)[0].strip()
    return token.split(".", 1)[0]


def parse_date(raw: str) -> tuple[str, datetime] | tuple[None, None]:
    raw = (raw or "").strip()
    if not raw:
        return None, None
    try:
        d = datetime.strptime(raw, "%Y-%m-%d")
        return d.strftime("%Y-%m-%d"), d
    except ValueError:
        pass
    if len(raw) == 4 and raw.isdigit():
        d = datetime(int(raw), 1, 1)
        return d.strftime("%Y-%m-%d"), d
    return None, None


def segment_name(row: dict[str, str]) -> str | None:
    raw_segment = (row.get("Segment") or "").strip()
    if raw_segment in SEGMENT_FROM_NUMBER:
        return SEGMENT_FROM_NUMBER[raw_segment]
    upper = raw_segment.upper()
    for segment in SEGMENT_ORDER:
        if segment in upper:
            return segment
    return None


def codon_proportions(seq: str) -> np.ndarray:
    counts = {codon: 0.0 for codon in [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"]}
    total_codons = 0
    sequence = seq.upper()
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]
        if all(base in "ACGT" for base in codon):
            counts[codon] += 1.0
            total_codons += 1
    if total_codons == 0:
        return np.zeros(len(CODON_ORDER), dtype=np.float64)
    return np.array([counts[codon] / total_codons for codon in CODON_ORDER], dtype=np.float64)


def load_h1n1_metadata(csv_path: Path) -> tuple[dict[str, list[dict[str, str]]], set[str], list[str]]:
    by_accession: dict[str, list[dict[str, str]]] = defaultdict(list)
    wanted_accessions: set[str] = set()
    fieldnames: list[str] = []
    with csv_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []
        for row in reader:
            if (row.get("Genotype") or "").strip().upper() != "H1N1":
                continue

            accession = (row.get("Accession") or "").strip().split(".", 1)[0]
            if not accession:
                continue

            parsed_date, parsed_dt = parse_date(row.get("Collection_Date") or "")
            if parsed_date is None:
                continue

            seg = segment_name(row)
            if seg is None:
                continue

            row = dict(row)
            row["_accession"] = accession
            row["_segment"] = seg
            row["_parsed_date"] = parsed_date
            row["_parsed_ord"] = str(parsed_dt.toordinal())
            by_accession[accession].append(row)
            wanted_accessions.add(accession)
    return by_accession, wanted_accessions, fieldnames


def load_sequences(fasta_path: Path, wanted_accessions: set[str]) -> dict[str, str]:
    sequences: dict[str, str] = {}
    header = None
    seq_chunks: list[str] = []

    def commit(current_header: str | None, chunks: list[str]) -> None:
        if current_header is None:
            return
        accession = parse_accession_from_header(current_header)
        if accession in wanted_accessions:
            sequences[accession] = "".join(chunks).upper()

    with fasta_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                commit(header, seq_chunks)
                header = line
                seq_chunks = []
            else:
                seq_chunks.append(line)
        commit(header, seq_chunks)

    return sequences


def primary_record_key(row: dict[str, str]) -> str:
    isolate = (row.get("Isolate") or "").strip()
    if isolate:
        return isolate
    assembly = (row.get("Assembly") or "").strip()
    if assembly:
        return assembly
    return row["_accession"]


def group_rows_by_strain(rows_by_accession: dict[str, list[dict[str, str]]]) -> dict[tuple[str, str], dict[str, dict[str, str]]]:
    grouped: dict[tuple[str, str], dict[str, dict[str, str]]] = {}
    for rows in rows_by_accession.values():
        for row in rows:
            key = (primary_record_key(row), row["_parsed_date"])
            if key not in grouped:
                grouped[key] = {}
            segment = row["_segment"]
            existing = grouped[key].get(segment)
            if existing is None:
                grouped[key][segment] = row
                continue
            old_len = int(existing.get("Length") or 0)
            new_len = int(row.get("Length") or 0)
            if new_len > old_len:
                grouped[key][segment] = row
    return grouped


def representative_field(segment_rows: dict[str, dict[str, str]], field: str) -> str:
    for segment in SEGMENT_ORDER:
        row = segment_rows.get(segment)
        if row is None:
            continue
        value = (row.get(field) or "").strip()
        if value:
            return value
    return ""


def build_matrix(
    grouped_rows: dict[tuple[str, str], dict[str, dict[str, str]]],
    sequences: dict[str, str],
    original_fields: list[str],
) -> tuple[np.ndarray, list[dict[str, str]]]:
    sorted_keys = sorted(grouped_rows.keys(), key=lambda x: (x[1], x[0]))
    matrix = np.zeros((len(sorted_keys), len(SEGMENT_ORDER), len(CODON_ORDER)), dtype=np.float64)
    metadata_rows: list[dict[str, str]] = []

    for row_index, key in enumerate(sorted_keys):
        strain_key, parsed_date = key
        segment_rows = grouped_rows[key]
        observed_segments = 0
        assembly = representative_field(segment_rows, "Assembly")

        for segment in SEGMENT_ORDER:
            row = segment_rows.get(segment)
            if row is None:
                continue
            seq = sequences.get(row["_accession"])
            if not seq:
                continue
            matrix[row_index, SEGMENT_INDEX[segment], :] = codon_proportions(seq)
            observed_segments += 1

        rowdict: dict[str, str] = {
            "row_index": str(row_index),
            "strain_key": strain_key,
        }
        # preserve original CSV columns by taking a representative value across segments
        for field in original_fields:
            # skip internal helper fields if present
            if field in ("_accession", "_segment", "_parsed_date", "_parsed_ord"):
                continue
            rowdict[field] = representative_field(segment_rows, field)

        rowdict.update(
            {
                "parsed_date": parsed_date,
                "date_ordinal": str(datetime.strptime(parsed_date, "%Y-%m-%d").toordinal()),
                "segments_present": str(observed_segments),
                "segment_list": ";".join(sorted(segment_rows.keys())),
            }
        )
        metadata_rows.append(rowdict)

    return matrix, metadata_rows


def save_metadata(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    root = Path(__file__).resolve().parent
    csv_path = root / "Genomes" / "Influenza" / "sequences.csv"
    fasta_path = root / "Genomes" / "Influenza" / "sequences.fasta"
    out_dir = root / "Results2" / "Influenza"
    out_dir.mkdir(parents=True, exist_ok=True)

    by_accession, wanted_accessions, original_fields = load_h1n1_metadata(csv_path)
    sequences = load_sequences(fasta_path, wanted_accessions)
    grouped = group_rows_by_strain(by_accession)
    matrix, metadata_rows = build_matrix(grouped, sequences, original_fields)

    np.save(out_dir / "codon_proportions.npy", matrix)
    # build output columns: include row_index, strain_key, all original fields, then parsed_date etc.
    columns = ["row_index", "strain_key"] + [f for f in original_fields if f not in (None, "") and f not in ("_accession", "_segment", "_parsed_date", "_parsed_ord")] + ["parsed_date", "date_ordinal", "segments_present", "segment_list"]
    save_metadata(out_dir / "codon_proportions_metadata.csv", metadata_rows, columns)

    print(f"H1N1 accessions in metadata: {len(wanted_accessions)}")
    print(f"H1N1 strain/date rows: {matrix.shape[0]}")
    print(f"Saved matrix shape: {matrix.shape} to {out_dir / 'codon_proportions.npy'}")
    if metadata_rows:
        covered = [int(row["segments_present"]) for row in metadata_rows]
        print(f"Segments present per row: min={min(covered)}, max={max(covered)}, avg={sum(covered)/len(covered):.2f}")


if __name__ == "__main__":
    main()