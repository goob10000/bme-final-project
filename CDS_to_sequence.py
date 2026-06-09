from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path


GENE_ID_RE = re.compile(r'gene_id\s+(?:"|"")([^"]+)(?:"|"")')
TRANSCRIPT_ID_RE = re.compile(r'transcript_id\s+(?:"|"")([^"]+)(?:"|"")')
TAG_RE = re.compile(r'tag\s+(?:"|"")([^"]+)(?:"|"")')


@dataclass(frozen=True)
class FastaRecordMeta:
	start_offset: int
	line_bases: int
	line_width: int
	length: int


@dataclass(frozen=True)
class CDSPart:
	chrom: str
	start: int
	end: int
	strand: str
	phase: int


@dataclass
class TranscriptCDS:
	gene_id: str
	transcript_id: str
	chrom: str
	strand: str
	parts: list[CDSPart] = field(default_factory=list)
	tags: set[str] = field(default_factory=set)


def parse_gene_id(attributes: str) -> str | None:
	match = GENE_ID_RE.search(attributes)
	return match.group(1) if match else None


def parse_transcript_id(attributes: str) -> str | None:
	match = TRANSCRIPT_ID_RE.search(attributes)
	return match.group(1) if match else None


def parse_tags(attributes: str) -> set[str]:
	return set(TAG_RE.findall(attributes))


def reverse_complement(seq: str) -> str:
	table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
	return seq.translate(table)[::-1]


def build_fasta_index(fasta_path: Path) -> dict[str, FastaRecordMeta]:
	index: dict[str, FastaRecordMeta] = {}

	with fasta_path.open("rb") as handle:
		current_name: str | None = None
		seq_start: int | None = None
		line_bases: int | None = None
		line_width: int | None = None
		seq_length = 0

		while True:
			line = handle.readline()
			if not line:
				if current_name is not None and seq_start is not None and line_bases and line_width:
					index[current_name] = FastaRecordMeta(
						start_offset=seq_start,
						line_bases=line_bases,
						line_width=line_width,
						length=seq_length,
					)
				break

			if line.startswith(b">"):
				if current_name is not None and seq_start is not None and line_bases and line_width:
					index[current_name] = FastaRecordMeta(
						start_offset=seq_start,
						line_bases=line_bases,
						line_width=line_width,
						length=seq_length,
					)

				current_name = line[1:].decode("utf-8").strip().split()[0]
				seq_start = handle.tell()
				line_bases = None
				line_width = None
				seq_length = 0
				continue

			stripped = line.rstrip(b"\r\n")
			if not stripped:
				continue

			if line_bases is None:
				line_bases = len(stripped)
				line_width = len(line)

			seq_length += len(stripped)

	return index


def fetch_subsequence(
	handle,
	index: dict[str, FastaRecordMeta],
	chrom: str,
	start_1based: int,
	end_1based: int,
) -> str:
	if chrom not in index:
		raise KeyError(f"Chromosome/scaffold '{chrom}' not found in FASTA.")

	meta = index[chrom]
	if start_1based < 1 or end_1based > meta.length or start_1based > end_1based:
		raise ValueError(
			f"Invalid coordinate range {chrom}:{start_1based}-{end_1based}; "
			f"record length is {meta.length}."
		)

	start0 = start_1based - 1
	end0 = end_1based - 1

	start_line = start0 // meta.line_bases
	start_col = start0 % meta.line_bases

	end_line = end0 // meta.line_bases
	end_col = end0 % meta.line_bases

	byte_start = meta.start_offset + start_line * meta.line_width + start_col
	byte_end = meta.start_offset + end_line * meta.line_width + end_col
	nbytes = byte_end - byte_start + 1

	handle.seek(byte_start)
	chunk = handle.read(nbytes)
	return chunk.replace(b"\n", b"").replace(b"\r", b"").decode("ascii")


def load_transcript_cds_parts(cds_gtf_path: Path) -> dict[str, TranscriptCDS]:
	# Some annotations reuse transcript_id values (for example "unknown_transcript_1")
	# across different genes/loci, so we key by a composite identity internally.
	transcripts: dict[tuple[str, str, str, str], TranscriptCDS] = {}

	with cds_gtf_path.open("r", encoding="utf-8") as handle:
		for raw_line in handle:
			if not raw_line or raw_line.startswith("#"):
				continue

			cols = raw_line.rstrip("\n").split("\t")
			if len(cols) < 9:
				continue

			if cols[2] != "CDS":
				continue

			chrom = cols[0]
			start = int(cols[3])
			end = int(cols[4])
			strand = cols[6]
			phase = int(cols[7])
			attributes = cols[8]

			gene_id = parse_gene_id(attributes)
			transcript_id = parse_transcript_id(attributes)
			if gene_id is None or transcript_id is None:
				continue

			transcript_key = (gene_id, transcript_id, chrom, strand)

			if transcript_key not in transcripts:
				transcripts[transcript_key] = TranscriptCDS(
					gene_id=gene_id,
					transcript_id=transcript_id,
					chrom=chrom,
					strand=strand,
				)

			transcript = transcripts[transcript_key]

			transcript.parts.append(CDSPart(chrom=chrom, start=start, end=end, strand=strand, phase=phase))
			transcript.tags.update(parse_tags(attributes))

	# Downstream code only needs the TranscriptCDS records.
	return {f"{i}": t for i, t in enumerate(transcripts.values())}


def ordered_parts(parts: list[CDSPart], strand: str) -> list[CDSPart]:
	reverse = strand == "-"
	return sorted(parts, key=lambda part: part.start, reverse=reverse)


def build_gene_transcript_sequence_dict(
	cds_gtf_path: Path,
	fasta_path: Path,
) -> dict[str, dict[str, dict[str, object]]]:
	transcripts = load_transcript_cds_parts(cds_gtf_path)
	fasta_index = build_fasta_index(fasta_path)

	result: dict[str, dict[str, dict[str, object]]] = defaultdict(dict)

	with fasta_path.open("rb") as fasta_handle:
		for transcript_id, transcript in transcripts.items():
			parts_in_order = ordered_parts(transcript.parts, transcript.strand)
			sequence_chunks: list[str] = []

			for part in parts_in_order:
				chunk = fetch_subsequence(
					fasta_handle,
					fasta_index,
					part.chrom,
					part.start,
					part.end,
				)
				if transcript.strand == "-":
					chunk = reverse_complement(chunk)
				sequence_chunks.append(chunk)

			sequence = "".join(sequence_chunks)
			first_phase = parts_in_order[0].phase
			first_complete_codon_offset = first_phase
			missing_start_bases = 0 if first_phase == 0 else 3 - first_phase

			result[transcript.gene_id][transcript_id] = {
				"sequence": sequence,
				"chromosome": transcript.chrom,
				"strand": transcript.strand,
				"cds_length": len(sequence),
				"tags": sorted(transcript.tags),
				"cds_start_nf": "cds_start_NF" in transcript.tags,
				"mrna_start_nf": "mRNA_start_NF" in transcript.tags,
				"first_cds_phase": first_phase,
				"first_complete_codon_offset": first_complete_codon_offset,
				"missing_start_bases_for_partial_codon": missing_start_bases,
			}

	return dict(result)


def main() -> None:
	parser = argparse.ArgumentParser(
		description=(
			"Build a dictionary of gene_id -> transcript_id -> CDS sequence from filtered "
			"CDS GTF rows using a reference FASTA."
		)
	)
	parser.add_argument(
		"cds_gtf",
		help="Path to the filtered CDS GTF output (from previous step).",
	)
	parser.add_argument(
		"--fasta",
		default="HumanGenome/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
		help="Reference FASTA path.",
	)
	parser.add_argument(
		"-o",
		"--output",
		default="gene_id_to_transcripts.json",
		help="Output JSON file for gene_id -> transcript_id -> metadata dictionary.",
	)

	args = parser.parse_args()
	cds_path = Path(args.cds_gtf)
	fasta_path = Path(args.fasta)
	output_path = Path(args.output)

	gene_to_transcripts = build_gene_transcript_sequence_dict(cds_path, fasta_path)

	with output_path.open("w", encoding="utf-8") as out:
		json.dump(gene_to_transcripts, out, indent=2)

	transcript_count = sum(len(transcripts) for transcripts in gene_to_transcripts.values())
	print(
		f"Wrote {len(gene_to_transcripts)} genes and {transcript_count} transcripts to {output_path}"
	)


if __name__ == "__main__":
	main()
