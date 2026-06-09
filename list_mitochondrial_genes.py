from __future__ import annotations

import argparse
import csv
import re
import urllib.error
import urllib.request
from pathlib import Path


GENE_ID_RE = re.compile(r'gene_id\s+"([^"]+)"')
GENE_NAME_RE = re.compile(r'gene_name\s+"([^"]+)"')
GENE_BIOTYPE_RE = re.compile(r'gene_biotype\s+"([^"]+)"')

NCBI_MITOCHONDRIAL_ACCESSIONS = {
	"Bat_MyotisLucifugus": "NC_029849.1",
	"Bat_RhinolophusFerrumequinum": "MK089295.1",
	"Pangolin_ManisJavanica": "NC_026781.1",
	"Pangolin_ManisPentadactyla": "NC_016008.1",
}

PROTEIN_GENE_MAP = {
	"ATP6": "MT-ATP6",
	"ATP8": "MT-ATP8",
	"COX1": "MT-CO1",
	"COX2": "MT-CO2",
	"COX3": "MT-CO3",
	"CYTB": "MT-CYB",
	"ND1": "MT-ND1",
	"ND2": "MT-ND2",
	"ND3": "MT-ND3",
	"ND4": "MT-ND4",
	"ND4L": "MT-ND4L",
	"ND5": "MT-ND5",
	"ND6": "MT-ND6",
}

TRNA_GENE_MAP = {
	"ALA": "MT-TA",
	"ARG": "MT-TR",
	"ASN": "MT-TN",
	"ASP": "MT-TD",
	"CYS": "MT-TC",
	"GLN": "MT-TQ",
	"GLU": "MT-TE",
	"GLY": "MT-TG",
	"HIS": "MT-TH",
	"ILE": "MT-TI",
	"LYS": "MT-TK",
	"MET": "MT-TM",
	"PHE": "MT-TF",
	"PRO": "MT-TP",
	"THR": "MT-TT",
	"TRP": "MT-TW",
	"TYR": "MT-TY",
	"VAL": "MT-TV",
}

RRNA_GENE_MAP = {
	"12S RIBOSOMAL RNA": "MT-RNR1",
	"16S RIBOSOMAL RNA": "MT-RNR2",
	"S-RRNA": "MT-RNR1",
	"L-RRNA": "MT-RNR2",
}

GENBANK_URL_TEMPLATE = (
	"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
	"?db=nuccore&id={accession}&rettype=gb&retmode=text"
)


def parse_attr(pattern: re.Pattern[str], attributes: str) -> str | None:
	match = pattern.search(attributes)
	return match.group(1) if match else None


def is_mitochondrial(seqname: str, gene_name: str | None, gene_biotype: str | None) -> bool:
	seq = seqname.strip().upper()
	if seq in {"MT", "CHRM", "M"}:
		return True

	name = (gene_name or "").upper()
	if name.startswith("MT-") or name in {"MTRNR1", "MTRNR2"}:
		return True

	biotype = (gene_biotype or "").lower()
	if biotype.startswith("mt_") or biotype == "mitochondrial":
		return True

	return False


def extract_qualifier(block_text: str, qualifier: str) -> str | None:
	match = re.search(rf'/{re.escape(qualifier)}="(.*?)"', block_text, re.S)
	if match:
		return match.group(1)

	match = re.search(rf'/{re.escape(qualifier)}=([^\s/]+)', block_text)
	return match.group(1) if match else None


def extract_gene_id(block_text: str) -> str | None:
	match = re.search(r"GeneID:?\s*\[?(\d+)\]?", block_text)
	return match.group(1) if match else None


def iter_genbank_feature_blocks(record_text: str):
	in_features = False
	current_key: str | None = None
	current_lines: list[str] = []

	for raw_line in record_text.splitlines():
		if raw_line.startswith("FEATURES             Location/Qualifiers"):
			in_features = True
			continue

		if not in_features:
			continue

		if raw_line.startswith("ORIGIN"):
			break

		feature_key = raw_line[5:21].strip() if len(raw_line) >= 21 else ""
		if feature_key:
			if current_key is not None:
				yield current_key, "\n".join(current_lines)
			current_key = feature_key
			current_lines = [raw_line]
		elif current_key is not None:
			current_lines.append(raw_line)

	if current_key is not None:
		yield current_key, "\n".join(current_lines)


def normalize_remote_mitochondrial_gene(feature_key: str, block_text: str) -> tuple[str, str] | None:
	key = feature_key.strip()
	if key == "gene":
		gene_name = extract_qualifier(block_text, "gene")
		if not gene_name:
			return None
		mapped = PROTEIN_GENE_MAP.get(gene_name.upper())
		return (mapped, "protein_coding") if mapped else None

	if key == "tRNA":
		product = extract_qualifier(block_text, "product") or ""
		amino_acid = product.replace("tRNA-", "", 1).strip().upper()
		hint_text = block_text.upper()

		if amino_acid == "LEU":
			if "UUR" in hint_text:
				return "MT-TL1", "tRNA"
			if "CUN" in hint_text:
				return "MT-TL2", "tRNA"
			return None

		if amino_acid == "SER":
			if "UCN" in hint_text:
				return "MT-TS1", "tRNA"
			if "AGY" in hint_text:
				return "MT-TS2", "tRNA"
			return None

		gene_name = TRNA_GENE_MAP.get(amino_acid)
		return (gene_name, "tRNA") if gene_name else None

	if key == "rRNA":
		product = (extract_qualifier(block_text, "product") or "").upper()
		note = (extract_qualifier(block_text, "note") or "").upper()
		hint_text = f"{product} {note}"

		for marker, gene_name in RRNA_GENE_MAP.items():
			if marker in hint_text:
				return gene_name, "rRNA"

		return None

	return None


def fetch_genbank_record(accession: str) -> str:
	url = GENBANK_URL_TEMPLATE.format(accession=accession)
	request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
	with urllib.request.urlopen(request, timeout=30) as response:
		return response.read().decode("utf-8")


def iter_remote_mitochondrial_genes(accession: str):
	record_text = fetch_genbank_record(accession)
	seen: set[str] = set()

	for feature_key, block_text in iter_genbank_feature_blocks(record_text):
		if feature_key not in {"gene", "tRNA", "rRNA"}:
			continue

		normalized = normalize_remote_mitochondrial_gene(feature_key, block_text)
		if normalized is None:
			continue

		gene_name, gene_biotype = normalized
		if not gene_name or gene_name in seen:
			continue

		seen.add(gene_name)
		gene_id = extract_gene_id(block_text)
		yield {
			"gene_id": f"GeneID:{gene_id}" if gene_id else f"{accession}:{gene_name}",
			"gene_name": gene_name,
			"gene_biotype": gene_biotype,
			"chromosome": "MT",
			"gtf": accession,
		}


def iter_mitochondrial_genes(gtf_path: Path):
	seen: set[tuple[str, str]] = set()

	with gtf_path.open("r", encoding="utf-8") as handle:
		for raw_line in handle:
			if not raw_line or raw_line.startswith("#"):
				continue

			cols = raw_line.rstrip("\n").split("\t")
			if len(cols) < 9 or cols[2] != "gene":
				continue

			seqname = cols[0]
			attributes = cols[8]
			gene_id = parse_attr(GENE_ID_RE, attributes)
			gene_name = parse_attr(GENE_NAME_RE, attributes)
			gene_biotype = parse_attr(GENE_BIOTYPE_RE, attributes)

			if gene_id is None:
				continue

			if not is_mitochondrial(seqname, gene_name, gene_biotype):
				continue

			key = (gene_id, gene_name or "")
			if key in seen:
				continue

			seen.add(key)
			yield {
				"gene_id": gene_id,
				"gene_name": gene_name or "",
				"gene_biotype": gene_biotype or "",
				"chromosome": seqname,
				"gtf": gtf_path,
			}


def find_gtf_files(genomes_root: Path):
	return sorted(p for p in genomes_root.rglob("*.gtf") if p.is_file())


def species_name_from_gtf(gtf_path: Path, genomes_root: Path) -> str:
	relative_parts = gtf_path.relative_to(genomes_root).parts
	if len(relative_parts) >= 2:
		return relative_parts[0]
	return gtf_path.stem


def iter_records_for_species(species_dir: Path, genomes_root: Path, source: str):
	species = species_dir.name
	accession = NCBI_MITOCHONDRIAL_ACCESSIONS.get(species)

	if source != "local" and accession:
		try:
			for record in iter_remote_mitochondrial_genes(accession):
				yield {
					**record,
					"species": species,
				}
			return
		except (urllib.error.URLError, TimeoutError, ValueError):
			if source == "ncbi":
				raise

	for gtf_path in sorted(species_dir.rglob("*.gtf")):
		if not gtf_path.is_file():
			continue
		for record in iter_mitochondrial_genes(gtf_path):
			yield {
				**record,
				"species": species_name_from_gtf(gtf_path, genomes_root),
			}


def main() -> int:
	parser = argparse.ArgumentParser(
		description="List mitochondrial genes from genome GTF files under Genomes/.",
	)
	parser.add_argument(
		"--genomes-root",
		default="Genomes",
		help="Root directory containing species genome annotations.",
	)
	parser.add_argument(
		"--species",
		help="Optional species folder name to restrict the scan.",
	)
	parser.add_argument(
		"--source",
		choices=("auto", "local", "ncbi"),
		default="auto",
		help="Use local GTFs, NCBI mitochondrial references, or auto-select NCBI for supported species.",
	)
	args = parser.parse_args()

	genomes_root = Path(args.genomes_root)

	writer = csv.writer(__import__("sys").stdout, delimiter="\t", lineterminator="\n")
	writer.writerow(["species", "gene_id", "gene_name", "gene_biotype", "chromosome", "gtf_file"])

	if args.species:
		species_dirs = [genomes_root / args.species]
	else:
		species_dirs = sorted(p for p in genomes_root.iterdir() if p.is_dir())

	for species_dir in species_dirs:
		if not species_dir.exists():
			continue
		for record in iter_records_for_species(species_dir, genomes_root, args.source):
			writer.writerow([
				record["species"],
				record["gene_id"],
				record["gene_name"],
				record["gene_biotype"],
				record["chromosome"],
				str(record["gtf"]),
			])

	return 0


if __name__ == "__main__":
	raise SystemExit(main())