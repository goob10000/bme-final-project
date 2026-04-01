# Codon Usage Pipeline

End-to-end workflow for extracting canonical coding sequences and analyzing codon proportions by expression level.



## Overview

Three sequential steps:
1. Filter GTF to CDS rows with `Ensembl_canonical` tag
2. Assemble CDS sequences from the reference genome  
3. Compute codon proportions weighted by gene expression

## Optional: Reconstruct Expression Matrix from Human2 (10x)

Use this when GEO provides single-cell 10x files (`.features.tsv/.gz`, `.barcodes.tsv/.gz`, `.matrix.mtx/.gz`) instead of one combined counts CSV.

**File:** `reconstruct_human2_raw_counts.py`

```powershell
python reconstruct_human2_raw_counts.py [--human2-dir DIR] [--needs-series-matrix] [--series-matrix FILE] [--group healthy|ipf|all] [--norm raw|cpm|tpm] [--gene-lengths FILE] [--output OUT.csv]
```

**Arguments:**
- `--human2-dir` – Directory with per-sample 10x files (default: `Expression/Human2`)
- `--needs-series-matrix` – Enable GEO metadata filtering (Healthy/IPF) via series matrix. If omitted, all samples in `--human2-dir` are processed.
- `--series-matrix` – GEO series matrix TXT used to map samples to Healthy/IPF (default: `Expression/Human2/GSE275264_series_matrix.txt`)
- `--group` – Which samples to include (`healthy`, `ipf`, or `all`; default: `healthy`). Ignored unless `--needs-series-matrix` is set.
- `--norm` – Output scale per sample (`raw`, `cpm`, or `tpm`; default: `raw`)
- `--gene-lengths` – Required for TPM; tab-delimited `gene_id<TAB>length_bp`
- `--output` – Output CSV path (default: `Expression/Human2/GSE275264_<norm>_<group>.csv`; when `--needs-series-matrix` is not set, default group becomes `all`)

**Output:** CSV with `Gene_ID` rows and one donor/sample column per selected GEO sample.

**Examples:**

All detected samples without series metadata:
```powershell
python reconstruct_human2_raw_counts.py --norm raw
```

Healthy-only raw counts:
```powershell
python reconstruct_human2_raw_counts.py --needs-series-matrix --group healthy --norm raw
```

Healthy-only CPM:
```powershell
python reconstruct_human2_raw_counts.py --needs-series-matrix --group healthy --norm cpm
```

Healthy-only TPM (requires gene lengths):
```powershell
python reconstruct_human2_raw_counts.py --needs-series-matrix --group healthy --norm tpm --gene-lengths HumanGenome/gene_lengths.tsv
```

All donors with explicit output path:
```powershell
python reconstruct_human2_raw_counts.py --needs-series-matrix --group all --norm raw --output Expression/Human2/GSE275264_raw_counts_all.csv
```

## Optional: Extract CDS Sequences from GBFF

Use this when your reference annotation is GenBank format (`.gbff` or `.gbff.gz`) instead of GTF+FASTA.

**File:** `gbff_to_sequence.py`

```powershell
python gbff_to_sequence.py <input.gbff|input.gbff.gz> [-o output.json]
```

**Arguments:**
- `input.gbff|input.gbff.gz` - GenBank file containing CDS feature annotations
- `-o, --output` - Output JSON path (default: `gene_id_to_sequence_transcripts.json`)

**Output:** JSON with structure `gene_id -> transcript_id -> metadata`, including `sequence`, `strand`, `cds_length`, and codon-phase fields compatible with downstream codon-usage analysis.

**Examples:**

From uncompressed GBFF:
```powershell
python gbff_to_sequence.py Genomes/Covid/reference.gbff -o Genomes/Covid/gene_id_to_sequence_transcripts.json
```

From gzipped GBFF:
```powershell
python gbff_to_sequence.py Genomes/Covid/reference.gbff.gz -o Genomes/Covid/gene_id_to_sequence_transcripts.json
```

## Step 1: Filter GTF for Canonical CDS


Open File
```bash
tar -xvzf filename.tar.gz
tar -xvf filename.tar
```

```bash
grep "        CDS.* tag \"Ensembl_canonical\"" *.gtf
```

## Step 2: Assemble CDS Sequences

**File:** `CDS_to_sequence.py`

```powershell
python311 CDS_to_sequence.py <cds_gtf> [--fasta ref_genome.fa] [-o output.json]
```

**Arguments:**
- `cds_gtf` – Filtered GTF from Step 1
- `--fasta` – Reference genome FASTA (default: `HumanGenome/Homo_sapiens.GRCh38.dna.primary_assembly.fa`)
- `-o, --output` – Output JSON path (default: `HumanGenome/gene_id_to_transcripts.json`)

**Output:** JSON with structure `{gene_id: {transcript_id: {sequence, cds_length, tags, cds_start_nf, mrna_start_nf, first_cds_phase, first_complete_codon_offset, ...}}}`

**Example:**
```powershell
python311 CDS_to_sequence.py canonical_cds.gtf -o HumanGenome/gene_id_to_transcripts.json
```

## Step 3: Analyze Codon Proportions by Expression

**File:** `topn_codon_proportions.py`

```powershell
python311 topn_codon_proportions.py [--sample SAMPLE] [--top-n N] [--top-n-series N1,N2,...] [-o output.csv] [--plot-file file.png] [--plot-by-aa-dir dir/]
```

**Arguments:**
- `--transcript-json` – CDS JSON from Step 2 (default: `HumanGenome/gene_id_to_sequence_transcripts.json`)
- `--counts-csv` – Expression counts CSV (default: `HumanGenome/GSE318035_raw_counts.csv`)
- `--sample` – Sample column name (default: first sample in CSV)
- `--top-n` – Primary gene count cutoff (default: 500)
- `--top-n-series` – Comma-separated extra cutoffs for curve generation (default: `100,200,500,1000,2000,3000,4000,5000`)
- `-o, --output` – Output table (.csv, .json, or .parquet)
- `--plot-file` – Save single matplotlib figure to this path (.png, .pdf, .svg)
- `--plot-by-aa-dir` – Save one plot per amino acid to this directory
- `--aa` – Filter to one amino acid (e.g., `L` for Leucine) in plots
- `--plot-dpi` – DPI for image output (default: 160)

**Output:** Long-format table with columns: `sample, top_n, aa, codon, weighted_codon_count, aa_total_weighted, codon_proportion`

**Integration note (GBFF + reconstructed counts):** If you generated transcript sequences with `gbff_to_sequence.py` and expression with `reconstruct_human2_raw_counts.py`, pass both explicitly to this step.

```powershell
python311 topn_codon_proportions.py --transcript-json Genomes/Covid/gene_id_to_sequence_transcripts.json --counts-csv Expression/Human2/GSE275264_raw_healthy.csv --sample H1 --top-n 500 -o Genomes/Covid/topn_codon_proportions_human2_healthy.csv
```

**Examples:**

Generate table and single plot:
```bash
python311 topn_codon_proportions.py --sample COVID_1 --top-n 500 --top-n-series 1,2,3,4,5,10,20,40,60,80,100,200,500,1000,2000,3000,4000,5000 -o HumanGenome/topn_codon_proportions_covid1.csv --plot-file Expression/Human1/topn_codon_proportions_covid1.png
```

Generate table and per-amino-acid plots:
```bash
python311 topn_codon_proportions.py --sample COVID_1 --top-n 500 --top-n-series 1,2,3,4,5,10,20,40,60,80,100,200,500,1000,2000,3000,4000,5000 -o HumanGenome/topn_codon_proportions_covid1.csv --plot-file Expression/Human1/topn_codon_proportions_covid1.png --plot-by-aa-dir Expression/Human1/aa_plots_small
```

Plot only Leucine (multi-codon):
```powershell
python311 topn_codon_proportions.py --sample COVID_1 --aa L --plot-file HumanGenome/codon_prop_L.png
```

## Full Pipeline Example

```powershell
# Step 1: Filter GTF
python311 temp.py t -o canonical_cds.gtf

# Step 2: Assemble sequences  
python311 CDS_to_sequence.py canonical_cds.gtf -o HumanGenome/gene_id_to_transcripts.json

# Step 3: Analyze codon usage with expression weighting
python311 topn_codon_proportions.py --sample COVID_1 --top-n 500 --top-n-series 100,200,500,1000,2000,3000,4000,5000 -o HumanGenome/topn_codon_proportions.csv --plot-by-aa-dir HumanGenome/aa_plots --plot-dpi 220
```

## Key Concepts

**first_complete_codon_offset** – Number of bases to skip at the start of a CDS before reading complete codons. Used for transcripts with `cds_start_NF` or `mRNA_start_NF` tags (incomplete 5' ends).

**codon_proportion** – For each amino acid, the weighted count of a specific codon divided by the total weighted count for that amino acid. Weights come from per-gene expression values.

**Singleton-AA codons excluded** – Codons `ATG` (M) and `TGG` (W) are dropped from codon frequency arrays since they have no synonymous alternatives.
