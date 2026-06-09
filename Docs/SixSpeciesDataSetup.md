# Six-Species Expansion Setup

This document records genome and expression sources for six additional species and maps them to repository paths already used by the existing codon pipeline.

## Added Species and Paths

- Ecoli_K12_MG1655
  - Genome path: Genomes/Ecoli_K12_MG1655
  - Expression path: Expression/Ecoli_K12_MG1655
  - Results path: Results2/Ecoli_K12_MG1655
- Yeast_S288C
  - Genome path: Genomes/Yeast_S288C
  - Expression path: Expression/Yeast_S288C
  - Results path: Results2/Yeast_S288C
- RaccoonDog_NyctereutesProcyonoides
  - Genome path: Genomes/RaccoonDog_NyctereutesProcyonoides
  - Expression path: Expression/RaccoonDog_NyctereutesProcyonoides
  - Results path: Results2/RaccoonDog_NyctereutesProcyonoides
- Mink_NeogaleVison
  - Genome path: Genomes/Mink_NeogaleVison
  - Expression path: Expression/Mink_NeogaleVison
  - Results path: Results2/Mink_NeogaleVison
- Bat_RhinolophusAffinis
  - Genome path: Genomes/Bat_RhinolophusAffinis
  - Expression path: Expression/Bat_RhinolophusAffinis
  - Results path: Results2/Bat_RhinolophusAffinis
- Bat_RhinolophusSinicus
  - Genome path: Genomes/Bat_RhinolophusSinicus
  - Expression path: Expression/Bat_RhinolophusSinicus
  - Results path: Results2/Bat_RhinolophusSinicus

## Source Summary

Full source details are in Docs/SixSpeciesManifest.csv. Genome assemblies selected:

- Escherichia coli K-12 MG1655: GCF_000005845.2
- Saccharomyces cerevisiae S288C: GCF_000146045.2
- Nyctereutes procyonoides: GCF_905146905.1
- Neogale vison: GCF_020171115.1
- Rhinolophus affinis: GCA_043728065.1
- Rhinolophus sinicus: GCF_036562045.2

Representative RNA-seq runs selected for initial processing:

- E. coli: DRR089608
- S. cerevisiae: DRR212838
- Nyctereutes procyonoides: SRR24920857
- Neogale vison: ERR12791496
- Rhinolophus affinis: SRR12145327
- Rhinolophus sinicus: DRR734818

## Required Files Per Species

To work with the existing programs, each species should populate the following files:

- Genomes/<SpeciesKey>/genomic.fna
- Genomes/<SpeciesKey>/genomic.gtf
- Genomes/<SpeciesKey>/t  (CDS-only GTF rows)
- Genomes/<SpeciesKey>/gene_id_to_transcripts.json
- Expression/<SpeciesKey>/ExpressionData.csv

The table format expected for expression transfer into Results2 is:

- Header column 1: Gene_ID
- Remaining columns: sample names (for example DRR089608)

## Integration With Existing Programs

1. Build CDS transcript JSON from the CDS-filtered annotation and FASTA.

- Script: CDS_to_sequence.py
- Output: Genomes/<SpeciesKey>/gene_id_to_transcripts.json

2. Build expression table from RNA-seq alignment/counting pipeline.

- Align reads to Genomes/<SpeciesKey>/genomic.fna
- Count genes with featureCounts using Genomes/<SpeciesKey>/genomic.gtf
- Filter to gene IDs present in gene_id_to_transcripts.json
- Write Expression/<SpeciesKey>/ExpressionData.csv

Reusable command wrapper in this repository:

bash Expression/run_species_expression_from_ena.sh <SpeciesKey> <RunAccession> <Threads>

Examples:

- bash Expression/run_species_expression_from_ena.sh Ecoli_K12_MG1655 DRR089608 4
- bash Expression/run_species_expression_from_ena.sh Yeast_S288C DRR212838 4
- bash Expression/run_species_expression_from_ena.sh RaccoonDog_NyctereutesProcyonoides SRR24920857 4
- bash Expression/run_species_expression_from_ena.sh Mink_NeogaleVison ERR12791496 4
- bash Expression/run_species_expression_from_ena.sh Bat_RhinolophusAffinis SRR12145327 4
- bash Expression/run_species_expression_from_ena.sh Bat_RhinolophusSinicus DRR734818 4

3. Copy expression table and run codon frequency export into Results2.

Per-species command template:

python codon_frequency.py --json_file Genomes/<SpeciesKey>/gene_id_to_transcripts.json -o Results2/<SpeciesKey>/stuff.npz

and copy expression table:

cp Expression/<SpeciesKey>/ExpressionData.csv Results2/<SpeciesKey>/ExpressionData.csv

Repository helper script (all six species, or one with --species):

bash run_results2_codon_frequency.sh

Example single-species run:

bash run_results2_codon_frequency.sh --species Ecoli_K12_MG1655

## Notes

- E. coli and yeast are valid negative controls but represent distinct biology compared with mammalian hosts.
- Rhinolophus affinis assembly is currently listed as GenBank (GCA) in the selected source.
- Use paired-end counting settings only when run layout is paired.

## Current Verification Status

An end-to-end smoke test has been completed for E. coli using the revised program flow (CDS extraction compatibility + Results2 codon_frequency export):

- Generated files in Genomes/Ecoli_K12_MG1655:
  - genomic.fna
  - genomic.gff
  - genomic.gtf
  - t
  - gene_id_to_transcripts.json
- Generated files in Expression/Ecoli_K12_MG1655:
  - ExpressionData.csv (smoke-test table, sample name `SMOKE_ECOLI`)
- Generated Results2 smoke-test artifacts:
  - Results2/Ecoli_K12_MG1655/ExpressionData.csv
  - Results2/Ecoli_K12_MG1655/stuff.npz

This confirms that the newly added species folder/path conventions are compatible with the existing codon analysis scripts.
