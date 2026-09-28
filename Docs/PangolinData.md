# Pangolin Genome and Expression Data (Methods-Style Summary)

## Data Sources

Reference genomes and annotations were obtained from NCBI assemblies for two pangolin species. For *Manis javanica*, genome sequence and annotation were taken from assembly GCF_014570555.1 (YNU_ManPten_2.0), using the genomic FASTA and GTF files in [Genomes/Pangolin_ManisJavanica](../Genomes/Pangolin_ManisJavanica). For *Manis pentadactyla*, genome sequence and annotation were taken from assembly GCF_030020395.1 (mManPen7.hap1), using the genomic FASTA and GTF files in [Genomes/Pangolin_ManisPentadactyla](../Genomes/Pangolin_ManisPentadactyla).

RNA-seq run candidates were compiled from NCBI SRA runinfo manifests and summarized in [Expression/Pangolin/PANGOLIN_DATA_SOURCING.md](../Expression/Pangolin/PANGOLIN_DATA_SOURCING.md). Species-specific expression tables were generated from representative RNA-seq runs: SRR8943548 for *M. javanica* and SRR29468387 for *M. pentadactyla*.

## Expression Processing

Expression processing was performed with the pangolin pipeline scripts in [Expression/Pangolin](../Expression/Pangolin), including [Expression/Pangolin/build_pangolin_expression_data.sh](../Expression/Pangolin/build_pangolin_expression_data.sh), [Expression/Pangolin/build_pangolin_expression_data_quick.sh](../Expression/Pangolin/build_pangolin_expression_data_quick.sh), and [Expression/Pangolin/run_full_pentadactyla.sh](../Expression/Pangolin/run_full_pentadactyla.sh).

For each species, paired-end FASTQ files were retrieved via ENA/SRA links, aligned to the corresponding reference genome with minimap2 (`-ax sr`), and sorted with samtools. Gene-level read counting was then performed with featureCounts in paired-end mode (`-p --countReadPairs`) using exon features (`-t exon`) and `gene_id` grouping (`-g gene_id`).

Count tables were converted to species-level expression matrices by filtering to valid gene identifiers present in the species CDS-derived gene set, sorting by gene ID, and writing CSV outputs with columns `Gene_ID` and sample run ID. Final outputs are stored at [Expression/Pangolin_ManisJavanica/ExpressionData.csv](../Expression/Pangolin_ManisJavanica/ExpressionData.csv) and [Expression/Pangolin_ManisPentadactyla/ExpressionData.csv](../Expression/Pangolin_ManisPentadactyla/ExpressionData.csv).

## Genome-to-CDS Data Used Downstream

To keep expression and codon analyses on a consistent gene universe, species CDS-based gene dictionaries generated in the genome preprocessing workflow (including the CDS-only annotation file `t` and gene-level JSON outputs in each species genome folder) were used as the gene-ID filter set during expression table construction.

## Reproducibility Commands

The expression tables reported above can be reproduced with the commands below from the repository root.

### Manis javanica

Run the multi-species builder (this generates the *M. javanica* table using SRR8943548):

```bash
bash Expression/Pangolin/build_pangolin_expression_data.sh
```

Primary output artifact:

- [Expression/Pangolin_ManisJavanica/ExpressionData.csv](../Expression/Pangolin_ManisJavanica/ExpressionData.csv)

Key intermediate artifacts:

- [Expression/Pangolin/work/ManisJavanica/SRR8943548.sorted.bam](../Expression/Pangolin/work/ManisJavanica/SRR8943548.sorted.bam)
- [Expression/Pangolin/work/ManisJavanica/SRR8943548.featureCounts.txt](../Expression/Pangolin/work/ManisJavanica/SRR8943548.featureCounts.txt)

### Manis pentadactyla

Run the dedicated full pipeline script (SRR29468387):

```bash
bash Expression/Pangolin/run_full_pentadactyla.sh
```

Primary output artifact:

- [Expression/Pangolin_ManisPentadactyla/ExpressionData.csv](../Expression/Pangolin_ManisPentadactyla/ExpressionData.csv)

Key intermediate artifacts:

- [Expression/Pangolin/work/ManisPentadactyla/SRR29468387.sorted.bam](../Expression/Pangolin/work/ManisPentadactyla/SRR29468387.sorted.bam)
- [Expression/Pangolin/work/ManisPentadactyla/SRR29468387.featureCounts.txt](../Expression/Pangolin/work/ManisPentadactyla/SRR29468387.featureCounts.txt)
- [Expression/Pangolin/work/ManisPentadactyla/SRR29468387.featureCounts.txt.summary](../Expression/Pangolin/work/ManisPentadactyla/SRR29468387.featureCounts.txt.summary)
