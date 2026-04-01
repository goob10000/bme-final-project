# Final Project Notes

## Notes from presentation
1. Ensure that I have at least 15-20 species worth of data for my PCA and clustering.
1. Ensure that when I'm doing PCA that I use a method specifically for data that adds up to 1.
1. Investigate what the genomic identity % is typically for conclusive results about where viruses and corona viruses came from.

## More Notes

Ensembl_canonical - Revealed when you want the canonical gene. Pick these
CDS - Coding sequence
ZNF695 - Good gene
GTF is 1 offset, not 0 offset

Figure out what type of cells COVID likes to proliferate in. Use these as reference for codon preference in highly expressed genes. Compare to an unrelated tissue as a control.

Compare a person with covid vs a person without covid

If a large portion of genes are missing then have to change something about the process

Do quality check to make sure that the genes that are highly expressed make sense

Quick test:
- Compare codon preferences of different cell types in the same species to see if it varies greatly

If I can't find expression data for a specific species, blast the sequences against human sequences.
If gene names are different, if I have the transcripts I can blast against the pangolin dataset.


gff/gtf files for genes and their locations.
Some kind of gene expression matrix to accompany it. Or just fish out the genes that are similar to highly expressed genes in humans (based on assumption of similar genetics to other mammals/vertebrates)

- [Human #1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318035) 
- [Human #2](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275264)
- [Pig #1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231494)
- [Mouse #1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287197)
- [Myotis Lucifugus](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262772) 
- [Rhinolophus Ferrumequinum](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274964)


Covid:
run_accession	sample_accession	scientific_name	library_strategy	fastq_ftp
SRR11412217	SAMN14444845	Homo sapiens	RNA-Seq	ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/017/SRR11412217/SRR11412217.fastq.gz
SRR11412218	SAMN14444845	Homo sapiens	RNA-Seq	ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/018/SRR11412218/SRR11412218.fastq.gz
SRR11412222	SAMN14444844	Homo sapiens	RNA-Seq	ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/022/SRR11412222/SRR11412222.fastq.gz
SRR11412223	SAMN14444843	Homo sapiens	RNA-Seq	ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/023/SRR11412223/SRR11412223.fastq.gz
SRR11412225	SAMN14444843	Homo sapiens	RNA-Seq	ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/025/SRR11412225/SRR11412225.fastq.gz
SRR11412226	SAMN14444843	Homo sapiens	RNA-Seq	ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/026/SRR11412226/SRR11412226.fastq.gz
SRR11412229	SAMN14444842	Homo sapiens	RNA-Seq	ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/029/SRR11412229/SRR11412229.fastq.gz
SRR11412232	SAMN14444841	Homo sapiens	RNA-Seq	ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/032/SRR11412232/SRR11412232.fastq.gz
SRR11412235	SAMN14444840	Homo sapiens	RNA-Seq	ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/035/SRR11412235/SRR11412235.fastq.gz
SRR11412239	SAMN14444839	Homo sapiens	RNA-Seq	ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/039/SRR11412239/SRR11412239.fastq.gz
SRR11412248	SAMN14444837	Homo sapiens	RNA-Seq	ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/048/SRR11412248/SRR11412248.fastq.gz
SRR11412250	SAMN14444837	Homo sapiens	RNA-Seq	ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/050/SRR11412250/SRR11412250.fastq.gz


## Ideas

1. Need a way to validate whether the genomes are actually the genome I think they are
    - https://github.com/jenniferlu717/KrakenTools
1. Check that top genes actually represent genes that make sense (like cellular respiration, etc.)
1. Compare distances of the two bats and two pangolins showing that their codon preferences are more similar than that of other species
1. racoon dogs as another relevant species

## Stretch Goals

1. Show how codon preference can be used to generate a phylogenetic tree for a group of species?



Look at immune response genes for COVID
- Codon bias influenza response -> Influenza causes certain genes to be much more highly expressed and likely alters the codon biases of the cell


## High Frequency Gene Validation

a simple check that coding sequence is multiple of 3 long starts with ATG and ends with stop codon will catch most the junk. you can check for ORF (no stop codon internal) as one more check

Also adjust for TPM

### Top Genes: Covid
1. N
1. ORF10
1. M
1. S
1. ORF1ab
### Top Genes: Human1
1. ENSG00000236637
    - Interferon Alpha 4
    - [ncbi](https://www.ncbi.nlm.nih.gov/gene?Cmd=DetailsSearch&Term=3441)
    - Seems likely that this is a response to the infection
1. ENSG00000074047
    - GLI family zinc finger 2
    - [ncbi](https://www.ncbi.nlm.nih.gov/gene/?term=ENSG00000074047)
    - Seems unlikely this should be a highly expressed gene (although it may be tissue dependent)
    - GPT:
        - Rationale: GLI2 is a Hedgehog-pathway transcription factor usually low/moderate and inducible; infections that alter developmental/signaling pathways or cause tissue repair/remodeling can change GLI2 expression.
        - Evidence: Some transcriptomic studies of SARS-CoV-2 show altered developmental/signaling pathways and transcription-factor changes, but I found no consistent, widely reported direct upregulation of GLI2 across COVID-19 datasets.
        - Conclusion: GLI2 could be elevated in specific cell types or at particular stages of infection (e.g., injured airway epithelium, fibrotic/remodeling responses), but you should check RNA-seq/GEO/GTEx/Single-cell COVID datasets for the tissue and timepoint of interest to confirm.
1. ENSG00000234829
    - interferon alpha 17
    - [ncbi](https://www.ncbi.nlm.nih.gov/gene/?term=ENSG00000234829)
1. ENSG00000160856
    - Fc receptor like 3
    - [ncbi](https://www.ncbi.nlm.nih.gov/gene/?term=ENSG00000160856)
1. ENSG00000179674
    - ARF like GTPase 14
    - [ncbi](https://www.ncbi.nlm.nih.gov/gene/?term=ENSG00000179674)
### Top Genes: Human2
1. ENSG00000149021
    - secretoglobin family 1A member 1
    - Lung tissue specific
    - [ncbi](https://www.ncbi.nlm.nih.gov/gene/?term=ENSG00000149021)
1. ENSG00000161055
    - secretoglobin family 3A member 1
    - [ncbi](https://www.ncbi.nlm.nih.gov/gene/?term=ENSG00000161055)
1. ENSG00000124107
    - secretory leukocyte peptidase inhibitor
    - [ncbi](https://www.ncbi.nlm.nih.gov/gene/?term=ENSG00000124107)
1. ENSG00000198712
    - mitochondrially encoded cytochrome c oxidase II
    - [ncbi](https://www.ncbi.nlm.nih.gov/gene/?term=ENSG00000198712)
1. ENSG00000198938
    - mitochondrially encoded cytochrome c oxidase III
    - [ncbi](https://www.ncbi.nlm.nih.gov/gene/?term=ENSG00000198938)
### Top Genes: Mouse
1. ENSMUSG00000024653
    - secretoglobin, family 1A, member 1
    - [ncbi](https://www.ncbi.nlm.nih.gov/gene/?term=ENSMUSG00000024653)
1. ENSMUSG00000052305
    - hemoglobin, beta adult s chain
    - [ncbi](https://www.ncbi.nlm.nih.gov/gene/?term=ENSMUSG00000052305)
1. ENSMUSG00000069919
    - hemoglobin alpha, adult chain 1
    - [ncbi](https://www.netflix.com/watch/81044688?trackId=155573558)
1. ENSMUSG00000064357
    - ATP synthase 6, mitochondrial
    - [ncbi](https://www.ncbi.nlm.nih.gov/gene/?term=ENSMUSG00000064357)
1. ENSMUSG00000064358
    - cytochrome c oxidase III, mitochondrial
    - [ncbi](https://www.ncbi.nlm.nih.gov/gene/?term=ENSMUSG00000064358)
### Top Genes: Pig
1. ENSSSCG00000014540
    - ferritin heavy chain 1
    - [ncbi](https://www.ncbi.nlm.nih.gov/gene/?term=ENSSSCG00000014540)
1. ENSSSCG00000008245
    - thymosin beta 10
    - [ncbi](https://www.ncbi.nlm.nih.gov/gene/?term=ENSSSCG00000008245)
1. ENSSSCG00000012119
    - thymosin beta 4 X-linked
    - [ncbi](https://www.ncbi.nlm.nih.gov/gene/?term=ENSSSCG00000012119)
1. ENSSSCG00000004687
    - beta-2-microglobulin
    - [ncbi](https://www.ncbi.nlm.nih.gov/gene/?term=ENSSSCG00000004687)
1. ENSSSCG00000032857
    - S100 calcium binding protein A12
    - [ncbi](https://www.ncbi.nlm.nih.gov/gene/?term=ENSSSCG00000032857)
### Top Genes: Bat_RhinolophusFerrumequinum
1. FN1
1. THBS1
1. COL1A1
1. ACTB
1. COL1A2
### Top Genes: Bat_MyotisLucifugus
1. FN1
1. EEF1A1
1. PABPC1
1. TUBA1B
1. EEF2