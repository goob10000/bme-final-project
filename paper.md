# Outline

## Introduction

In this paper I demonstrate how codon usage bias in highly expressed genes can be used to illucidate the path viruses take between hosts and how they adapt to the codon preferences of said hosts

## Hypothesis

1. Virus should adapt to human codon bias over time
1. Codon bias should show traces of its most recent species via its codon usage
1. Codon bias should be most striking in highly expressed genes where the codon bias is most evident

## Background

Cite other papers
- Other virus transfers and showing adaptation to humans
1. Negative control showing how codon bias between a mouse virus and a human are different.
1. Positive control showing how codon bias between of HIV and humans are very similar.
1. Positive control showing how codon bias between mouse virus and mouse are very similar
- Highest expressed genes exhibit different codon biases to the rest of the genes
1. Graphs demonstrating codon bias and how the highest expressed genes in an organism exaggerate that relationship

## Methods

1. Codon bias by filtering down gene expression data to one CDS per gene.
1. TPM for gene expression data
1. Filter out ribosomal genes (although they may represent an additional fingerprint as ribosomes must differ slightly between organisms << low priority)
1. Try index weighted by expression rather than a cutoff.

1. Validation involving looking up gene names to make sure they represent accurate genes for the organisms.

## Result

1. Graphs demonstrating the similarity of codon preferences across the tissues of a single organism. 
1. Graphs demonstrating codon bias similarity or differences between infected organisms with high interferon expresion 
1. Plots showing the codon bias similarity between covid and humans over the past few years
    1. Worth focusing on silent mutations and making a transition matrix of shfiting codon to codon to see if it moves in a given direction (towards humans) or if it's random
1. Plots showing codon bias similatiry between related organisms. Could be used to generate a phylogenetic tree as an example.
1. Plots demonstrating KL divergence between organisms for various different levels of highest expressed genes
1. Graph of the similarity to a prvious species and humanity over time
1. Worth looking at the same analysis but on all mRNA not just protein coding regions. 

## Conclusion

Covid's path can be traced based on its similarity to the codon preference of other organisms. 