# Outline

## Title

Codon bias in highly expressed genes contains information on viral adaptation and lifestyle

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

https://www.scientificamerican.com/article/new-evidence-supports-animal-origin-of-covid-virus-through-raccoon-dogs/

## Methods

1. Codon bias by filtering down gene expression data to one CDS per gene.
1. TPM for gene expression data
1. Filter out ribosomal genes (although they may represent an additional fingerprint as ribosomes must differ slightly between organisms << low priority)
1. Try index weighted by expression rather than a cutoff.
- This would have the effect of basically "codon load"
- and if it remains as proportions (weighted average rather than count) can be compared in similarity to the organism

1. Validation involving looking up gene names to make sure they represent accurate genes for the organisms.
- Viruses typically ship with their own polymerase so effects on tRNA concentration and translation generally seem more important

## Result

1. Graphs demonstrating the similarity of codon preferences across the tissues of a single organism. 
1. Graphs demonstrating codon bias similarity or differences between infected organisms with high interferon expresion 
1. Plots showing the codon bias similarity between covid and humans over the past few years
    1. Worth focusing on silent mutations and making a transition matrix of shfiting codon to codon to see if it moves in a given direction (towards humans) or if it's random
    1. Need some kind of null hypothesis (sortting genes randomly or compare to lowest expression genes instead of highest)
1. Plots showing codon bias similatiry between related organisms. Could be used to generate a phylogenetic tree as an example.
1. Plots demonstrating KL divergence between organisms for various different levels of highest expressed genes
1. Graph of the similarity to a prvious species and humanity over time
1. Worth looking at the same analysis but on all mRNA not just protein coding regions. 

## Conclusion

Covid's path can be traced based on its similarity to the codon preference of other organisms. 

## References:

Hadfield et al, Nextstrain: real-time tracking of pathogen evolution, Bioinformatics (2018)
Sagulenko et al, TreeTime: Maximum-likelihood phylodynamic analysis, Virus Evolution (2017)

1919 flu
- need to separate mutation bias from pressure to adapt to codon bias
- segmented genome. Look in aggregate but consider looking at them individually

Look at HIV too

More organisms to check:
- drosophila (bit of an aside)

1. Show how this feature contains enough information that it can be used to generate a phylogenetic tree consistently and accurately
1. Track the path of SARS-CoV-2 as it moves between species based on the hypothesis that viral genomes must adapt to the codon bias of their host
1. Prove that this analysis can be generalized beyond humans and SARS-CoV-2 to other viruses and species

1. Compare codon bias of different organisms

1. Ecoli and lambda phage

Future Ideas:
- Can you use a similar analysis with codon bias or other mutational pressures to estimate the likelihood of viral jumps for viruses that haven't entered a human population before?
    - Is there some metric of similarity that is required for jumping or increases likelihood for jumping?

I think as long as it’s a high-quality sequence, you can simplify your world to that example and you are certain all of your code is perfectly correct? ;) one sanity check I like to do is visualize the genome you can print it out as wrapped lines of DNA with amino acids under them and then under that the number of the codon preference as a rank so you can see how many 1, 2 verses 3 4 t are being picked as a sanity check if something is further away, you should see a lot more higher numbered pics and if it’s close, you should see a lot of first or second choices by frequency