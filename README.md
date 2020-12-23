# Analysis of population structure

This repository contains scripts and code to carry out the analyses in:

[Zaidi AA, Mathieson I. Demographic history mediates the effect of stratification on polygenic scores. Elife 2020;9:](https://elifesciences.org/articles/61548)

Please cite this paper if you use any of the code in your own analyses.

This site describes most of the post-simulation analyses carried out in the paper and how each figure was generated.

1. Code for the plot of PCA on genotype data simulated under the recent and perpetual structure model is detailed [here](plt_PCA.html):

2. Simulation of phenotypes (both non-heritable and heritable) are described [here](Simulating_heritable_phenotypes.html).

3. Code for QQ-plot of GWAS on non-heritable phenotypes under different demographic histories is detailed [here]()

4. Code for QQ-plot of gene burden tests is detailed [here](plt_burden_association.html)

5. Code for Gini curves showing the geographically clustered common and rare variants and gene burden can be found [here](plt_burden_clustering.html)

6. Code for studying the population structure in polygenic scores simulated under the 'grid' model of population structure is [here](plottingprs_distribution_gridt.html)

7. Code for studying the population structure in the more complex model of population structure can be found [here](plt_ukb_unrelated_prs.html)

8. Code for studying the population structure in polygenic scores generated from sibling GWAS can be found [here](biasvaccuracy_prsascertainment.html)

9. Code for studying the effect of fine-mapping vs clumping+threshold on stratification and prediction accuracy of polygenic scores can be found [here](prs_wt_finemapping.html)

10. Code for studying the effect of variant ascertainment schemes and effect size (re-)estimation can be found [here](biasvaccuracy_prsascertainment.html)

The [wiki page](https://github.com/Arslan-Zaidi/popstructure/wiki) contains additional information on how scripts were chained etc. for each analysis.

/optimize_migration_rate : scripts used to choose migration rate between demes such that the Fst and genomic inflation (on birthplace) matches that observed in the UK Biobank.

/simulating_genotypes : scripts to carry out msprime simulations, output the genotypes from these simulations, and run PCA on those genotypes. The folder is organized further by different demographic models used:
  - /grid : Simulations carried out for individuals sampled in a 6x6 square grid.
    - tau-9 : Individuals sampled in a 6x6 grid where all demes remain structured back in time (i.e. perpetual model)
    - tau100 : Individuals sampled in a 6x6 grid where all demes become panmictic 100 generations in the past.
  - /ukb : Simulations carried out for individuals sampled from 35 demes arranged according to the geographic structure in England and Wales.

/pca_plots : scripts used to generate scatter plots of eigenvectors from PCA carried out on simulated genotype data

/gwas : scripts used to simulate phenotypes, carry out GWAS, and construct polygenic scores.

/burden_msprime: scripts used to carry out burden tests
