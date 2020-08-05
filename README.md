# Population structure simulation and anlysis

This repository contains scripts and code to carry out the analyses in:

Demographic history impacts stratification in polygenic scores
Arslan A. Zaidi, Iain Mathieson

bioRxiv 2020.07.20.212530
doi: https://doi.org/10.1101/2020.07.20.212530

Please cite this paper if you use any of the code in your own analyses.

The following is a brief description of each folder. The [wiki page|https://github.com/Arslan-Zaidi/popstructure/wiki] contains additional information on how scripts were chained etc. for each analysis.

/optimize_migration_rate : scripts used to choose migration rate between demes such that the Fst and genomic inflation (on birthplace) matches that observed in the UK Biobank.

/simulating_genotypes : scripts to carry out msprime simulations, output the genotypes from these simulations, and run PCA on those genotypes. The folder is organized further by different demographic models used:
  - /grid : Simulations carried out for individuals sampled in a 6x6 square grid.
    - tau-9 : Individuals sampled in a 6x6 grid where all demes remain structured back in time (i.e. perpetual model)
    - tau100 : Individuals sampled in a 6x6 grid where all demes become panmictic 100 generations in the past.
  - /ukb : Simulations carried out for individuals sampled from 35 demes arranged according to the geographic structure in England and Wales.

/pca_plots : scripts used to generate scatter plots of eigenvectors from PCA carried out on simulated genotype data

/gwas : scripts used to simulate phenotypes, carry out GWAS, and construct polygenic scores.

/burden_msprime: scripts used to carry out burden tests

