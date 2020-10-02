#!/usr/bin/env Rscript

#script that takes a set of variants ascertained in one sample set and
#re-estimates effects in a different sample set (finds from a gwas output file)

args = commandArgs(TRUE)
variants_file = args[1] #list of variants for which to get se
effects_file = args[2] #pre-computed GWAS results file with BETA and SE columns
output_dir = args[3] # output directory path
gwas_format = args[4] # gwas format

if( length(args) != 4 ){stop("Usage: get_se.R <variants list> <gwas effects> <output dir> <gwas_format (either 'plink' or 'sibling')>")}

library(ggplot2)
library(dplyr)
library(data.table)
library(rprojroot)
F = is_rstudio_project$make_fix_file()

print("reading variant list")
#read list of ascertained variants
variants = fread(F(variants_file))

#this list usually contains the following columns:
#1. effect size ID
#2. effect allele
#3. - n: effect sizes for different methods of correction or clumping (e.g. causal, lead)

print("reading gwas results")

#read gwas effects file
effects = fread(F(effects_file))

#assign column names based on whether the format is plink or sibling
if(gwas_format == "plink"){

colnames(effects) = c("CHROM","POS","ID","REF","ALT","A1","TEST","OBS_CT","BETA","SE","T_STAT","P")

}

if(gwas_format == "sibling"){

colnames(effects) = c("ID","POS","REF","A1","BETA","SE","P","CHI","CHI_P")

}

effects = effects[,.(ID,SE)] ##keep relevant columns


#merge with variant list
variants = merge(variants, effects, by.x="V1", by.y="ID",sort=FALSE)

variants[is.na(SE),SE:=0]

#write to file
fwrite(variants, F(output_dir),
       col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
