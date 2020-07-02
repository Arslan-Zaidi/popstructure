

#!/usr/bin/env Rscript

args=commandArgs(TRUE)

if( length(args) != 3){stop("Usage: <causal effects file> <prefix for glm.linear file (w/o phenotype or correction)> <output file prefix>") }

library(data.table)
library(dplyr)
library(ggplot2)

gwas_file = args[1]
effects_file = args[2]
output_file=args[3]

gwas1 = fread(gwas_file)
colnames(gwas1) = c("CHROM","SNP","POS","REF","A1","BETA","SE","P","CHI","CHI_P")

# #function to get the effect size for the T allele
# flip_effect = function(gwas_df,beta_colname){
#   #gwas_df = gwas_df[ ID %in% causal$rsid, .(CHROM,POS,ID,A1,BETA,P)]
#   gwas_df = gwas_df[ A1=="A", beta_colname := -BETA]
#   gwas_df = gwas_df[ A1=="T", beta_colname := BETA]
#   gwas_df = gwas_df[,.(CHROM, POS,SNP,A1,beta_colname,P,BETA)]
#   colnames(gwas_df)[5] = beta_colname
#   return(gwas_df)
# }
#
# gwas1 = flip_effect(gwas1,beta_colname = "BETA1")

effects.nc.asc = fread(effects_file)
colnames(effects.nc.asc) = c("SNP","A1","pcs0","cm","re","cmre")

gwas.red = gwas1[SNP%in%effects.nc.asc$SNP]
gwas.red = gwas.red[,.(SNP,A1,BETA)]
gwas.red[is.na(BETA),BETA:=0]
gwas.red[,BETA1:=BETA]
gwas.red[,BETA2:=BETA]
gwas.red[,BETA3:=BETA]
gwas.red[,BETA4:=BETA]

list_pcs0=effects.nc.asc$SNP[which(effects.nc.asc$pcs0==0)]
list_cm=effects.nc.asc$SNP[which(effects.nc.asc$cm==0)]
list_re=effects.nc.asc$SNP[which(effects.nc.asc$re==0)]
list_cmre=effects.nc.asc$SNP[which(effects.nc.asc$cmre==0)]

gwas.red[ SNP%in%list_pcs0, BETA1:=0 ]
gwas.red[ SNP%in%list_cm, BETA2:=0 ]
gwas.red[ SNP%in%list_re, BETA3:=0 ]
gwas.red[ SNP%in%list_cmre, BETA4:=0 ]

gwas.red = gwas.red[,.(SNP,A1,BETA1,BETA2,BETA3,BETA4)]

fwrite(gwas.red,
       output_file,
       col.names=FALSE,
       row.names=FALSE,
       quote = FALSE,
       sep="\t")
