#script that takes a set of variants ascertained in one sample set and
#re-estimates effects in a different sample set (finds from a gwas output file)

args = commandArgs(TRUE)
variants_file = args[1] #list of LD-clumped variants
effects_file = args[2] #pre-computed GWAS results file with BETA column
output_dir = args[3] # output directory path

if( length(args) != 3 ){stop("Usage: Rscript plot_prs_all.R <variants list> <gwas effects> <output dir>")}

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

mvariants = reshape2::melt(variants,id.vars=c("V1","V2"))

print("reading gwas results")
#read gwas file
effects = fread(F(effects_file))
#keep relevant columns
effects = effects[,.(ID,A1,BETA)]

#merge with variant list
mvariants = merge(mvariants, effects, by.x="V1", by.y="ID",sort=FALSE)

mvariants = mvariants%>%
  mutate(esize = case_when(value==0~0,
                           TRUE~BETA))

#dcast
mvariants = reshape2::dcast(mvariants,V1+V2~variable,value.var="esize")

#write to file
fwrite(mvariants, F(output_dir),
       col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
