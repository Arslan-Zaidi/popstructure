

#!/usr/bin/env Rscript

args=commandArgs(TRUE)

if(length(args)<3){
  stop("At least three arguments required: [ldsc file] [gwas file] [output prefix]")
}else{

ldsc_file = args[1] #file with LD scores
gwas_file = args[2] #file with GWAS effect sizes
frq_file = args[3] # frequency of the variants
out_prefix = args[4] #prefix of output file
}

#load libraries
library(data.table)
library(ggplot2)
library(here)
library(cowplot)
library(dplyr)
library(tidyr)
library(readr)

#Load frequency file
print("loading input files")
freq100=fread(frq_file)
colnames(freq100)=c("CHROM","SNP","REF","ALT","ALT_FREQS","OBS_CT")

#calculate MAF from Alt allele frequency
freq100=freq100%>%
  mutate(MAF=case_when(ALT_FREQS<0.5~ALT_FREQS,
                       TRUE~1-ALT_FREQS))

#load LD scores
ldsc=fread((ldsc_file))
colnames(ldsc)=c("CHR","SNP","POS","LDScore")

#add MAF to LD scores
ldsc=merge(ldsc,freq100[,c("SNP","MAF")],by="SNP")

#bin SNPs by LD scores - helpful for plotting
ldsc=ldsc%>%
  filter(MAF>0.01)%>%
  mutate(ldbin=cut(LDScore,
                   breaks = quantile(LDScore,probs=seq(0,1,0.01)),
                   include.lowest = T))

#function to read GWAS results and process them
read_gwas100<-function(gwas_file){

  df=fread(gwas_file) # read GWAS file
  snp_column=which(colnames(df)%in%c("SNP","ID")) #check for column named 'SNP' or 'ID'
  snp_col_name=colnames(df)[snp_column]
  df=df[,c(snp_col_name,"BETA","SE","P"),with=F]

  df=merge(df,ldsc[,c("SNP","CHR","POS","LDScore","ldbin","MAF")],
           by.x=snp_col_name,by.y="SNP",all.y=T)

  #calculate regression weights and x2 stats from beta and se
  df=df%>%
    mutate(WT = 1/(1+LDScore),
           CHI = (BETA/SE)^2)

  return(df)
}

#read GWAS effect sizes
gwas_effects=read_gwas100(gwas_file)

#function to calculate intercept and slope
flm.pt=function(df){
  l1=lm(data=df, CHI~LDScore, weights = WT)
  s1=summary(l1)
  d1=data.table(Variable=c("Intercept","Slope"),
                Estimate=s1$coefficients[,1])
  return(d1)
}

print("Carrying out LDSC regression and generating SEs")
#calculate point estimate for intercept and slope
ldsc.pt=gwas_effects%>%
  do(flm.pt(.))

#save point estimates to file
fwrite(ldsc.pt,
  paste(out_prefix,
        ".ldsc.pt",
        sep=""),
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)
