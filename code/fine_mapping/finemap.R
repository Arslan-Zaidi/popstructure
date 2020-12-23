#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(rprojroot)
library(susieR)
library(matrixcalc)
library(Matrix)

F = is_rstudio_project$make_fix_file()

args = commandArgs(TRUE)
if( length(args) != 4 ){ stop("Usage: finemap.R <gwas_file> <pvalue cutoff> <path to ld files> <output_file> ") }

gwas_file = args[1]
pvalue = args[2]
ldpath = args[3]
output_file = args[4]

print("reading GWAS results")

gwas = fread( F(gwas_file) )
colnames(gwas)[1] = "CHROM"

#function to get the effect size for the T allele
flip_effect = function(gwas_df,beta_colname){
  #gwas_df = gwas_df[ ID %in% causal$rsid, .(CHROM,POS,ID,A1,BETA,P)]
  gwas_df = gwas_df[ A1=="A", beta_colname := -BETA]
  gwas_df = gwas_df[ A1=="T", beta_colname := BETA]
  gwas_df$A1="T"
  gwas_df = gwas_df[,.(CHROM,POS,ID,A1,beta_colname,P,T_STAT)]
  colnames(gwas_df)[5] = beta_colname
  return(gwas_df)
}

gwas = flip_effect(gwas,beta_colname = "BETA1")

fclump=function(gwas_df, pcutoff = pvalue){
  if ( missing(pcutoff) ){
    df.red = gwas_df
  }else{
    df.red = gwas_df[ P < pcutoff ]
  }

  for(i in 1:100){
    if(i == 1){
      start = 1
      stop = (start + 1e5) - 1
    }else{
      start=((i-1)*1e5) +1
      stop=start + 1e5 -1
    }
    #print(i)
    df.red[ (POS>=start & POS<=stop) , window:=i]
    df.red[(POS>=start & POS<=stop) ,start.bp:=start]
    df.red[(POS>=start & POS<=stop) ,stop.bp:=stop]
  }
  df.red[,window_name:= paste(CHROM,window,sep="_") ]
  return(df.red)

}

gwas = fclump(gwas)
windows = gwas%>%
  filter(P < pvalue)%>%
  distinct(window_name)

gwas.red = gwas[window_name %in% windows$window_name]

#function to fine-map each window
fine.map = function(dat){
  chrom = dat$CHROM[1]
  start = format(dat$start[1],scientific = FALSE)
  stop = format(dat$stop[1],scientific = FALSE)

  print( paste( "fine-mapping :- ", chrom, ":", start, "-", stop, sep = "" ))
  ldfile = paste(ldpath,"/ld_chr", chrom, "_", start, "_", stop, ".ld.gz",
                  sep = "")

  #load ld file
  R = as.matrix(fread(F(ldfile)))

  #ensure ld matrix is positive definite
  npd = nearPD(R)
  R2 = as.matrix(npd$mat)

  #run susie
  fitted_rss = susie_rss(dat$T_STAT,
                         R2, L=1, check_z = FALSE,
                         estimate_prior_variance = TRUE,
                         estimate_residual_variance = TRUE)

  #get variant with highest posterior probability
  pips = fitted_rss$pip
  lead_index = which(pips == max(pips))

  #sometimes all variants have pip=0.
  #we don't want to save all
  if( pips[lead_index] == 0 ){
    index_snp = data.table(ID = "none", 
                           A1 = "none", 
                           BETA1 = 0, 
                           pip = 0)
  }else{
    #sometimes multiple variants have the same pip
    #select one at random if so
    sample.vec <- function(x, ...) x[sample(length(x), ...)]
    lead_index = sample.vec(lead_index, 1) 
    index_snp = dat[lead_index, .(ID,A1,BETA1)]
    index_snp$pip = pips[lead_index]}

  #check if credible set is empty
  cs = summary(fitted_rss)$cs

  if( is.null(cs) == "FALSE" ){
    pindicator = 1}else{
      pindicator = 0}

  index_snp$in_cs = pindicator
  
  #now calculate marginal beta using original effect sizes
  #we're going to assume effect in the same direction as the lead snp
  marg.b1 = t(abs(dat$BETA1)) %*% pips
  if(index_snp$BETA1 < 0 ){
    marg.b1 = -marg.b1}else{
      marg.b1 = marg.b1
    }
  
  index_snp$marginal_b = marg.b1
  
  return(index_snp)
}

est.betas = gwas.red[, fine.map(.SD) , by = .(window_name)]

fwrite( est.betas %>% select( ID, A1, BETA1, pip, in_cs, marginal_b, window_name ) ,
        F(output_file), sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
