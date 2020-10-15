
#!/usr/bin/env Rscript

args=commandArgs(TRUE)

if( length(args) != 3){stop("Usage: <causal effects file> <prefix for glm.linear file (w/o phenotype or correction)> <phenotype name> <output file prefix>") }

library(data.table)
library(dplyr)
library(tidyr)

gwas_file_prefix=args[1]
phenotype=args[2]
output_file_prefix=args[3]


print("reading gwas files")
#gwas effects
gwas1=fread(paste(gwas_file_prefix,".pcs0.",phenotype,".glm.linear.gz",sep=""),fill=T)

colnames(gwas1)[1] = "CHROM"

#function to get the effect size for the T allele
flip_effect = function(gwas_df,beta_colname){
  #gwas_df = gwas_df[ ID %in% causal$rsid, .(CHROM,POS,ID,A1,BETA,P)]
  gwas_df = gwas_df[ A1=="A", beta_colname := -BETA]
  gwas_df = gwas_df[ A1=="T", beta_colname := BETA]
  gwas_df$A1="T"
  gwas_df = gwas_df[,.(CHROM,POS,ID,A1,beta_colname,P)]
  colnames(gwas_df)[5] = beta_colname
  return(gwas_df)
}

gwas1 = flip_effect(gwas1,beta_colname = "BETA1")

print("ld clumping")
#write function to assign each SNP to window, then find a single hit within each
fclump=function(gwas_df,pcutoff=5e-04){
 if(missing(pcutoff)){
   df.red=gwas_df
 }else{
   df.red = gwas_df[P<pcutoff]
 }


 for(i in 1:101){
   if(i == 1){
     start=0
     stop=start+1e5
   }else{

   start=((i-1)*1e5) +1
   stop=start + 1e5 -1
   }

   df.red[ (POS>=start &POS<stop) , window:=i]
 }
 df.red[,window_name:= paste(CHROM,window,sep="_") ]
 return(df.red)

}

flead=function(df){
  df.red = df[P == min(P)]
 return(df.red)
}

#causal=fclump(causal)
gwas1.red = fclump(gwas1,5e-04)
gwas1.red = gwas1.red[,flead(.SD),by=.(window_name)]
gwas1.red = unique(gwas1.red, by ="window_name")
gwas1.red = gwas1.red[,.(ID,A1,BETA1)]

gwas1.red$A1 = "T"

gwas1.red[is.na(BETA1),BETA1:=0]

fwrite(gwas1.red,
       paste(output_file_prefix,".nc.betas",sep=""),
       col.names=F,row.names=F,quote=F,sep="\t")
