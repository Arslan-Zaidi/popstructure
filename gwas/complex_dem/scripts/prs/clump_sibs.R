
#!/usr/bin/env Rscript

args=commandArgs(TRUE)

if( length(args) != 4){stop("Usage: <causal effects file> <prefix for glm.linear file (w/o phenotype or correction)> <phenotype name> <output file prefix>") }

library(data.table)
library(dplyr)
library(tidyr)

geffects_file=args[1]
gwas_file=args[2]
phenotype=args[3]
output_file_prefix=args[4]



#read in list of causal variants (true simulated effects)
causal=fread(geffects_file)
colnames(causal)=c("rsid","allele","esize")

causal=causal%>%
  separate(rsid,into=c("CHROM","POS","ref","alt"),sep="_",remove=F)

causal$POS=as.numeric(causal$POS)
causal$CHROM=as.numeric(causal$CHROM)

print("reading gwas files")
#gwas effects
gwas1=fread(gwas_file,fill=T)


colnames(gwas1) = c("CHROM","ID","POS","REF","A1","BETA","SE","P","CHI","CHI_P")

#function to get the effect size for the T allele
flip_effect = function(gwas_df,beta_colname){
  #gwas_df = gwas_df[ ID %in% causal$rsid, .(CHROM,POS,ID,A1,BETA,P)]
  gwas_df = gwas_df[ A1=="A", beta_colname := -BETA]
  gwas_df = gwas_df[ A1=="T", beta_colname := BETA]
  gwas_df = gwas_df[,.(CHROM, POS,ID,A1,beta_colname,P)]
  colnames(gwas_df)[5] = beta_colname
  return(gwas_df)
}

gwas1 = flip_effect(gwas1,beta_colname = "BETA1")

print("filtering causal variants")
#select effect sizes for causal variants
gwas1.1 = gwas1[ID%in%causal$rsid]
gwas1.1 = gwas1.1[,.(ID,A1,BETA1,P)]


print("filtering variants under a pvalue threshold")

#write function to select causal variants below some p-value threshold
fcausal_p = function(gwas_df,beta_colname,pvalue=5e-04){

    gwas_df[ P >= pvalue, (beta_colname) := 0]
    gwas_df[ P < pvalue, (beta_colname) := BETA1]
    return(gwas_df)
}

gwas1.2 = fcausal_p( gwas1.1, "BETA2" )



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

gwas1.red=fclump(gwas1,5e-04)
gwas1.red=gwas1.red[,flead(.SD),by=.(window_name)]
gwas1.red=gwas1.red[,.(ID,A1,BETA1)]

gwas1.red$A1="T"
colnames(gwas1.red)[3]="BETA3"

gwas.all = merge(gwas1.1[,.(ID,A1,BETA1)],
                 gwas1.2[,.(ID,BETA2)], by="ID",all=TRUE)
gwas.all = merge(gwas.all,
                 gwas1.red[,.(ID,BETA3)], by="ID",all=TRUE)

gwas.all[is.na(BETA1),BETA1:=0]
gwas.all[is.na(BETA2),BETA2:=0]
gwas.all[is.na(BETA3),BETA3:=0]

gwas.all$A1="T"


fwrite(gwas.all,
       paste(output_file_prefix,".all.betas",sep=""),
       col.names=F,row.names=F,quote=F,sep="\t")
