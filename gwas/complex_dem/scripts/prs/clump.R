
#!/usr/bin/env Rscript

args=commandArgs(TRUE)

if( length(args) != 4){stop("Usage: <causal effects file> <prefix for glm.linear file (w/o phenotype or correction)> <phenotype name> <output file prefix>") }

library(data.table)
library(dplyr)
library(tidyr)

geffects_file=args[1]
gwas_file_prefix=args[2]
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
gwas1=fread(paste(gwas_file_prefix,".pcs0.",phenotype,".glm.linear.gz",sep=""),fill=T)
gwas2=fread(paste(gwas_file_prefix,".cm.",phenotype,".glm.linear.gz",sep=""),fill=T)
gwas3=fread(paste(gwas_file_prefix,".re.",phenotype,".glm.linear.gz",sep=""),fill=T)
gwas4=fread(paste(gwas_file_prefix,".cmre.",phenotype,".glm.linear.gz",sep=""),fill=T)


colnames(gwas1)[1]=colnames(gwas2)[1]=colnames(gwas3)[1]=colnames(gwas4)[1]="CHROM"

#function to get the effect size for the T allele
flip_effect = function(gwas_df,beta_colname){
  #gwas_df = gwas_df[ ID %in% causal$rsid, .(CHROM,POS,ID,A1,BETA,P)]
  gwas_df = gwas_df[ A1=="A", beta_colname := -BETA]
  gwas_df = gwas_df[ A1=="T", beta_colname := BETA]
  gwas_df$A1=="T"
  gwas_df = gwas_df[,.(CHROM,POS,ID,A1,beta_colname,P)]
  colnames(gwas_df)[5] = beta_colname
  return(gwas_df)
}

gwas1 = flip_effect(gwas1,beta_colname = "BETA1")
gwas2 = flip_effect(gwas2,beta_colname = "BETA2")
gwas3 = flip_effect(gwas3,beta_colname = "BETA3")
gwas4 = flip_effect(gwas4,beta_colname = "BETA4")

print("filtering causal variants")
#select effect sizes for causal variants
gwas1.1 = gwas1[ID%in%causal$rsid]
gwas2.1 = gwas2[ID%in%causal$rsid]
gwas3.1 = gwas3[ID%in%causal$rsid]
gwas4.1 = gwas4[ID%in%causal$rsid]

gwas.causal=cbind( gwas1.1[,c("ID","A1","BETA1")],
                    gwas2.1[,c("BETA2")],
                    gwas3.1[,c("BETA3")],
                    gwas4.1[,c("BETA4")])

fwrite(gwas.causal,
       paste(output_file_prefix,".c.betas",sep=""),
       col.names=F,row.names=F,quote=F,sep="\t")

print("filtering variants under a pvalue threshold")

#write function to select causal variants below some p-value threshold
fcausal_p = function(gwas_df,beta_colname,pvalue=5e-04){
    
    gwas_df[ P > pvalue, (beta_colname) := 0]
    return(gwas_df)
}

gwas1.2 = fcausal_p( gwas1.1, "BETA1" )
gwas2.2 = fcausal_p( gwas2.1, "BETA2" )
gwas3.2 = fcausal_p( gwas3.1, "BETA3" )
gwas4.2 = fcausal_p( gwas4.1, "BETA4" )


gwas.causal.p = cbind( gwas1.2[,c("ID","A1","BETA1")],
                    gwas2.2[,c("BETA2")],
                    gwas3.2[,c("BETA3")],
                    gwas4.2[,c("BETA4")])


fwrite(gwas.causal.p,
       paste(output_file_prefix, ".c.p.betas" , sep=""),
       col.names=F, row.names=F , quote=F , sep="\t")

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

gwas2.red=fclump(gwas2,5e-04)
gwas2.red=gwas2.red[,flead(.SD),by=.(window_name)]
gwas2.red=gwas2.red[,.(ID,A1,BETA2)]

gwas3.red=fclump(gwas3,5e-04)
gwas3.red=gwas3.red[,flead(.SD),by=.(window_name)]
gwas3.red=gwas3.red[,.(ID,A1,BETA3)]

gwas4.red=fclump(gwas4,5e-04)
gwas4.red=gwas4.red[,flead(.SD),by=.(window_name)]
gwas4.red=gwas4.red[,.(ID,A1,BETA4)]

gwas.red = merge(gwas1.red[,c("ID","A1","BETA1")], gwas2.red[,c("ID","BETA2")], by="ID", all=TRUE)
gwas.red = merge(gwas.red, gwas3.red[,c("ID","BETA3")] , by="ID", all=TRUE)
gwas.red = merge(gwas.red, gwas4.red[,c("ID","BETA4")] , by="ID", all=TRUE)

gwas.red$A1="T"

gwas.red[is.na(BETA1),BETA1:=0]
gwas.red[is.na(BETA2),BETA2:=0]
gwas.red[is.na(BETA3),BETA3:=0]
gwas.red[is.na(BETA4),BETA4:=0]

fwrite(gwas.red,
       paste(output_file_prefix,".nc.betas",sep=""),
       col.names=F,row.names=F,quote=F,sep="\t")
