
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


colnames(gwas1)[1]=colnames(gwas2)[1]=colnames(gwas3)[1]="CHROM"

print("filtering causal variants")
#select effect sizes for causal variants
gwas1.1 = gwas1%>%
    filter(ID %in% causal$rsid)%>%
    mutate(BETA1 = case_when(A1=="A"~-BETA,
                             TRUE~BETA))

gwas2.1 = gwas2%>%
    filter(ID %in% causal$rsid)%>%
  mutate(BETA2 = case_when(A1=="A"~-BETA,
                           TRUE~BETA))

gwas3.1 = gwas3%>%
    filter(ID %in% causal$rsid)%>%
  mutate(BETA3 = case_when(A1=="A"~-BETA,
                           TRUE~BETA))

gwas4.1 = gwas4%>%
   filter(ID %in% causal$rsid)%>%
 mutate(BETA4 = case_when(A1=="A"~-BETA,
                          TRUE~BETA))

gwas.causal=cbind( gwas1.1[,c("ID","A1","BETA1")],
                    gwas2.1[,c("BETA2")],
                    gwas3.1[,c("BETA3")],
                    gwas4.1[,c("BETA4")])

fwrite(gwas.causal,
       paste(output_file_prefix,".c.betas",sep=""),
       col.names=F,row.names=F,quote=F,sep="\t")

print("filtering variants under a pvalue threshold")

#write function to select causal variants below some p-value threshold
fcausal_p = function(df,pvalue=5e-04){

    df=df%>%
    mutate(BETA = case_when(P > pvalue ~ 0,
                          TRUE ~ BETA))
    return(df)
}

gwas1.2 = fcausal_p( gwas1.1 )
gwas2.2 = fcausal_p( gwas2.1 )
gwas3.2 = fcausal_p( gwas3.1 )
gwas4.2 = fcausal_p( gwas4.1 )


gwas.causal.p = cbind( gwas1.2[,c("ID","A1","BETA")],
                    gwas2.2[,c("BETA")],
                    gwas3.2[,c("BETA")],
                    gwas4.2[,c("BETA")])


fwrite(gwas.causal.p,
       paste(output_file_prefix, ".c.p.betas" , sep=""),
       col.names=F, row.names=F , quote=F , sep="\t")

print("ld clumping")
#write function to assign each SNP to window, then find a single hit within each
fclump=function(df,pcutoff){
 if(missing(pcutoff)){
   df.red=df
 }else{
   df.red=df%>%
     filter(P<pcutoff)
 }

 df.red$window=NA
 for(i in 1:101){
   start=((i-1)*1e5 - 5e4) +1
   stop=start+1e5
   df.red$window[which((df.red$POS>=start) & (df.red$POS<stop))]=i
 }

 df.red$window_name=paste(df.red$CHROM,df.red$window,sep="_")

 return(df.red)

}

flead=function(df){
 df.red=df%>%
   filter(P==min(P))
 return(df.red)
}

causal=fclump(causal)

gwas1.red=fclump(gwas1,
                 0.0005)%>%
 group_by(window_name)%>%
 flead(.)%>%
 ungroup()%>%
 select(ID,A1,BETA)

gwas2.red=fclump(gwas2,
                0.0005)%>%
  group_by(window_name)%>%
  flead(.)%>%
  ungroup()%>%
  select(ID,A1,BETA)

gwas3.red=fclump(gwas3,
                 0.0005)%>%
 group_by(window_name)%>%
 flead(.)%>%
 ungroup()%>%
 select(ID,A1,BETA)

gwas4.red=fclump(gwas4,
                  0.0005)%>%
  group_by(window_name)%>%
  flead(.)%>%
  ungroup()%>%
  select(ID,A1,BETA)

#flip effect size if A1 == "A"
#this way, the effects are consistent across all gwas1.reds and they can be cbinded

gwas1.red = gwas1.red %>%
  mutate(BETA1 = case_when(A1 == "A"~ -BETA,
                           A1 == "T"~ BETA))
gwas2.red = gwas2.red %>%
  mutate(BETA2 = case_when(A1 == "A"~ -BETA,
                           A1 == "T"~ BETA))
gwas3.red = gwas3.red %>%
  mutate(BETA3 = case_when(A1 == "A"~ -BETA,
                           A1 == "T"~ BETA))

gwas4.red = gwas4.red %>%
 mutate(BETA4 = case_when(A1 == "A"~ -BETA,
                          A1 == "T"~ BETA))

gwas.red = merge(gwas1.red[,c("ID","A1","BETA1")], gwas2.red[,c("ID","BETA2")], by="ID", all=TRUE)
gwas.red = merge(gwas.red, gwas3.red[,c("ID","BETA3")] , by="ID", all=TRUE)
gwas.red = merge(gwas.red, gwas4.red[,c("ID","BETA4")] , by="ID", all=TRUE)

gwas.red$A1="T"

gwas.red = gwas.red%>%
            replace_na(list( BETA1=0 , BETA2 = 0, BETA3 = 0, BETA4=0))

fwrite(gwas.red,
       paste(output_file_prefix,".nc.betas",sep=""),
       col.names=F,row.names=F,quote=F,sep="\t")
