
args=commandArgs(TRUE)

if(length(args)!=2){stop("usage: Rscript subsample4qq.R [gwas_file] [output_file]")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
}))

gwas_file=args[1]
output=args[2]


freq[,fcat:=cut(ALT_FREQS, breaks=c(0,0.01,0.05,1), labels=c("rare","medium","common"))]
fcat_df=freq$fcat

snp.counts=freq%>%group_by(fcat)%>%count
snp.cm=as.numeric(snp.counts[2,2])
snp.re=as.numeric(snp.counts[1,2])

rm(freq)

print("reading GWAS results")

#function that:
#1. reads in GWAS association
#2. calculates expected -log10pvalues
#3. subsamples for plotting
read_gwas<-function(){
    df=fread(gwas_file)

  colnames(df)=c("CHROM","SNP","position","A1","A2","ALT_FREQS","BETA","SE","P")
  df=df[,c("ID","P")]

  df<-cbind(df,fcat_df)
  colnames(df)[ncol(df)]="fcat"

  df=df[fcat%in%c("common","rare"),]


  #calculating expected Pvalue distribution and lambda
  df=df[order(P),.SD,by="fcat"]

  df[ , ix := 1:.N, by = fcat]

  #df[,c("exp.p","obs.chi"):=list(ppoints(length(ID)), qchisq(P,df=1,lower.tail = F)), by="fcat"]
  #df[,"exp.chi":=qchisq(exp.p,df=1,lower.tail = F)]
  #df[,"lambda":=obs.chi/exp.chi]
  #df=df[,chi.percentile:=(length(exp.p):1)/length(exp.p),by="fcat"]

  #df[,lower.ci:=qbeta(0.025,shape1=ix,shape2 = max(ix)-ix),by=fcat]
  #df[,upper.ci:=qbeta(0.975,shape1=ix,shape2 = max(ix)-ix),by=fcat]

  #reducing the number of rows for plotting
  #keep the first 1000 SNPs (low-p tails) and subsample from the rest (only for rare variants, keep all common)
  #calculate step size so that only 100k SNPs are sampled
  #function to sample
  f.subsample=function(df1){
    if(nrow(df1)>(1e4+1000)){
      step_size=round(nrow(df1)/1e5)
      df2 = rbind(df1[c(1:1000),],
                      df1[seq(1001,nrow(df1),step_size),])
    }else{df2=df1}
    return(df2)
  }
  df.common = f.subsample(df[fcat == "common",])
  df.rare = f.subsample(df[fcat == "rare",])
  df = rbind(df.common,df.rare)

  return(df)
}

gwas_results = read_gwas()

#gwas_results = gwas_results[,.(fcat,exp.p,P,chi.percentile,lower.ci,upper.ci,lambda),]

fwrite(gwas_results,
      output,
      col.names=T,row.names=F,quote=F,sep="\t")
