

#!/usr/bin/env Rscript

args=commandArgs(TRUE)

if(length(args)<4){
  stop("At least three arguments required: [ldsc file] [gwas file] [freq_file] [output prefix]")
}else{

ldsc_file = args[1] #file with LD scores
gwas_file = args[2] #file with GWAS effect sizes
frq_file = args[3]
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

#sort gwas effects by chromosome and position
gwas_effects=gwas_effects%>%arrange(CHR,POS)

#function to calculate SE for intercept and slope using block Jacknife
#currently supports 1,000 bootstraps
flm.se=function(df,nboot){
  mat=matrix(NA,nrow=nboot,ncol=2)
  #total number of SNPs
  nsnps.total=nrow(df)
  for(i in 1:nboot){
    #Select a window to drop
    rand.snp=sample(nsnps.total,1) # first, select random SNP

    #construct window around this SNP
    if((rand.snp+2000)<=nrow(df)){
      start=rand.snp
      stop=start+2000
    }else{
      start=rand.snp-2000
      stop=rand.snp
    }

    df2=df[-(start:stop),] #drop this window

    #calculate LDSC intercept and slope
    l1=lm(data=df2, CHI~LDScore, weights = WT)
    s1=summary(l1)
    mat[i,1]=s1$coefficients[1,1]
    mat[i,2]=s1$coefficients[2,1]
  }

  #calculate 95%CI for intercept and slope
  lower.intercept=quantile(mat[,1],probs=0.025)
  upper.intercept=quantile(mat[,1],probs=0.975)

  lower.slope=quantile(mat[,2],probs=0.025)
  upper.slope=quantile(mat[,2],probs=0.975)

  #put this into dataframe
  d1=data.table(variable=c("intercept","slope"),
                lower=c(lower.intercept,lower.slope),
                upper=c(upper.intercept,upper.slope))
  return(d1)
}

print("plotting LDSC regression")
#function to generate weighted regression line of best fit for plotting
lmpredict=function(df){
  l1=lm(data=df, CHI~LDScore, weights = WT)
  df$PCHI=predict(l1)
  return(df)
}

#generate weighted best fit line
gwas_effects=gwas_effects%>%
  do(lmpredict(.))

#calculate mean LD score and x2 statistic per LDSC bin for plotting
gwas.binned=gwas_effects%>%
  group_by(ldbin)%>%
  summarize(MAF=mean(MAF),
            WT=mean(WT),
            CHI=mean(CHI),
            LDScore=mean(LDScore),
            PCHI=mean(PCHI))

f.ldsc.plt=function(binned.df,pt.estimates){
  ggplot()+
    geom_point(data=binned.df,
               aes(LDScore,CHI,color=WT),
               show.legend=T)+
    theme_bw()+
    geom_line(data=binned.df,
              aes(LDScore,PCHI),
              color="black",
              linetype="dashed")+
    geom_hline(yintercept=1,color="red")+
    theme(panel.grid = element_blank())+
    scale_color_gradient(low="blue",high="red")+
    geom_hline(data=pt.estimates%>%filter(Variable=="Intercept"),
               aes(yintercept=Estimate),
               color="black",
               linetype="dashed")+
    labs(x="LD score bin",
         y=bquote("Mean"~chi^2),
         color="Weight")}

#plot
plt=f.ldsc.plt(gwas.binned,ldsc.pt)

#save plot
ggsave(
  paste(out_prefix,
        ".pdf",
        sep=""),
       plt,
       height=4,
       width=2)


#save ggplot object for later manipulation
write_rds(plt,
	paste(out_prefix,".rds",sep=""),compress="none")

#calculate SE for LDSC intercept and slope
ldsc.se=gwas_effects%>%
  do(flm.se(.,nboot=1000))

#ldsc estimates (both point and SE)
ldsc.est=cbind(ldsc.pt,ldsc.se[,c(2,3)])

#save point estimates and SE to file
fwrite(ldsc.est,
  paste(out_prefix,
        ".estimates",
        sep=""),
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)
