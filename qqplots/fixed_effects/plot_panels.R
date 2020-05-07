
args=commandArgs(TRUE)

if(length(args)!=3){stop("usage: Rscript plot_panels.R [phenotype] [tau] [output prefix]")} 

library(data.table)
library(ggplot2)
library(here)
library(cowplot)
library(patchwork)
library(dplyr)
library(tidyr)

phenotype=args[1]
tau=args[2]
output=args[3]

if(tau==100){migrate=0.05}
if(tau==-9){migrate=0.07}

genotype_dir=paste("gwas/grid/genotypes/tau",tau,"/ss250/gwas_results/",sep="")

print("Reading allele frequency information") ##########################

freq100<-fread(here(
  paste("gwas/grid/genotypes/tau",tau,
        "/ss250/genotypes/genos_grid_d36_m",
        migrate,
        "_s250_t",tau,"_gl200.frq.afreq",sep="")))
freq100[,fcat:=cut(ALT_FREQS, breaks=c(0,0.01,0.05,1), labels=c("rare","medium","common"))]
fcat_df=freq100$fcat

snp.counts=freq100%>%group_by(fcat)%>%count
snp.cm=as.numeric(snp.counts[2,2])
snp.re=as.numeric(snp.counts[1,2])

rm(freq100)

rlunif <- function(n, min, max, base=exp(1)) {
  if (mode(n) != "numeric")
    stop("'n' must be a non-empty numeric vector")
  if (any(missing(min), missing(max)))
    stop("'min' and 'max' not provided, without default.\n")
  ifelse(base == exp(1),
         return(exp(runif(n, log(min, base), log(max, base)))),
         return(base ^ (runif(n, log(min, base), log(max, base)))))
}

#generate CI for -log10(P) vaues
lambda1=data.table(common=c(1:1000,round(log(rlunif(4000,min=exp(1e-08),max=exp(0.99)))*snp.cm)),
                   rare=c(1:1000,round(log(rlunif(4000,min=exp(1e-08),max=exp(0.99)))*snp.re)))

mlambda=melt(lambda1)
colnames(mlambda)=c("fcat","ix")

mlambda=mlambda%>%
  group_by(fcat)%>%
  arrange(ix)%>%
  mutate(exp.p=ix/max(ix),
         lower.ci=qbeta(0.025,shape1=ix,shape2 = max(ix)-ix),
         upper.ci=qbeta(0.975,shape1=ix,shape2 = max(ix)-ix))

mlambda$fcat=factor(mlambda$fcat,levels=c("rare","common"))

#############################################################

print("reading GWAS results")

#function that:
#1. reads in GWAS association
#2. calculates expected -log10pvalues
#3. subsamples for plotting
read_gwas100<-function(pcs="pcs0"){
  #pcs can either be "pcs0","cm", or "re"
    df=fread(here(
      paste(genotype_dir,
            "gwas_grid_d36_m",
            migrate,"_s250_t",
            tau,
            "_gl200_e2_",pcs,".all.",phenotype,".glm.linear",
        sep="")))
    
  colnames(df)[1]="CHROM"
  df=df[,c("ID","P")]
  
  df<-cbind(df,fcat_df)
  colnames(df)[ncol(df)]="fcat"
  
  df=df[fcat%in%c("common","rare"),]

  
  #calculating expected Pvalue distribution and lambda
  df=df[order(P),.SD,by="fcat"]
  df[,c("exp.p","obs.chi"):=list(ppoints(length(ID)), qchisq(P,df=1,lower.tail = F)), by="fcat"]
  df[,"exp.chi":=qchisq(exp.p,df=1,lower.tail = F)]
  df[,"lambda":=obs.chi/exp.chi]
  df=df[,chi.percentile:=(length(exp.p):1)/length(exp.p),by="fcat"]
  
  #reducing the number of rows for plotting
  #keep the first 1000 SNPs (low-p tails) and subsample from the rest (only for rare variants, keep all common)
  df.common=df[fcat=="common",]
  df=df[fcat=="rare",]
  df=rbind(df.common,
           df[c(1:1000),],
           df[seq(1001,nrow(df),1000)]
  )
  
  return(df)
}

sm_fcommon_pc0_t100=read_gwas100(pcs = "pcs0")
sm_fcommon_pc100_t100=read_gwas100(pcs = "cm")
sm_frare_pc100_t100=read_gwas100(pcs = "re")

#get max lambda value for plotting
max.lambda=max(sapply(list(sm_fcommon_pc0_t100,sm_fcommon_pc100_t100,sm_frare_pc100_t100),
       function(x){
         max( x[ which(x$chi.percentile>0.999), "lambda"])
         }))

print("Making plots")

##function for plotting
fplt=function(df,
              tit=element_blank()){
  
  plt1<-ggplot(data=df)+
    geom_ribbon(data=mlambda,aes(x=-log10(exp.p),
                                 ymin=-log10(lower.ci),
                                 ymax=-log10(upper.ci),
                                 fill=fcat),
                alpha=0.2,
                show.legend = F)+
    geom_line(aes(-log10(exp.p),
                  -log10(P),
                  color=fcat),
              size=0.7,
              alpha=0.5,
              show.legend = F)+
    geom_abline(intercept=0,slope=1,color="black")+
    scale_color_manual(values=c("#ff7f00","#377eb8"))+
    scale_fill_manual(values=c("#ff7f00","#377eb8"))+
    theme_bw()+
    theme(panel.grid=element_blank(),
          axis.text=element_text(size=10),
          axis.title=element_blank(),
          plot.title = element_text(hjust=0.5),
          plot.background = element_blank(),
          plot.margin = unit(rep(0.5,4), "pt"))+
    labs(color="Freq.",
         title=tit)+
    xlim(c(0,8))+
    ylim(c(0,8))
  
  plt.inset=ggplot()+
    geom_line(data=df,aes(chi.percentile,
                          lambda,
                          color=fcat),
              show.legend = F,
              size=0.5)+
    annotate(geom="text",
             x=0.9993,
             y=0.9*max.lambda,
             label="lambda[p]",parse=T)+
    theme_bw()+
    theme(panel.grid.major.x = element_blank(),
          legend.position="none",
          axis.title=element_blank(),
          panel.grid=element_blank(),
          plot.background = element_blank(),
          axis.text.x = element_text(hjust=0,size=9),
          axis.text.y = element_text(size=9))+
    scale_x_log10(limits=c(0.999,1),
                  breaks=c(0.999,1),
                  labels=c("0.999","1"),
                  position="top")+
    scale_y_continuous(limits=c(0.99, round(max.lambda,2)),
                       breaks=c(1, round(max.lambda,2)),
                       position="left")+
    labs(x="p")+
    scale_color_manual(values=c("#ff7f00","#377eb8"))
  
  plt.wt.inset<-ggdraw(plt1) +
    draw_plot(plt.inset, x=0.3, y=0.08, height=0.4,width=0.7)
  
  return(plt.wt.inset)
  
}

plt_sm_fcommon_pc0_t100=fplt(sm_fcommon_pc0_t100, 
                             tit="No correction")

plt_sm_fcommon_pc100_t100=fplt(sm_fcommon_pc100_t100, 
                               tit="Common PCA")

plt_sm_frare_pc100_t100=fplt(sm_frare_pc100_t100, 
                             tit="Rare PCA")

plt_grid=plt_sm_fcommon_pc0_t100 +
  plt_sm_fcommon_pc100_t100 +
  plt_sm_frare_pc100_t100

ggsave(here(paste("analyses/qqplots/fixed_effects/",
                  output,
                  ".png",sep="")),
       plt_grid,
       height = 6,
       width = 14,
       units = "cm")

ggsave(here(paste("analyses/qqplots/fixed_effects/",
                  output,
                  ".pdf",sep="")),
       plt_grid,
       height = 6,
       width = 14,
       units = "cm")

