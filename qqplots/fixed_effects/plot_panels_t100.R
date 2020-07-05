
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)



#function to summarize
fsummarize=function(phenotype="smooth",correction="pcs0"){
  dat = fread(paste("gwas/grid/genotypes/tau100/ss500/train/gwas_results/fixed_effects/noge/genos_gridt100_",
                    phenotype,".",
                    correction,
                    ".all.txt.gz",sep=""))
  colnames(dat) = c("fcat","exp.p","P","chi.percentile","lower.ci","upper.ci","lambda")
  dat = dat[, lapply(.SD,mean), by = .(fcat,chi.percentile)]
  
  return(dat)
}

dat1=fsummarize("smooth","pcs0")
dat2=fsummarize("smooth","cm")
dat3=fsummarize("smooth","re")

max.lambda=max(sapply(list(dat1,dat2,dat3),
                      function(x){
                        max( x[ which(x$chi.percentile>0.999), "lambda"])
                      }))

fplot=function(dat,title="No correction"){
  
  plt1<-ggplot(data=dat)+
    geom_abline(intercept=0,slope=1,color="black")+
    geom_ribbon(aes(x=-log10(exp.p),
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
         title=title)+
    xlim(c(0,8.1))+
    ylim(c(0,8.1))
  
  plt.inset=ggplot()+
    geom_line(data=dat[chi.percentile>0.999,],
              aes(chi.percentile,
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

plt1=fplot(dat1,"No correction")
plt2=fplot(dat2,"Common PCA")
plt3=fplot(dat3,"Rare PCA")

plt_combined.sm=plt1+plt2+plt3


dat1=fsummarize("sharp","pcs0")
dat2=fsummarize("sharp","cm")
dat3=fsummarize("sharp","re")

max.lambda=max(sapply(list(dat1,dat2,dat3),
                      function(x){
                        max( x[ which(x$chi.percentile>0.999), "lambda"])
                      }))

plt1=fplot(dat1,"No correction")
plt2=fplot(dat2,"Common PCA")
plt3=fplot(dat3,"Rare PCA")

plt_combined.shp=plt1+plt2+plt3

plt_combined.both=plt_combined.sm/plt_combined.shp

ggsave("plots/qqplots/plt_qq_t100_07032020.pdf",
       plt_combined.both,
       height=12,
       width=14,
       units="cm")

