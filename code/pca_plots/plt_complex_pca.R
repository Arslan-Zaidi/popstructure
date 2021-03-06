

library(data.table)
library(dplyr)
library(ggplot2)
library(here)
library(patchwork)


cm = fread(here("gwas/complex_dem/train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.cmpruned.pca.eigenvec"))
re = fread(here("gwas/complex_dem/train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.repruned.pca.eigenvec"))

pop = fread(here("gwas/complex_dem/genos_complex_l1e7_ss500_m0.07.train.pop"))

#define bivariate color scheme for pca
d<-expand.grid(lon=1:6,lat=1:6)
d$deme<-seq(1:36)-1


r1 <- with(d,(lon-min(lon))/(max(lon)-min(lon)))
r2 <- with(d,(lat-min(lat))/(max(lat)-min(lat)))
r3 <- 1-sqrt((r1^2+r2^2)/2)
cc <- rbind(r1, r2, r3, 0.5)
rownames(cc) <- c("red", "green", "blue", "alpha")
cols <- rgb(t(cc))
d$cols=cols

pop<-merge(pop,d,by="deme")

cm = merge(cm, pop, by="IID")
re = merge(re, pop, by="IID")

plt_cm<-ggplot(cm,aes(PC1,PC2,color=as.factor(deme)))+
  geom_point(alpha=0.5, size=0.1)+
  theme_bw()+
  scale_color_manual(values=cols)+
  theme(legend.position="none",
        axis.text=element_text(size=10),
        panel.grid = element_blank(),
        plot.margin = unit(rep(0.5,4), "pt"),
        plot.background = element_blank())

plt_re<-ggplot(re,aes(PC1,PC2,color=as.factor(deme)))+
  geom_point(alpha=0.5, size=0.1)+
  theme_bw()+
  scale_color_manual(values=cols)+
  theme(legend.position="none",
        axis.text=element_text(size=10),
        panel.grid = element_blank(),
        plot.margin = unit(rep(0.5,4), "pt"),
        plot.background = element_blank())

plt_combined = plt_cm / plt_re

ggsave(here("plots//pca_plots/plt_complex_pca.pdf"),
       plt_combined,
       height=100,
       width=70,
       units="mm")
