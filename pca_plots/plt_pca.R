library(data.table)
library(dplyr)
library(ggplot2)
library(here)
library(patchwork)

com_9<-fread(
  here("gwas/grid/genotypes/tau-9/ss500/train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.subsampled.cmpruned.pca.eigenvec"),
  header=T)

rare_9<-fread(
  here("gwas/grid/genotypes/tau-9/ss500/train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.subsampled.repruned.pca.eigenvec"),
  header=T)

com100<-fread(
  here("gwas/grid/genotypes/tau100/ss500/train/genotypes/genos_grid_d36_m0.05_s500_t100.rmdup.train.cmpruned.pca.eigenvec"),
  header=T)

rare100<-fread(
  here("gwas/grid/genotypes/tau100/ss500/train/genotypes/genos_grid_d36_m0.05_s500_t100.rmdup.train.repruned.pca.eigenvec"),
  header=T)

#multiply PC2 of rare100 by -1 to flip axis to match orientation of demes in other plots
#rare100$PC2=rare100$PC2*-1

colnames(com_9)<-colnames(rare_9)<-colnames(com100)<-colnames(rare100)<-c("FID","IID",paste("PC",seq(1,100),sep=""))

pop=fread("gwas/grid/genotypes/tau100/ss500/genos_grid_d36_m0.05_s500_t100.train.pop")
colnames(pop)=c("FID","IID","deme","longitude","latitude")

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

com_9<-merge(com_9,pop,by=c("FID","IID"))
com100<-merge(com100,pop,by=c("FID","IID"))
rare_9<-merge(rare_9,pop,by=c("FID","IID"))
rare100<-merge(rare100,pop,by=c("FID","IID"))



plt_common_9<-ggplot(com_9,aes(PC1,PC2,color=as.factor(deme)))+
  geom_point(alpha=0.5, size=0.1)+
  theme_bw()+
  scale_color_manual(values=cols)+
  theme(legend.position="none",
        axis.text=element_text(size=10),
        panel.grid = element_blank(),
        plot.margin = unit(rep(0.5,4), "pt"),
        plot.background = element_blank())

plt_rare_9<-ggplot(rare_9,aes(PC1,PC2,color=as.factor(deme)))+
  geom_point(alpha=0.5, size=0.1)+
  theme_bw()+
  scale_color_manual(values=cols)+
  theme(legend.position="none",
        axis.text = element_text(size=10),
        panel.grid=element_blank(),
        plot.margin = unit(rep(0.5,4), "pt"))

plt_common100<-ggplot(com100,aes(PC1,PC2,color=as.factor(deme)))+
  geom_point(alpha=0.5, size=0.1)+
  theme_bw()+
  scale_color_manual(values=cols)+
  theme(legend.position="none",
        axis.text=element_text(size=10),
        panel.grid=element_blank(),
        plot.margin = unit(rep(0.5,4), "pt"),
        plot.background=element_blank())

plt_rare100<-ggplot(rare100,aes(PC1,PC2,color=as.factor(deme)))+
  geom_point(alpha=0.5, size=0.1)+
  theme_bw()+
  scale_color_manual(values=cols)+
  theme(legend.position="none",
        axis.text=element_text(size=10),
        panel.grid=element_blank(),
        plot.margin = unit(rep(0.5,4), "pt"),
        plot.background = element_blank())

plt_combined <- (plt_common_9 + plt_common100) / (plt_rare_9 + plt_rare100)

ggsave(here("analyses/pca_plots/plt_grid_pca.pdf"),
            plt_combined,
            height=100,
            width=120,
            units="mm")


################# IBD sharing plot ################

ibd=fread(
  here("gwas/grid/genotypes/tau100/ss250/germline/germ.cm_ibd.eigenvec"))

ibd2=ibd[,c("msp","iid"):=tstrsplit(V2,split="_")]
ibd2[,iid:=paste("msp",as.numeric(iid)-1,sep="_")]

ibd2=merge(pop,ibd2,by.x="IID",by.y="iid",all.x=T,sort=F)

plt_ibd=ggplot(ibd2,aes(V4,V5,color=as.factor(deme)))+
  geom_point(show.legend = F,
             alpha=0.5, size=0.1)+
  scale_color_manual(values=cols)+
  theme_bw()+
  theme(legend.position="none",
        axis.text=element_text(size=10),
        plot.title = element_text(size = 11),
        axis.title = element_text(size = 11),
        panel.grid=element_blank())+
  labs(x="PC2",y="PC3",
       title="PCA on IBD-sharing")

ggsave(here("analyses/pca_plots/plt_grid_ibd_pca.pdf"),
       plt_ibd,
       height=70,
       width=70,
       units="mm")


       
       
       
