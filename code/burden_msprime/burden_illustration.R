
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

burden=fread("burden_msprime/burden_msprime_t100_rho0_clustering_1gene.txt")

colnames(burden)=paste("g",1:8,sep="")
burden = burden%>%
  replace(is.na(.),0)
pop = fread("gwas/grid/genotypes/tau100/ss500/genos_grid_d36_m0.05_s500_t100.train.pop")

burden=cbind(pop,burden)

mburden=burden%>%
  melt(id.vars=c("FID","IID","deme","longitude","latitude"),
       variable.name="length",
       value.name="burden")

mburden.sum = mburden%>%
  group_by(length,longitude,latitude)%>%
  summarize(sum.burden=sum(burden))

plt_rho0=ggplot(mburden.sum,aes(longitude,latitude,fill=sum.burden))+
  geom_tile()+
  facet_wrap(length~.)+
  scale_fill_viridis_c()+
  theme_void()

burden=fread("burden_msprime/burden_msprime_t100_rho1e-08_clustering_1gene.txt")
colnames(burden)=paste("g",1:8,sep="")
burden = burden%>%
  replace(is.na(.),0)
burden=cbind(pop,burden)

mburden=burden%>%
  melt(id.vars=c("FID","IID","deme","longitude","latitude"),
       variable.name="length",
       value.name="burden")

mburden.sum = mburden%>%
  group_by(length,longitude,latitude)%>%
  summarize(sum.burden=sum(burden))

plt_rho=ggplot(mburden.sum,aes(longitude,latitude,fill=sum.burden))+
  geom_tile()+
  facet_wrap(length~.)+
  scale_fill_gradient(low="white",high="red")+
  theme_bw()


ggsave("plots/burden_msprime/plt_burden_illustration_spatial_rho0.pdf",plt_rho0)
ggsave("plots/burden_msprime/plt_burden_illustration_spatial_rho.pdf",plt_rho)
