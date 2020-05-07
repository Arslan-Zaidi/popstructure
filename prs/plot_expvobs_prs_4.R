#!/usr/bin/env Rscript


args=commandArgs(TRUE)

library(data.table)
library(ggplot2)
library(dplyr)

pgen = args[1]
iter=args[2]
phenotype=args[3]

plt_output=paste("plt_",iter,".pdf",sep="")
causal_prs = fread(paste( pgen,".",phenotype,".c.sscore" , sep=""))
ncausal_prs = fread(paste( pgen,".",phenotype,".nc.sscore" , sep=""))
gvalue = fread(paste(pgen, ".test.gvalue.sscore", sep=""))


colnames(causal_prs) = colnames(ncausal_prs) = c("IID","dosage","PRS")
colnames(gvalue) = c("IID","dosage","GVALUE")

causal_prs = causal_prs%>%select(IID,PRS)
ncausal_prs = ncausal_prs%>%select(IID,PRS)

causal_prs$ascertainment = "causal_all"
ncausal_prs$ascertainment = "lead"
eprs = rbind(causal_prs,ncausal_prs)

eprs = merge(gvalue%>%c("IID","GVALUE"), eprs, by="IID")

plt = ggplot(eprs,aes(GVALUE,PRS))+
  geom_point(alpha=0.5)+
  theme_bw()+
  geom_abline(intercept=0,slope=1,color="red")+
  stat_smooth(method="lm")+
  facet_wrap(~ascertainment,
    scale="free")+
  labs(x="Expected PRS",
       y="Observed PRS",
       title=paste("Iteration",iter,sep="_"))

ggsave(plt_output, plt)
