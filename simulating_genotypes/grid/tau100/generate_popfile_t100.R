library(data.table)
library(dplyr)
library(ggplot2)
library(rprojroot)

F=is_rstudio_project$make_fix_file()

pop = fread("gwas/grid/genotypes/tau100/ss500/genos_gridt100_l1e7_ss750_m.pop")
colnames(pop)=c("FID","IID","deme","longitude","latitude")

pop$FID=pop$IID=paste("tsk_",seq(0,26999,1),sep="")

train_sno=seq(1,27000,3)
test_sno=seq(2,27000,3)
sib_sno=seq(3,27000,3)

pop.train = pop[train_sno,]
pop.sib = pop[sib_sno,]
pop.test = pop[test_sno,]

fwrite(pop.train,
       F("gwas/grid/genotypes/tau100/ss500/iid_train.txt"),
       col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

fwrite(pop.test,
       F("gwas/grid/genotypes/tau100/ss500/iid_test.txt"),
       col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

fwrite(pop.sib,
       F("gwas/grid/genotypes/tau100/ss500/iid_sib.txt"),
       col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")





