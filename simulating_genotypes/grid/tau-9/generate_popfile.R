
library(data.table)
library(dplyr)
library(rprojroot)

F=is_rstudio_project$make_fix_file()

pop = fread(F("gwas/grid/genotypes/tau-9/ss500/genos_grid_d36_m0.07_s500_t9.test.pop"))

train = data.table(FID=paste("tsk_",seq(0,17998,2),sep=""),
                   IID=paste("tsk_",seq(0,17998,2),sep=""),
                   deme=rep(seq(0,35),each=250),
                   longitude=pop$latitude,
                   latitude=pop$longitude)

test = data.table(FID=paste("tsk_",seq(1,17999,2),sep=""),
                   IID=paste("tsk_",seq(1,17999,2),sep=""),
                   deme=rep(seq(0,35),each=250),
                   longitude=pop$latitude,
                   latitude=pop$longitude)

fwrite(test,
       F("gwas/grid/genotypes/tau-9/ss500/genos_grid_d36_m0.07_s500_t9.test.pop"),
       sep="\t",
       col.names=TRUE,
       row.names=FALSE,
       quote=FALSE)

fwrite(train,
       F("gwas/grid/genotypes/tau-9/ss500/genos_grid_d36_m0.07_s500_t9.train.pop"),
       sep="\t",
       col.names=TRUE,
       row.names=FALSE,
       quote=FALSE)
