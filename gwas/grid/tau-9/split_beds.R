
library(data.table)
library(rprojroot)

F = is_rstudio_project$make_fix_file()

dat = fread(F("gwas/grid/genotypes/tau-9/ss500/train/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.frq.afreq"))

colnames(dat)[1]="CHROM"
dat = dat[order(CHROM),]

bims=seq(1,nrow(dat),15e4)

for(i in 1:length(bims)){
  fwrite(dat[bims[i]:(bims[i+1]-1),.(ID)],
         F(paste("gwas/grid/genotypes/tau-9/ss500/train/genotypes/bim_",i,
                 ".snplist",sep="")),
         col.names=FALSE,row.names=FALSE,quote=FALSE
         )
}
