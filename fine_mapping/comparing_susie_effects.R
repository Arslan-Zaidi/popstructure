
gwas_file = "gwas/grid/genotypes/tau100/ss500/train/gwas_results/fixed_effects/ge/gwas_gridt100_train.ge.1.pcs0.smooth.glm.linear.gz"
ldpath = "gwas/grid/genotypes/tau100/ss500/train/ldmats/"

F = is_rstudio_project$make_fix_file()
betas.pcs0 = fread(F("gwas/grid/genotypes/tau100/ss500/train/susie_betas/susie_all.pcs0.betas"))
colnames(betas.pcs0) = c("rep","SNP","A1","BETA","pip","in_cs","marginal","window_name")
betas.pcs0 = betas.pcs0[SNP!="none"]

betas.re = fread(F("gwas/grid/genotypes/tau100/ss500/train/susie_betas/susie_all.re.betas"))
colnames(betas.re) = c("rep","SNP","A1","BETA","pip","in_cs","marginal","window_name")
betas.re = betas.re[SNP!="none"]

ggplot(betas.re,aes(BETA,marginal,color=pip))+
  geom_point()+
  scale_color_gradient(low="yellow",high="blue")+
  geom_abline(intercept=0,slope=1,color="red")+
  theme_classic()
