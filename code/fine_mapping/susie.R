library(susieR)
library(data.table)
library(rprojroot)
library(dplyr)
library(Matrix)
library(tidyr)

F = is_rstudio_project$make_fix_file()

R = as.matrix(fread(F("gwas/grid/genotypes/tau100/ss500/train/ldmats/ld_chr1_5000_5100.ld.gz")))

npd = nearPD(R)
R2 = as.matrix(npd$mat)


gwas = fread(F("gwas/grid/genotypes/tau100/ss500/train/gwas_results/fixed_effects/ge/gwas_gridt100_train.ge.1.pcs0.smooth.glm.linear.gz"))
colnames(gwas)[1]="CHROM"
gwas = gwas[CHROM==1]
gwas = gwas[POS>=5e6 & POS<5.1e6]

causal = fread(F("gwas/grid/genotypes/tau100/ss500/train/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.all.thinned_100kb.effects"))
colnames(causal) = c("rep","SNP","allele","esize")

causal = causal%>%
  filter(rep==1)%>%
  separate(SNP, sep="_",
           into=c("chrom","position","ref","alt"),
           remove = FALSE)

causal$chrom = as.numeric(causal$chrom)
causal$position = as.numeric(causal$position)

causal.red = causal%>%
  filter(chrom==1 & position>=5e6 & position<5.1e6)

gwas = gwas%>%
  mutate(b = case_when(ID%in%causal.red$SNP ~ 1,
                       TRUE ~ 0))

fitted_rss = susie_rss(gwas$T_STAT, 
                       R2, L=10, check_z = FALSE)

summary(fitted_rss)$cs

susie_plot(fitted_rss, y="PIP", b=gwas$b)

ggplot(gwas,aes(POS,-log10(P)))+
  geom_point()+
  theme_minimal()+
  geom_point(data=gwas%>%filter(b==1),
             aes(POS,-log10(P)),
             color="red")

#now with rare PC correction
gwas.re = fread(F("gwas/grid/genotypes/tau100/ss500/train/gwas_results/fixed_effects/ge/gwas_gridt100_train.ge.1.cmre.smooth.glm.linear.gz"))
colnames(gwas.re)[1]="CHROM"
gwas.re = gwas.re[CHROM==1]
gwas.re = gwas.re[POS<=1e5]

gwas.re = gwas.re%>%
  mutate(b = case_when(ID%in%causal.red$SNP ~ 1,
                       TRUE ~ 0))

fitted_rss.re = susie_rss(gwas.re$T_STAT, 
                       R2, L=1, check_z = FALSE)

summary(fitted_rss.re)$cs

susie_plot(fitted_rss.re, y="PIP", b=gwas$b)

############## let's try with common

R_cm = as.matrix(fread(F("gwas/grid/genotypes/tau100/ss500/train/genotypes/ldmats/ldmat_chr1_1_100kb_common.ld.gz")))
npd_cm = nearPD(R_cm)
R2_cm = as.matrix(npd_cm$mat)

ids_cm = fread(F("gwas/grid/genotypes/tau100/ss500/train/genotypes/ldmats/genos_chr1_1_100kb_common.bim"))
colnames(ids_cm) = c("rep","SNP","cm","position","ref","alt")

gwas_cm = gwas%>%
  filter(ID %in% ids_cm$SNP)

fitted_rss_cm = susie_rss(gwas_cm$T_STAT, 
                          R2_cm, L=1, check_z = FALSE)

summary(fitted_rss_cm)$cs

susie_plot(fitted_rss_cm, y="PIP", b=gwas_cm$b)

################## correction with rare PCs

gwas.re = gwas.re%>%
  filter(ID %in% ids_cm$SNP)

fitted_rss.re = susie_rss(gwas.re$T_STAT, 
                          R2_cm, L=1, check_z = FALSE)

summary(fitted_rss.re)$cs

susie_plot(fitted_rss.re, y="PIP", b=gwas_cm$b)
