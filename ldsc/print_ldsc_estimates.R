

library(data.table)
library(dplyr)
library(ggplot2)

cm = fread("gwas/complex_dem/train/gwas_complex_train.ge.all.cm.grsmshp.ldsc.pt")
re = fread("gwas/complex_dem/train/gwas_complex_train.ge.all.re.grsmshp.ldsc.pt")
cmre = fread("gwas/complex_dem/train/gwas_complex_train.ge.all.cmre.grsmshp.ldsc.pt")

colnames(cm) = colnames(re) = colnames(cmre)= c("phenotype","iteration","variable","value")

cm$correction = "cm"
re$correction = "re"
cmre$correction = "cmre"

ldsc = rbind(cm,re,cmre)
ldsc.sum = ldsc%>%
  filter(variable=="Intercept")%>%
  group_by(phenotype,correction)%>%
  summarize(lambda = mean(value))

