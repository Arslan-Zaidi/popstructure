
prs.pcs0 = fread(F("gwas/grid/genotypes/tau100/ss500/test/susie_prs/prs_susie_smooth.all.pcs0.sscore"))
prs.re = fread(F("gwas/grid/genotypes/tau100/ss500/test/susie_prs/prs_susie_smooth.all.re.sscore"))

colnames(prs.pcs0) = colnames(prs.re)=c("rep","IID","lead","marginal")

gvalues = fread(F("gwas/grid/genotypes/tau100/ss500/test/gvalue/genos_grid_d36_m0.05_s500_t100.rmdup.test.all.gvalue.sscore.gz"))
colnames(gvalues) = c("rep","IID","dosage","gvalue")
#gvalues = gvalues%>%filter(rep==1)

pop = fread(F("gwas/grid/genotypes/tau100/ss500/iid_test.txt"))

prs.pcs0$correction="pcs0"
prs.re$correction="re"
prs.both= rbind(prs.pcs0,prs.re)

prs.both = merge(prs.both,pop)
prs.both = merge(prs.both,gvalues,by=c("IID","rep"))

mprs.both = reshape2::melt(prs.both,
                  id.vars=c("IID","rep","correction","FID","deme","longitude","latitude","dosage","gvalue"),
                  value.name="prs",
                  variable.name="estimate")

mprs.summary = mprs.both%>%
  group_by(correction,rep,estimate)%>%
  summarize(rlat = cor(prs,latitude),
            r2 = cor(prs,gvalue)^2)


plt.bias = ggplot(mprs.summary)+
  geom_histogram(aes(rlat),bins=10)+
  facet_grid(correction~estimate)+
  theme_classic()

plt.r2 = ggplot(mprs.summary)+
  geom_histogram(aes(r2),bins=10)+
  facet_grid(correction~estimate)+
  theme_classic()

