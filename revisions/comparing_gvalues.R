
gvalues2 = fread(F("gwas/grid/genotypes/tau100/ss500/revisions/prs_prediction/gvalues/gvalue.p2.all.sscore"))
gvalues3 = fread(F("gwas/grid/genotypes/tau100/ss500/revisions/prs_prediction/gvalues/gvalue.p3.all.sscore"))

colnames(gvalues2)=colnames(gvalues3)=c("rep","IID","gvalue")

prs.a1_r2_p3 = fread(F("gwas/grid/genotypes/tau100/ss500/revisions/prs_prediction/prs/a1_r2_p3.smooth.pcs0.all.sscore"))
prs.a1_r3s_p2 = fread(F("gwas/grid/genotypes/tau100/ss500/revisions/prs_prediction/prs/a1_r3s_p2.smooth.ncorr.all.sscore"))

colnames(prs.a1_r2_p3) = colnames(prs.a1_r3s_p2) = c("rep","IID","pcs0","cm","re","cmre")

prs.a1_r3s_p2 = merge(gvalues2,prs.a1_r3s_p2,by=c("rep","IID"))
prs.a1_r2_p3 = merge(gvalues3,prs.a1_r2_p3,by=c("rep","IID"))


prs.a1_p3 = fread(F("gwas/grid/genotypes/tau100/ss500/revisions/prs_prediction/prs1sample/prs.smooth.a1_p3.all.nc.sscore"),
                  header=FALSE)

colnames(prs.a1_p3) = c("rep","IID","pcs0","cm","re","cmre")


prs.a1_r2_p3 = merge(prs.a1_r2_p3, prs.a1_p3, by=c("rep","IID"))


ggplot(prs.a1_r2_p3, aes(pcs0.x, pcs0.y))+
  geom_point()



gvalues2.r2 = gvalues2%>%
  group_by(rep)%>%
  summarize(r2 = cor(gvalue,pcs0)^2)

gvalues3.r2 = gvalues3%>%
  group_by(rep)%>%
  summarize(r2 = cor(gvalue,pcs0)^2)

pop2 = fread(F("gwas/grid/genotypes/tau100/ss500/iid_test.txt"))
pop3 = fread(F("gwas/grid/genotypes/tau100/ss500/iid_sib.txt"))

pop2 = merge(gvalues2,pop2,by="IID")
pop3 = merge(gvalues3,pop3,by="IID")

pop2.mean = pop2%>%
  group_by(rep)%>%
  mutate(prs.centered = pcs0 - mean(pcs0))%>%ungroup()%>%
  group_by(rep,deme,longitude,latitude)%>%
  summarize(prs.c.mean = mean(prs.centered - gvalue))

pop3.mean = pop3%>%
  group_by(rep)%>%
  mutate(prs.centered = pcs0 - mean(pcs0))%>%ungroup()%>%
  group_by(rep,deme,longitude,latitude)%>%
  summarize(prs.c.mean = mean(prs.centered - gvalue))

ggplot(pop3.mean,aes(longitude,latitude,fill=prs.c.mean))+
  geom_tile()+
  

pop2.rlat = pop2%>%
  group_by(rep)%>%
  summarize(rlat = cor(pcs0,latitude))

pop3.rlat = pop3%>%
  group_by(rep)%>%
  summarize(rlat = cor(pcs0,latitude))

  
ggplot(pop2,aes(pcs0,latitude))+
  geom_point()

