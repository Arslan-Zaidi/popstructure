

args = commandArgs(TRUE)
pheno = args[1]
plot_title = args[2]

if( length(args) != 2 ){stop("Usage: Rscript plot_prs_all.R <phenotype> <plot title>")}

library(ggplot2)
library(dplyr)
library(data.table)
library(rprojroot)

#specify root of the directory
F = is_rstudio_project$make_fix_file()

#load pop file which contains latitude,longitude information
pop.test=fread(F("gwas/grid/genotypes/tau100/ss500/genos_grid_d36_m0.05_s500_t100.test.pop"))
pop.test$FID=pop.test$IID=paste("tsk_",seq(1,17999,2),sep="")

#load PRS files
#1: prs constructued using all causal variants
#2. prs using causal variants (pvalue< 5e-04)
#3. prs using lead SNPs
#prs1=fread(F(paste("gwas/grid/genotypes/tau100/ss500/test/prs/gridt100_prs_",pheno,".all.c.sscore.gz",sep="")))
prs2=fread(F(paste("gwas/grid/genotypes/tau100/ss500/test/prs/gridt100_prs_",pheno,".all.c.p.sscore.gz",sep="")))
prs3=fread(F(paste("gwas/grid/genotypes/tau100/ss500/test/prs/gridt100_prs_",pheno,".all.nc.sscore.gz",sep="")))

colnames(prs2)=colnames(prs3)=c("rep","IID","dosage_sum","pcs0","cm","re","cmre")
#prs1$ascertainment = "all_causal"
prs2$ascertainment = "causal_p"
prs3$ascertainment = "lead_snp"

#rbind so all analyses can be performed simultaneously
prs=rbind(prs2,prs3)

#removing cmre for now
#prs = prs%>%
#    select(-cmre)

prs=merge(prs,pop.test,by="IID")

gvalue = fread(F("gwas/grid/genotypes/tau100/ss500/test/gvalue/genos_grid_d36_m0.05_s500_t100.rmdup.test.all.gvalue.sscore.gz"))
colnames(gvalue) = c("rep","IID","dosage","gvalue")

prs = merge(prs, gvalue, by=c("rep","IID"))

#melt to long format
mprs=melt(prs%>%
            select(-c(dosage_sum,dosage,FID)),
          id.vars=c("rep","IID","gvalue","ascertainment","deme","longitude","latitude"),
          variable.name="correction",
          value.name="prs")

#remove variation due to simulated genetic value
mprs.adj = mprs%>%
  group_by(rep,correction,ascertainment)%>%
  mutate(prs.adjusted = resid(lm(prs~gvalue)),
         rlat = cor(prs.adjusted, latitude),
         rlong = cor(prs.adjusted, longitude))%>%
  ungroup()

#calculate mean prs adjusted for each deme
mprs.sum = mprs.adj%>%
  group_by(correction,ascertainment,longitude,latitude)%>%
  summarize(mean.prs = mean(prs.adjusted))%>%
  ungroup()

#calculate mean of rlat and rlong across reps
mprs.r = mprs.adj %>%
  group_by(correction,ascertainment)%>%
  summarize(rlat = mean(rlat),
            rlong = mean(rlong))%>%
  ungroup()


mprs.r = mprs.r %>%
  mutate(xlat = 0,
         ylat = 5,
         label.lat = paste("r-lat: ", round(rlat,3), sep=""),
         xlong = 0,
         ylong = 3,
         label.long = paste("r-long: ", round(rlong,3),sep=""))


labels_prs=c(
  all_causal="All causal",
  causal_p="Causal\n(P<5e-04)",
  lead_snp="Lead SNP",
  pcs0="No correction",
  cm="Common PCA",
  re="Rare PCA",
  cmre="Common+rare PCA"
)

if(pheno %in% c("smooth","smooth_long","grandom")){

  plt_prs_phe=ggplot() +
    geom_tile(data = mprs.sum,
              aes(longitude, latitude, fill = mean.prs),
              show.legend = FALSE) +
    geom_text(data = mprs.r,
              aes(xlat, ylat, label = label.lat),
              hjust = 0,
              vjust=1,
              size = 9/(14/5)) +
    geom_text(data = mprs.r,
              aes(xlong, ylong, label = label.long),
              hjust = 0,
              vjust=1,
              size = 9/(14/5))+
    theme_bw()+
    facet_grid(correction ~ ascertainment,
               labeller=as_labeller(labels_prs)) +
    scale_fill_gradient2(high = "#fc8d59", mid = "#ffffbf", low = "#91bfdb")+
    labs(x="Longitude", y="Latitude",title=plot_title)+
    theme(plot.title=element_text(hjust=0,size=11),
          strip.text = element_text(size=10),
          panel.grid = element_blank())

}

if(pheno %in% c("sharp")){

  plt_prs_phe=ggplot() +
    geom_tile(data = mprs.sum,
              aes(longitude, latitude, fill = mean.prs),
              show.legend = FALSE) +
    annotate(geom="text",
             x=2, y=0, label = "*", vjust = 0.7) +
    theme_bw()+
    geom_text(data = mprs.r,
              aes(xlat, ylat, label = label.lat),
              hjust = 0,
              vjust=1,
              size = 9/(14/5)) +
    geom_text(data = mprs.r,
              aes(xlong, ylong, label = label.long),
              hjust = 0,
              vjust=1,
              size = 9/(14/5))+
    facet_grid(correction ~ ascertainment,
               labeller=as_labeller(labels_prs)) +
    scale_fill_gradient2(high = "#fc8d59", mid = "#ffffbf", low = "#91bfdb")+
    labs(x="Longitude", y="Latitude",title=plot_title)+
    theme(plot.title=element_text(hjust=0,size=11),
          strip.text = element_text(size=10),
          panel.grid = element_blank())
}


ggsave(F(paste("plots/prs/grid/tau100/plt_prs_",pheno,".pdf",sep="")),
       plt_prs_phe,
       height=95,
       width=95,
       units="mm")
