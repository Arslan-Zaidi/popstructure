

args = commandArgs(TRUE)
pheno = args[1]
plot_title = args[2]

if( length(args) != 2 ){stop("Usage: Rscript plot_prs_all.R <phenotype> <plot title>")}

library(ggplot2)
library(dplyr)
library(data.table)
library(here)


pop.test=fread(here("gwas/grid/genotypes/tau100/ss500/genos_grid_d36_m0.05_s500_t100.test.pop"))
pop.test$FID=pop.test$IID=paste("tsk_",seq(1,17999,2),sep="")

prs1=fread(here(paste("gwas/nostructure/test/prs/nostructure_prs_",pheno,".all.c.sscore",sep="")))
prs2=fread(here(paste("gwas/nostructure/test/prs/nostructure_prs_",pheno,".all.c.p.sscore",sep="")))
prs3=fread(here(paste("gwas/nostructure/test/prs/nostructure_prs_",pheno,".all.nc.sscore",sep="")))

colnames(prs1)=colnames(prs2)=colnames(prs3)=c("IID","dosage_sum","pcs0","cm","re")
prs1$ascertainment = "all_causal"
prs2$ascertainment = "causal_p"
prs3$ascertainment = "lead_snp"

prs=rbind(prs1,prs2,prs3)
prs=merge(prs,pop.test,by="IID")

mprs=melt(prs,
          id.vars=c("IID","dosage_sum","ascertainment","FID","deme","longitude","latitude"),
          variable.name="correction",
          value.name="prs")

mprs.sum = mprs%>%
  group_by(correction,ascertainment,latitude,longitude)%>%
  summarize(prs=mean(prs))

phe.prs = mprs.sum%>%
  group_by(correction,ascertainment)%>%
  mutate(prs.scaled = prs - mean(prs))

phe.prs.r = phe.prs %>%
  group_by(correction, ascertainment) %>%
  summarize(rlat = cor(prs.scaled, latitude),
            rlong = cor(prs.scaled, longitude),
            plat = summary(lm(prs.scaled ~ latitude))$coefficients[2, 4],
            plong = summary(lm(prs.scaled ~ longitude))$coefficients[2,4])%>%
  ungroup()


# phe.prs.r = phe.prs.r %>%
#   mutate(x = 0,
#          y = 5,
#          label = paste("rlat:", round(rlat,2),
#                     "\nplat:", formatC(plat, format = "e", digits = 2),
#                     "\nrlong:", round(rlong,2),
#                     "\nplong:", formatC(plong, format = "e", digits = 2), sep=""))

phe.prs.r = phe.prs.r %>%
  mutate(xlat = 0,
         ylat = 5,
         label.lat = case_when(plat < 0.01/9 ~ paste("r-lat: ", round(rlat,2),"**", sep=""),
                           plat < 0.01 ~ paste("r-lat:", round(rlat,2),"*", sep=""),
                           TRUE ~ paste("r-lat:", round(rlat,2), sep="")),
         xlong = 0,
         ylong = 3,
         label.long = case_when(plong < 0.01/9 ~ paste("r-long: ", round(rlong,2),"**", sep=""),
                                plong < 0.01 ~ paste("r-long:", round(rlong,2),"*", sep=""),
                                TRUE ~ paste("r-long:", round(rlong,2), sep="")))


labels_prs=c(
  all_causal="All causal",
  causal_p="Causal\n(P<5e-04)",
  lead_snp="Lead SNP",
  pcs0="No correction",
  cm="Common PCA",
  re="Rare PCA"
)

if(pheno %in% c("smooth","smooth_long","grandom")){

  plt_prs_phe=ggplot() +
    geom_tile(data = phe.prs,
              aes(longitude, latitude, fill = prs.scaled),
              show.legend = F) +
    geom_text(data = phe.prs.r,
              aes(xlat, ylat, label = label.lat),
              hjust = 0,
              vjust=1,
              size = 9/(14/5)) +
    geom_text(data = phe.prs.r,
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
    geom_tile(data = phe.prs,
              aes(longitude, latitude, fill = prs.scaled),
              show.legend = F) +
    annotate(geom="text",
             x=2, y=0, label = "*", vjust = 0.7) +
    theme_bw()+
    geom_text(data = phe.prs.r,
              aes(xlat, ylat, label = label.lat),
              hjust = 0,
              vjust=1,
              size = 9/(14/5)) +
    geom_text(data = phe.prs.r,
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

ggsave(here(paste("plots/prs/nostructure/plt_prs_",pheno,".pdf",sep="")),
       plt_prs_phe,
       height=95,
       width=95,
       units="mm")
