

args = commandArgs(TRUE)
pheno = args[1]
plot_title = args[2]

if( length(args) != 2 ){stop("Usage: Rscript plot_prs_all.R <phenotype> <plot_title>")}

library(ggplot2)
library(dplyr)
library(data.table)
library(here)


pop.test=fread(here("gwas/complex_dem/genos_complex_l1e7_ss500_m0.07.test.pop"))
#pop.test$FID=pop.test$IID=paste("tsk_",seq(1,17999,2),sep="")

#prs1=fread(here(paste("gwas/complex_dem/test/prs/complexdem_prs_",pheno,".all.c.sscore.gz",sep="")))
prs2=fread(here(paste("gwas/complex_dem/test/prs/complexdem_prs_",pheno,".all.c.p.sscore.gz",sep="")))
prs3=fread(here(paste("gwas/complex_dem/test/prs/complexdem_prs_",pheno,".all.nc.sscore.gz",sep="")))

colnames(prs2)=colnames(prs3)=c("rep","IID","dosage_sum","pcs0","cm","re","cmre")
#prs1$ascertainment = "all_causal"
prs2$ascertainment = "causal_p"
prs3$ascertainment = "lead_snp"

prs_df=rbind(prs2,prs3)

prs_df=merge(prs_df,pop.test,by="IID")

gvalue_df = fread(here("gwas/complex_dem/test/gvalue/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.test.all.gvalue.sscore.gz"))
colnames(gvalue_df) = c("rep","IID","dosage","gvalue")
gvalue_df = gvalue_df[,c('rep','IID','gvalue')]

prs_df = merge(prs_df, gvalue_df, by=c("rep","IID"))

#melt to long format
mprs_df=melt(prs_df%>%
               select(-c(dosage_sum,FID)),
             id.vars=c("rep","IID","gvalue","ascertainment","deme","longitude","latitude"),
             variable.name="correction",
             value.name="prs")

#remove variation due to simulated genetic value
mprs.adj = mprs_df%>%
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


labels_prs=c(
  causal_p="Causal",
  lead_snp="Lead SNP",
  pcs0="No correction",
  cm="Common PCA",
  re="Rare PCA",
  cmre="Common + rare"
)

if(pheno %in% c("smooth","smooth_long","grandom")){
  
  plt_prs_phe=ggplot() +
    geom_tile(data = mprs.sum,
              aes(longitude, latitude, fill = mean.prs),
              show.legend = T) +
    theme_bw()+
    facet_grid(correction ~ ascertainment,
               labeller=as_labeller(labels_prs)) +
    scale_fill_gradient2(high = "#fc8d59", mid = "#ffffbf", low = "#91bfdb")+
    labs(x="Longitude", y="Latitude", title=plot_title, fill="Mean\nPRS")+
    theme(plot.title=element_text(hjust=0,size=11),
          strip.text = element_text(size=10),
          panel.grid = element_blank())
  
}

if(pheno %in% c("sharp")){
  
  plt_prs_phe=ggplot() +
    geom_tile(data = mprs.sum,
              aes(longitude, latitude, fill = mean.prs),
              show.legend = T) +
    annotate(geom="text",
             x=2, y=0, label = "*", vjust = 0.7) +
    theme_bw()+
    facet_grid(correction ~ ascertainment,
               labeller=as_labeller(labels_prs)) +
    scale_fill_gradient2(high = "#fc8d59", mid = "#ffffbf", low = "#91bfdb")+
    labs(x="Longitude", y="Latitude", title=plot_title, fill="Mean\nPRS")+
    theme(plot.title=element_text(hjust=0,size=11),
          strip.text = element_text(size=10),
          panel.grid = element_blank())
}

ggsave(here(paste("plots/prs/complex_dem/plt_prs_",pheno,".pdf",sep="")),
       plt_prs_phe,
       height=110,
       width=95,
       units="mm")
