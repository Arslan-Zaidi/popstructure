
library(data.table)
#library(ggplot2)
library(dplyr)
library(here)
#library(ggrepel)
library(xtable)

fst100 = fread(here('optimize_migration_rate/grid/tau100/fst/genos_grid_t100_l1e7_ss250_mall_chr1.cmpruned.scikit.fst'))
fst9 = fread(here('optimize_migration_rate/grid/tau-9/fst/genos_grid_t-9_l1e7_ss250_mall_chr1.cmpruned.scikit.fst'))
fst.complex = fread(here('optimize_migration_rate/complex_dem/fst/genos_complex_l1e7_ss250_mall.cmpruned.scikit.fst.gz'))

colnames(fst9) = colnames(fst100) = colnames(fst.complex) =  c("m","a","b","c","fst")
fst100$tau = "Recent"
fst9$tau = "Perpetual"
fst.complex$tau = "Complex"

fst.all = rbind(fst100,fst9,fst.complex)

#function to calculate 95%CI by bootstrapping
boot.fst = function(fst.df,nboot=1000){
  fst.mat = matrix(NA,nrow=nboot,ncol=1)
  for(i in 1:nboot){
    fst.df2 = sample_n(fst.df,size=400,replace=T)
    fst.mat[i,1] = with(fst.df2,sum(a)/sum(a+b+c))
  }
  return(data.table(fst.mean=fst.mat[,1]))
}

fst.boot= fst.all%>%
  group_by(m,tau)%>%
  do(boot.fst(.))%>%
  ungroup()

fst.summary = fst.boot%>%
  group_by(m,tau)%>%
  summarize(median = median(fst.mean),
            lower = quantile(fst.mean, probs=0.025),
            upper = quantile(fst.mean, probs=0.975))%>%
  mutate_at(.vars=c("median","lower","upper"),
            function(x){format(x,digits=2,scientific=T)})%>%
  mutate(fst.ci = paste(median," (",lower," - ",upper,")",sep=""))%>%
  ungroup()%>%
  arrange(tau,m)


#read bplace gwas results


#write function to calculate expected pvalues (for qqplot)

#function to read glm results for longitude
read.gwas<-function(model_prefix,frq_file,maf){
  
  mlong = fread(
    here(paste(model_prefix,".longitude.glm.linear.gz",sep = "")),
    header=F,
    sep="\t")
  
  mlat = fread(
    here(paste(model_prefix,".latitude.glm.linear.gz",sep = "")),
    header=F,
    sep="\t")
  
  colnames(mlong) = colnames(mlat) = c("m","chrom","position","ID","REF","ALT","A1","TEST","nsamples","BETA","SE","TSTAT","P")
  mlong = mlong[,c("m","ID","P")]
  mlat = mlat[, c("m","ID","P")]
  mlong$phenotype = "longitude"
  mlat$phenotype = "latitude"
  mpheno = rbind(mlong,mlat)
  
  mpheno = mpheno%>%
    group_by(m,ID,phenotype)%>%
    mutate(duplicate_no = row_number())%>%
    ungroup()

  
  frq = fread(frq_file,header=F)
  colnames(frq) = c("m","chromosome","ID","ALT","REF","DAF","COUNT")
  frq = frq[,c("m","ID","DAF")]
  frq = frq%>%
    group_by(m,ID)%>%
    mutate(duplicate_no = row_number())%>%
    ungroup()%>%
    filter(DAF>=maf & DAF<=1-maf)
    
  
  mpheno = merge(mpheno,frq,by=c("m","ID","duplicate_no"),sort=F)
  
  
  mpheno = mpheno%>%
    group_by(m,phenotype)%>%
    mutate(CHI=qchisq(P,df=1,lower.tail = F))%>%
    summarize(lambda = median(CHI)/qchisq(0.5,df=1,ncp=0))
  
  return(mpheno)
}


gwas100 = read.gwas("optimize_migration_rate/grid/tau100/bplace_gwas/genos_grid_t100_l1e7_ss250_mall_chr1",
                    "optimize_migration_rate/grid/tau100/genotypes/genos_grid_t100_l1e7_ss250_mall_chr1.frq.afreq.gz",0.03)

gwas9 = read.gwas("optimize_migration_rate/grid/tau-9/bplace_gwas/genos_grid_t-9_l1e7_ss250_mall_chr1",
                  "optimize_migration_rate/grid/tau-9/genotypes/genos_grid_t-9_l1e7_ss250_mall_chr1.frq.afreq.gz",0.03)

gwas.complex = read.gwas("optimize_migration_rate/complex_dem/bplace_gwas/genos_complex_l1e7_ss250_mall",
                         "optimize_migration_rate/complex_dem/genotypes/genos_complex_l1e7_ss250_mall.frq.afreq.gz",0.03)


gwas100$tau = "Recent"
gwas9$tau = "Perpetual"
gwas.complex$tau = "Complex"

gwas.all = rbind(gwas100,gwas9,gwas.complex)

gwas.all = reshape2::dcast(gwas.all, m+tau~phenotype,value.var="lambda")

full_table = merge(fst.summary%>%
                     select(m,tau,fst.ci),
                   gwas.all%>%
                     select(m,tau,latitude,longitude),
                   by=c("tau","m"))

full_table$m = as.factor(full_table$m)
full_table$tau = factor(full_table$tau,levels=c("Recent","Perpetual","Complex"))

full_table = full_table%>%
  arrange(tau,m)


print(xtable(full_table,
             digits=4,
             caption="Mean Fst "),
      include.rownames = F,
      type="latex",
      hline.after = c(10,13),
      file=here("optimize_migration_rate/opt_results_05122020.txt"))



#no longer plotting fst distributions. 
# ggplot(fst100.boot)+
#   geom_density(aes(fst.mean,fill=m))+
#   facet_grid(m~.)+
#   geom_vline(data=fst100.mean,
#              aes(xintercept=fst),
#              color="black")+
#   geom_vline(xintercept=0,color="red")+
#   geom_text(data=fst100.mean,
#              aes(x=0.007,
#                  y=1000,
#                  label=paste("Fst :",round(fst,4),sep="")),
#                  hjust=0,
#              color="black")+
#   theme_minimal()+
#   theme(axis.text.y=element_blank(),
#         panel.grid = element_blank(),
#         strip.text = element_blank(),
#         strip.background = element_blank())+
#   labs(x="Fst", y ="Density",fill='Migration\n rate')+
#   scale_fill_gradient(low="orange",high="white")
# 
# ggplot(fst9)+
#   geom_density(aes(fst,fill=m))+
#   facet_grid(m~.)+
#   geom_vline(data=fst9.mean,
#              aes(xintercept=fst),
#              color="black")+
#   geom_vline(xintercept=0,color="red")+
#   geom_text(data=fst9.mean,
#             aes(x=0.0025,
#                 y=1000,
#                 label=paste("Fst :",round(fst,4),sep="")),
#             hjust=0,
#             color="black")+
#   theme_minimal()+
#   theme(axis.text.y=element_blank(),
#         panel.grid = element_blank(),
#         strip.text = element_blank(),
#         strip.background = element_blank()
#   )+
#   labs(x="Fst", y ="Density",fill='Migration\n rate')+
#   scale_fill_gradient(low="orange",high="white")
# 
