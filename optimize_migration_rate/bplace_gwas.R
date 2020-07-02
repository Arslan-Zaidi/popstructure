
library(data.table)
library(ggplot2)
library(dplyr)
library(here)



#write function to calculate expected pvalues (for qqplot)
fexp=function(x){
  #x is datframe
  x2=x%>%
    arrange(P)
  nsnps=nrow(x2)
  x2=x2%>%
    mutate(EXP_P=seq(1,nsnps)/nsnps,
           CHI=qchisq(P,df=1,lower.tail = F))
  return(x2)
}


#function to read glm results for longitude
read.long<-function(m){
  
  mlong = fread(
    here(paste("optimize_migration_rate/grid/archived/tau100/bplace_gwas/genos_grid_d36_m",m,"_s500_t100.bplace.Longitude.glm.linear",sep = "")),
    header=T,
    sep="\t")
  
  mlat = fread(
    here(paste("optimize_migration_rate/grid/archived/tau100/bplace_gwas/genos_grid_d36_m",m,"_s500_t100.bplace.Latitude.glm.linear",sep = "")),
    header=T,
    sep="\t")
  
  mlong$phenotype = "longitude"
  mlat$phenotype = "latitude"
  mpheno = rbind(mlong,mlat)
  mpheno = mpheno%>%
    group_by(phenotype)%>%
    do(fexp(.))%>%
    mutate(m = m)
  return(mpheno)
}


m1 = read.long(m = "0.0005")
m2 = read.long(m = "0.001")
m3 = read.long(m = "0.0025")
m4 = read.long(m = "0.005")
m5 = read.long(m = "0.0075")
m6 = read.long(m = "0.01")
m7 = read.long(m = "0.02")
m8 = read.long(m = "0.025")
m9 = read.long(m = "0.03")  
m10 = read.long(m = "0.05")

mall = rbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
mall.mean = mall%>%
  group_by(m,phenotype)%>%
  summarize(lambda = median(CHI)/qchisq(0.5,df=1),
            lambda.c = (lambda-1)*(300/18) +1)

fwrite(mall.mean,sep=" & ")

