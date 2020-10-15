
#library(ggplot2)
library(data.table)
library(dplyr)
library(rprojroot)
library(susieR)
library(matrixcalc)

F = is_rstudio_project$make_fix_file()


fclump=function(gwas_df,pcutoff=5e-04){
  if(missing(pcutoff)){
    df.red=gwas_df
  }else{
    df.red = gwas_df[P<pcutoff]
  }
  
  
  for(i in 1:101){
    if(i == 1){
      start=0
      stop=start+1e5
    }else{
      
      start=((i-1)*1e5) +1
      stop=start + 1e5 -1
    }
    
    df.red[ (POS>=start &POS<stop) , window:=i]
  }
  df.red[,window_name:= paste(CHROM,window,sep="_") ]
  return(df.red)
  
}

flead=function(df){
  df.red = df[P == min(P)]
  return(df.red)
}

gwas2 = fclump(gwas,5e-04)
gwas.red2 = gwas2[,flead(.SD),by=.(window_name)]
gwas.red2 = unique(gwas.red2, by ="window_name")

ggplot(gwas.red,
       aes(POS,-log10(P)))+
  geom_point()+
  geom_hline(yintercept=-log10(5e-04))

ggplot(gwas.red%>%
         filter(P<5e-04),
       aes(POS,-log10(P)))+
  geom_point()+
  geom_hline(yintercept=-log10(5e-04))+
  xlim(c(0,1e7))+
geom_point(data=gwas.red2%>%
             filter(P<5e-06),
           aes(POS,-log10(P)),
           color="red")+
  geom_segment(x=0,xend=1e5,
               y=40,yend=40)+
  theme_classic()

freq = fread(F("gwas/grid/genotypes/tau100/ss500/train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.snps.frq.afreq"))
gwas.red2 = merge(gwas.red2,freq,by="ID")


coordinates = matrix(NA, 100,1)
for(i in 1:100){
    start=((i-1)*1e5) +1
    #stop=start+1e5
  coordinates[i,1] = start
}

all.coordinates = expand.grid(c(1:20),
                              coordinates)

all.coordinates = cbind(all.coordinates,
                        all.coordinates[,2]+(1e5-1))

colnames(all.coordinates) = c("chrom","start","stop")
all.coordinates = as.data.table(all.coordinates)
all.coordinates = all.coordinates%>%
  arrange(chrom,start)

fwrite(all.coordinates,
       F("gwas/grid/genotypes/tau100/ss500/train/ldmats/window_coordinates.txt"),
       col.names = FALSE,
       row.names = FALSE,
       quote=FALSE,
       sep="\t")





