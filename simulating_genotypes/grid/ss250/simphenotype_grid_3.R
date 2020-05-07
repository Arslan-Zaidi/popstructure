# simulate phenotypes of individuals for the grid model
#simulate both 'smooth' north-south effect and 'sharp' effect

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript simphenotype_grid.R <pop file from msprime sim> <e_sd environmental effect size - of standard deviation of the phenotype> <sample size per deme>")}

#load libraries
library(dplyr)
library(data.table)
library(here)

popfile=args[1] # pop file from msprime simulation
e_sd=as.numeric(args[2]) # environmental effect size - of standard deviation of the phenotype
n=as.numeric(args[3])

#read pop file
pop<-fread(here(paste("gwas/grid/",popfile,sep="")))

colnames(pop)<-c("FID","IID","deme","latitude","longitude")
#number of latitude steps (dimension of grid)
lats=length(unique(pop$latitude))

#demes
demes=unique(pop$deme)

# create a 'smooth clinal phenotype' from north-south
# the effect size between the northern-most and southern-most deme is e_sd
# thus, the effect at each step of latitude is e_sd/lats
pop$smooth=sapply(pop$latitude,function(x){
  rnorm( 1, mean = x*e_sd/lats, sd=1)
})

#select one deme randomly to apply a 'sharp' environmental effect
random.deme=sample(demes,1)
pop$sharp=sapply(pop$deme,function(x){
  if(x==random.deme){
    y=rnorm(1,mean=e_sd,sd=1)}
  else{
    y=rnorm(1,mean=0,sd=1)}
  return(y)
})

#write to phenotype file for GWAS
filename=paste("grid_d36_s",n,"_e",e_sd,".pheno",sep="")

fwrite(pop,
       here(paste("gwas/grid/phenotypes/",filename,sep="")),
       sep="\t",
       row.names=F,
       col.names=T,
       quote=F)
