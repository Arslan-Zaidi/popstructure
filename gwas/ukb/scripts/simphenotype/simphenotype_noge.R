# simulate phenotypes of individuals for the UKB model
# simulate both 'smooth' north-south effect and 'sharp' effect

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript simphenotype_grid.R <pop file from msprime sim> <output_file> <seed>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

popfile=args[1] # pop file from msprime simulation
output=args[2] #path/name of output file

set.seed(args[3])

#read pop file
pop<-fread(popfile)

#colnames(pop)<-c("FID","IID","deme","latitude","longitude")
#number of latitude steps (dimension of grid)
#lats=length(unique(pop$latitude))

#demes
#demes=unique(pop$deme)

# create a 'smooth clinal phenotype' from north-south
# the effect size between the northern-most and southern-most deme is e_sd
# thus, the effect at each step of latitude is e_sd/lats
pop$smooth = sapply(pop$latitude,
                      function(x){
                        rnorm(n = 1,
                              mean = (x + 1)/3,
                              sd = 1)})

# create a 'smooth clinal phenotype' from west to east
pop$smooth_long = sapply(pop$longitude,
                      function(x){
                        rnorm(n = 1,
                              mean = (x + 1)/3,
                              sd = 1)})

#select one deme randomly to apply a 'sharp' environmental effect
#sharp environmental effect
pop$sharp = sapply(pop$deme,
                 function(x){
                   if(x == "UKJ1"){
                     rnorm(n = 1,
                           mean = 2,
                           sd = 1) }else{
                             rnorm(n = 1,
                                   mean = 0,
                                   sd = 1)
                           }})

#phenotype with no structure
pop$random=rnorm(nrow(pop), mean = 0, sd = 1)

#remove deme and longitude/latitude columns
pop=pop[,c("FID","IID","smooth","smooth_long","sharp","random")]


fwrite(pop,
       output,
       sep="\t",
       row.names=F,
       col.names=T,
       quote=F)
