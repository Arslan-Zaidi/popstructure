
#script takes in genetic values and generates phenotypes with and without stratification

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript simphenotype_ge.R <genetic value file> <pop file> <output_file> <seed>")}

library(data.table)
library(dplyr)

gvalue_file = args[1] #genetic values
popfile = args[2] # pop file from msprime simulation
output_file = args[3] #name of output file
set.seed(args[4])

prs=fread(gvalue_file)
colnames(prs)<-c("IID","dosage","prs")

sample_size=nrow(prs)


#load file containing the deme id and latitude and longitude for each individual
pop=fread(popfile)

#add this info to prs file
prs=merge(prs, pop, by="IID", sort=F)

##### simulate environmental effects
#no 'environmental' effect
#NOTE: SIGMA_G SET AT 0.8
prs$grandom = rnorm(sample_size,0, sqrt(1 - 0.8))

#smooth effect - latitude (OG effect)
prs$smooth = sapply(prs$latitude,
                      function(x){
                        rnorm(n = 1,
                              mean = (x + 1)/3,
                              sd = sqrt(1 - 0.8))})

#smooth effect  - longitude
prs$smooth_long = sapply(prs$longitude,
                      function(x){
                        rnorm(n = 1,
                              mean = (x + 1)/3,
                              sd = sqrt(1 - 0.8))})


#sharp environmental effect
prs$sharp = sapply(prs$deme,
                 function(x){
                   if(x == 2){
                     rnorm(n = 1,
                           mean = 2,
                           sd = sqrt(1 - 0.8)) }else{
                             rnorm(n = 1,
                                   mean = 0,
                                   sd = sqrt(1 - 0.8))
                           }})

#scale each so that their variances are 0.2 - to adjust heritability to 0.8
prs$grandom = scale( prs$grandom, scale = T) * sqrt( 1 - 0.8)
prs$sharp = scale( prs$sharp , scale = T) * sqrt( 1 - 0.8)
prs$smooth = scale( prs$smooth, scale = T) * sqrt(1 - 0.8)
prs$smooth_long = scale( prs$smooth_long, scale = T) * sqrt(1 - 0.8)

#add prs to each of the environmental effects
prs = prs %>%
  mutate(grandom = prs + grandom,
         smooth = prs + smooth,
         smooth_long = prs + smooth_long,
         sharp = prs + sharp)

fwrite(
  prs%>%
  mutate(FID=IID)%>%
  select(FID,IID,grandom,smooth,smooth_long,sharp),
  output_file,
  col.names=T,row.names=F,quote=F,sep="\t")
