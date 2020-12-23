
#script takes in genetic values and generates phenotypes with and without stratification

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript simphenotype_ge.R <gvalue file> <sample_file file> <output_file> <seed>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

gvalue_file = args[1] #gvalue file should contain gvalues for all chromosomes
samfile = args[2] # .sample file carrying deme and family info for each sib pair
#popfile = args[3] # pop file from msprime simulation
output_file = args[3] #name of output file
set.seed(args[4])

prs=fread(gvalue_file,header=F)

colnames(prs)<-c("iteration","FID","IID","chromosome","prs")


#sum genetic value across chromosomes to get one for each individual
prs = prs%>%
  group_by(FID,IID,iteration)%>%
  summarize(prs = sum(prs))%>%
  ungroup()

sample_size=nrow(prs)

#load file containing sample ID, deme id and latitude and longitude for each individual
sample_file = fread(samfile,header=T,sep=" ")
colnames(sample_file) = c("FID","IID","population","longitude","latitude","parents")

#add this info to prs file
prs=merge(prs, sample_file, by=c("FID","IID"), sort=F)

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
prs$sharp = sapply(prs$population,
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
prs$grandom = as.numeric(scale( prs$grandom, scale = T) * sqrt( 1 - 0.8))
prs$sharp = as.numeric(scale( prs$sharp , scale = T) * sqrt( 1 - 0.8))
prs$smooth = as.numeric(scale( prs$smooth, scale = T) * sqrt(1 - 0.8))
prs$smooth_long = as.numeric(scale( prs$smooth_long, scale = T) * sqrt(1 - 0.8))

#add prs to each of the environmental effects
prs = prs %>%
  mutate(grandom = prs + grandom,
         smooth = prs + smooth,
         smooth_long = prs + smooth_long,
         sharp = prs + sharp)

fwrite(
  prs%>%
    select(FID,IID,grandom,smooth,smooth_long,sharp),
  output_file,
  col.names=T,row.names=F,quote=F,sep="\t")
