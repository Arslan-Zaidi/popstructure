
#script will take in frequency file for the entire genotype data and:
#1. select 2000 causal variants spread 100kb apart
#2. generate effect sizes for these variants

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript simphenotype_ge_3.R <frequency file> <output_file> <seed>")}

suppressWarnings(suppressMessages({
library(data.table)
library(dplyr)
library(tidyr)
}))

#set.seed(123)

#frequency file
freq_file=args[1]

#output file - genetic effects
effects_file = args[2]

#random seed
set.seed(args[3])


# load variant frequency file
p = fread(freq_file)

colnames(p)=c("CHROM","ID","REF","ALT","ALT_FREQS","COUNT")
p=p[,c("CHROM","ID","ALT_FREQS")]
p[, c("chr", "position","ref","alt") := tstrsplit(ID, "_", fixed=TRUE)]
p = p[,c("CHROM","ID","position","ALT_FREQS")]
p$position = as.numeric(p$position)

#for each chromosome, sample the first variant from the first 100kb
#then, select every other variant to be at least 100kb apart

#write function to do this for each chromosome separately

sample.variant=function(df1){
  #sample first variant
  position1 = as.numeric( sample_n( df1[ position < 1e5, 'position' ], 1 ))
  #select all other variants to be at least 100kb apart
  #minimum positions for each window
  positions = position1 + seq(0,99)*1e5
  #pick variants that are further than these
  positions.adj = lapply( positions, function(x){
    ix = min( df1[position > x, which =TRUE ] )
    return(df1[ix])
  })
  #return datatable
  positions.adj = bind_rows(positions.adj)
  return(positions.adj)
}

#carry this out grouped by chromosome
causal.variants <- p[, sample.variant(.SD), by=CHROM]

#for some reason, sometimes the final window does not have a variant. let's remove NAs here
causal.variants = causal.variants%>%drop_na(ID)

#Now generate the effect sizes from these variants

#calculate the independent component of variance required
sigma2_l = 0.8 / sum( sapply( causal.variants$ALT_FREQS,function(x){
  beta= ( 2*x*(1-x)) ^ (1-0.4)
  return(beta)
}))

#sample maf-dependent effects using the model above
causal.variants$beta = sapply( causal.variants$ALT_FREQS , function(x){
  beta = rnorm( 1 , mean = 0, sd = sqrt(sigma2_l * (2*x*(1-x))^-0.4 ))
})

#let's calculate sigma2_g to confirm that the total genetic variance is indeed 0.8
sigma2_g = sum( mapply(function(b,p){ b^2* 2*p*(1-p) }, causal.variants$beta, causal.variants$ALT_FREQS))


#save the effect sizes to file and use plink2 to generate PRS
fwrite(causal.variants%>%
         mutate(ALT = "T")%>%
         select(ID,ALT,beta),
           effects_file,
       row.names=F,col.names=F,quote=F,sep="\t")
