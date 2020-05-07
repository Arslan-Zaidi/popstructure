

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript simphenotype_ge_3.R <frequency file> <output_file")}

library(data.table)
library(dplyr)

set.seed(123)

#frequency file
#colnames in frequency file must be ALT_FREQS
freq_file=args[1]

#output file - genetic effects
effects_file = args[2]


# load variant frequency file
p = fread(freq_file)

#calculate the independent component of variance required
sigma2_l=0.8/sum(sapply(p$ALT_FREQS,function(x){
  beta=(2*x*(1-x))^(1-0.4)
  return(beta)
}))

#sample maf-dependent effects using the model above
p$beta=sapply(p$ALT_FREQS,function(x){
  beta=rnorm( 1 , mean=0, sd=sqrt(sigma2_l * (2*x*(1-x))^-0.4 ))
})

#let's calculate sigma2_g to confirm that the total genetic variance is indeed 0.8
sigma2_g=sum(mapply(function(b,p){ b^2* 2*p*(1-p) }, p$beta, p$ALT_FREQS))


#save the effect sizes to file and use plink2 to generate PRS
fwrite(p%>%
         select(ID,ALT,beta),
           effects_file,
       row.names=F,col.names=F,quote=F,sep="\t")
