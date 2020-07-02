
library(ggplot2)
library(dplyr)
library(data.table)
library(rprojroot)

F=is_rstudio_project$make_fix_file()

bplace_summary = fread(F("popstruct_scripts/simulating_genotypes/ukb/ukb_bplace_summary.txt"))

#sample individuals from each deme based on their frequency in the UKB
bplace_summary = bplace_summary%>%
  arrange(-North)%>%
  mutate(ninds2 = round(ninds/sum(ninds)*17499))

#generate popfile such that individuals are sampled uniformly from each deme
popfile = data.table(
  IID = paste("tsk_",0:17499,sep=""),
  deme = rep( bplace_summary$nuts215cd , each = 500),
  longitude=rep(bplace_summary$East, each = 500),
  latitude=rep(bplace_summary$North, each = 500)
)


#convert no. of individuals in each deme to an even number so they can be split evenly between test and training set
is.even = function(x){
  if(x%%2==0){return(x)}else{return(x+1)}
}

bplace_summary$ninds2=sapply(bplace_summary$ninds2,is.even)


popfile.w = data.table(
  IID = paste("tsk_",0:(sum(bplace_summary$ninds2)-1),sep=""),
  deme = rep( bplace_summary$nuts215cd ,  bplace_summary$ninds2),
  longitude=rep(bplace_summary$East,  bplace_summary$ninds2),
  latitude=rep(bplace_summary$North,  bplace_summary$ninds2)
)


fwrite(popfile,
       F("gwas/ukb/popfiles/ukb_ss500_d35.uniform.pop"),
       col.names = TRUE,
       row.names=FALSE,
       quote=FALSE,
       sep="\t")

fwrite(popfile.w,
       F("gwas/ukb/popfiles/ukb_ss500_d35.weighted.pop"),
       col.names = TRUE,
       row.names=FALSE,
       quote=FALSE,
       sep="\t")

#function to split file into training and test sets and write ids to file
split.tt = function(x,w="weighted"){
  x.train = x%>%
    filter(IID%in%paste("tsk_",seq(0 , nrow(popfile)-1,2 ),sep=""))%>%
    mutate(FID=IID)%>%
    select(FID,IID)
  
  x.test = x%>%
    filter(IID%in%paste("tsk_",seq(1 , nrow(popfile)-1,2 ),sep=""))%>%
    mutate(FID=IID)%>%
    select(FID,IID)
  
  fwrite(x.train,
         F(paste("gwas/ukb/popfiles/ukb_ss500_d35.",
                 w,".train_iid.txt",sep="")),
         col.names = FALSE,
         row.names=FALSE,
         quote=FALSE,
         sep="\t")
  
  fwrite(x.test,
         F(paste("gwas/ukb/popfiles/ukb_ss500_d35.",
                 w,".test_iid.txt",sep="")),
         col.names = FALSE,
         row.names=FALSE,
         quote=FALSE,
         sep="\t")
}

split.tt(popfile,"uniform")
split.tt(popfile,"weighted")
split.tt(popfile.w,"uniform")
split.tt(popfile.w,"weighted")

