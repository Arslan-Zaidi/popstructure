
#!/usr/bin/env Rscript

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript edit_fam.R <fam_file> <.sample file> <output_file>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
}))

fam_file=args[1] # plink .fam file
sample_file=args[2] # .sample file containing family info

fam = fread(fam_file, header=F)
colnames(fam) = c("FID","IID","Mother","Father","Sex","Phenotype")
fam = fam %>%
  select(IID,Sex,Phenotype)


sam = fread(sample_file, header=T, sep=" ")

# create unique list of parents
# this will be used to give the same family ID to siblings
family_ids = data.table( group = unique(sam$group) )
family_ids$FID = seq( 1, nrow(family_ids) )

#add the newly generated family IDs to sam file
sam = merge(sam, family_ids, by="group", all.x=T, sort=F)
sam = sam %>%
  separate(group, into=c("Mother","Father"), sep=",")

#add family info to fam file
fam = merge(fam,
            sam %>%
              select(FID,sample,Mother,Father),
            all.x = T, sort = F, by.x = "IID", by.y = "sample")

#reorder columns so that they are in the same order as the plink .fam file
fam = fam %>%
  mutate(FID=)
  select(FID,IID,Mother,Father,Sex,Phenotype)

fwrite(fam, fam_file, col.names = F, row.names=F, quote=F,sep="\t")
fwrite(sam, sample_file, col.names = F, row.names=F, quote=F,sep="\t")
