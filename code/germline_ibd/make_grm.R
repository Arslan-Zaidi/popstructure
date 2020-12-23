
#write grm from ibd match file
library(data.table)
args=commandArgs(TRUE)
infile = args[1]
popfile = args[2]
nsnps = args[3]
outfile = args[4]

dat = fread(infile)
colnames(dat) = c("chr","iid1","iid2","total_length")
dat[, total_length := total_length/200]

pop = fread(popfile, header = FALSE) 
colnames(pop) = c("FID","IID")
ninds = nrow(pop)
pop$iid = seq(1,ninds)

dat[, pair := paste(iid1,iid2,sep="_")]
dat2 = dat[,sum(total_length),by=c("pair","iid1","iid2")]

all_combs=t(combn(pop$iid,2))
all_combs=data.table(iid1 = all_combs[,1],
                     iid2 = all_combs[,2])
all_combs[, pair := paste(iid1,iid2,sep="_")]

all_combs = merge(all_combs,
                  dat2[, c("pair","V1")],
                  by = "pair",
                  sort = FALSE,
                  all.x = TRUE)

all_combs[is.na(all_combs$V1),
          V1 := 0]
all_combs[,iid1 := iid1+1]
all_combs[,iid2 := iid2+1]
all_combs[,nsnp := nsnps]

samz=data.table(pair = paste(pop$iid+1,
                             pop$iid+1,
                             sep="_"),
                iid1 = pop$iid+1,
                iid2 = pop$iid+1,
                V1 = 1,
                nsnp = nsnps)

all_combs=rbind(all_combs,samz)

all_combs=all_combs[order(iid1,iid2)]

#write grm
fwrite(all_combs[,.(iid1,iid2,nsnp,V1)],
       paste(outfile,".grm", sep=""),
       col.names=FALSE,
       row.names=FALSE,
       quote=FALSE,sep="\t")

#write grm id file
fwrite(data.table(fid=pop$FID,
                  iid=pop$IID),
       paste(outfile,".grm.id", sep=""),
       sep="\t",
       col.names=FALSE,
       row.names=FALSE,
       quote=FALSE)



