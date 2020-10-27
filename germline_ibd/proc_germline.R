library(data.table)

args=commandArgs(TRUE)

infile = args[1]
popfile = args[2]
nsnps = args[3]
outfile=args[4]

dat<-fread(infile)

colnames(dat)=c("f1","f2","chr","position","id2","total_snps","glength","unit","mismatch","ind1hom","ind2hom")

pop = fread(popfile, header=FALSE)
colnames(pop) = c("FID","IID")
ninds = nrow(pop)
pop$sno = seq(1,ninds)

chromosome=dat$chr[1]
dat[,c("fid1", "iid1") := tstrsplit(f1, " ", fixed=TRUE)][,f1:=NULL]
dat[,c("fid2", "iid2") := tstrsplit(f2, " ", fixed=TRUE)][,f2:=NULL]
#dat[,c("tsk", "iid1") := tstrsplit(iid1, "_", fixed=TRUE)]
#dat[,c("tsk", "iid2") := tstrsplit(iid2, "_", fixed=TRUE)]
dat[,c("start_bp", "stop_bp") := tstrsplit(position, " ", fixed=TRUE)][,position:=NULL]
#dat[,c("start_rsid", "stop_rsid") := tstrsplit(id2, " ", fixed=TRUE)][,id2:=NULL]
dat = merge(dat, pop, by.x =c("fid1","iid1"),
            by.y = c("FID","IID"))
colnames(dat)[15]="sno1"

dat = merge(dat, pop, by.x =c("fid2","iid2"),
            by.y = c("FID","IID"))
colnames(dat)[16]="sno2"

dat[,"pair":=paste(sno1,sno2,sep="_")]

dat2 = dat[,.(total_length=sum(glength)/200),
           by=c("pair",'sno1','sno2')]
dat2=dat2[,c('sno1','sno2','pair','total_length'),with=FALSE]


all_combs=t(combn(pop$sno,2))
all_combs=data.table(sno1 = all_combs[,1],
                     sno2 = all_combs[,2])

all_combs[, pair := paste(sno1,sno2,sep="_")]
all_combs = merge(all_combs,
                  dat2[, c("pair","total_length")],
                  by = "pair",
                  sort = FALSE,
                  all.x = TRUE)

all_combs[is.na(all_combs$total_length),
          total_length := 0]

all_combs[,nsnp := nsnps]

samz=data.table(pair = paste(pop$sno,
                             pop$sno,
                             sep="_"),
                sno1 = pop$sno,
                sno2 = pop$sno,
                total_length = 1,
                nsnp = nsnps)

all_combs=rbind(all_combs,samz)

all_combs=all_combs[order(sno2,sno1)]

#write grm
fwrite(all_combs[,.(sno2,sno1,nsnp,total_length)],
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



