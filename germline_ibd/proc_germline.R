library(data.table)

args=commandArgs(TRUE)

infile=args[1]
outfile=args[2]

dat<-fread(infile)

colnames(dat)=c("f1","f2","chr","position","id2","total_snps","glength","unit","mismatch","ind1hom","ind2hom")

chromosome=dat$chr[1]
dat[,c("fid1", "iid1") := tstrsplit(f1, " ", fixed=TRUE)][,f1:=NULL]
dat[,c("fid2", "iid2") := tstrsplit(f2, " ", fixed=TRUE)][,f2:=NULL]
dat[,c("tsk", "iid1") := tstrsplit(iid1, "_", fixed=TRUE)]
dat[,c("tsk", "iid2") := tstrsplit(iid2, "_", fixed=TRUE)]
dat[,c("start_bp", "stop_bp") := tstrsplit(position, " ", fixed=TRUE)][,position:=NULL]
dat[,c("start_rsid", "stop_rsid") := tstrsplit(id2, " ", fixed=TRUE)][,id2:=NULL]
dat[,"pair":=paste(iid1,iid2,sep="_")]

dat[,total_length:=sum(glength),by=c("iid1","iid2","pair")]
dat2=dat[,c('chr','iid1','iid2','total_length'),with=F]

fwrite(dat2,
       outfile,sep="\t",col.names=F,row.names=F,quote=F)
