
#write grm from ibd match file
library(data.table)
args=commandArgs(TRUE)
infile=args[1]
outfile=args[2]

dat=fread(infile)
colnames(dat)=c("chr","iid1","iid2","total_length")
dat[,total_length:=total_length/200]


dat[,pair:=paste(iid1,iid2,sep="_")]
dat2=dat[,sum(total_length),by=c("pair","iid1","iid2")]

all_combs=t(combn(0:8999,2))
all_combs=data.table(iid1=all_combs[,1],iid2=all_combs[,2])
all_combs[,pair:=paste(iid1,iid2,sep="_")]

all_combs=merge(all_combs,dat2[,c("pair","V1")],by="pair",sort=F,all.x=T)
all_combs[is.na(all_combs$V1),V1:=0]
all_combs[,iid1:=iid1+1]
all_combs[,iid2:=iid2+1]
all_combs[,nsnp:=15e4]

samz=data.table(pair=paste((1:9000),c(1:9000),sep="_"),iid1=1:9000,iid2=1:9000,V1=1,nsnp=150e3)

all_combs=rbind(all_combs,samz)

all_combs=all_combs[order(iid2,iid1)]

#write grm
fwrite(all_combs[,.(iid2,iid1,nsnp,V1)],
       paste(outfile,".grm"),
       col.names=F,row.names=F,quote=F,sep="\t")

#write grm id file
fwrite(data.table(fid=paste("msp_",(1:9000),sep=""),
                  iid=paste("msp_",(1:9000),sep="")),
       paste(outfile,".grm.id"),
       sep="\t",
       col.names=F,row.names=F,quote=F)



