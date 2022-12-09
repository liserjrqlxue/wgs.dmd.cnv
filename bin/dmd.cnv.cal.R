library(DNAcopy)
args<-commandArgs(T)
sampleID<-args[1]
gender<-args[2]
depth.X<-as.numeric(args[3])
depth<-args[4]
outdir<-args[5]

bin<-50
threshold<-50 # 50%
fixPath<-"fix.mean.txt"
fix<-read.table(fixPath,stringsAsFactors=F,header=T)

raw<-read.table(depth,stringsAsFactors=F)
colnames(raw)<-c("chr","pos","depth")

raw$ratio=raw$depth/depth.X
raw$factor=fix$mean
raw$fixRatio=raw$ratio/raw$factor

calCV<-function(x){
	x.mean<-mean(x)
	x.sd<-sd(x)
	x.cv<-x.sd/x.mean
	y<-c(x.mean,x.sd,x.cv)
	return(y)
}
qc.raw<-calCV(raw$ratio)
qc.fix<-calCV(raw$fixRatio)
message(paste(sampleID,"cv:",qc.raw[3],qc.fix[3],sep="\t"))

#cna<-segment(smooth.CNA(CNA(log2(raw$fixRatio),raw$chr,raw$pos,data.type = "logratio",sampleid=paste0("ID.",sampleID))))

n<-floor(nrow(raw)/bin)
bins<-raw[1:n,]
bins$end<-bins$pos
for(j in 1:n){
	index=(j-1)*bin+1
	region=(j-1)*bin+1:bin
	bins[j,]$chr<-raw[index,]$chr
	bins[j,]$pos<-raw[index,]$pos
	bins[j,]$depth<-mean(raw[region,]$depth)
	bins[j,]$ratio<-mean(raw[region,]$ratio)
	bins[j,]$factor<-mean(raw[region,]$factor)
	bins[j,]$fixRatio<-mean(raw[region,]$fixRatio)
	bins[j,]$end<-raw[j*bin,]$pos
}

cna<-segment(smooth.CNA(CNA(log2(bins$fixRatio),bins$chr,bins$pos,data.type = "logratio",sampleid=paste0("ID.",sampleID))))

output<-cna$output
output$ID=sampleID
output$ratio<-2^output$seg.mean

write.table(output,paste0(outdir,"/",sampleID,".DMD.bin",bin,".txt"),col.names=T,row.names=F,sep="\t",quote=F)
