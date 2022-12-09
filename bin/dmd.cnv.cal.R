library(DNAcopy)
args<-commandArgs(T)
sampleID<-args[1]
input<-args[2]
output<-args[3]

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

a<-read.table(input,stringsAsFactors=T,sep="\t",fill=T)
cna<-segment(smooth.CNA(CNA(log2(a$fixRatio),a$chr,a$start,data.type = "logratio",sampleid=paste0("ID.",sampleID))))
cnv<-cna$output
cnv$ID=sampleID
cnv$ratio<-2^output$seg.mean
write.table(cnv,output,col.names=T,row.names=F,sep="\t",quote=F)
