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
exon<-read.table("dmd.exon.bed",stringsAsFactors=F)
colnames(exon)<-c("chr","start","end","length","gene","trans","strand","exon")
exon$start<-exon$start+1

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
if(gender=="M"){
	output$percent<-abs(output$ratio-1)*100
}else{
	output$percent<-abs(output$ratio-1)*2*100
}

write.table(output,paste0(outdir,"/",sampleID,".DMD.bin",bin,".txt"),col.names=T,row.names=F,sep="\t",quote=F)

output$length<-0
output$factor<-0
output$anno<-""
output$coverage<-""

filter<-output[output$percent>=threshold,]
# annotation
if(nrow(filter)>0){
	for(j in 1:nrow(filter)){
		filter[j,]$percent<-min(filter[j,]$percent,100)
		filter[j,]$loc.end<-max(bins[bins$pos<=filter[j,]$loc.end,]$end)
		start<-filter[j,]$loc.start
		end<-filter[j,]$loc.end
		filter[j,]$length<-end-start+1
		filter[j,]$factor<-mean(bins[bins$pos>=start&bins$pos<=end,]$factor)

		t.exon<-exon[exon$end>=start&exon$start<=end,]
		anno<-c()
		coverage<-c()
		if(nrow(t.exon)>0){
			for(k in 1:nrow(t.exon)){
				anno<-append(anno,t.exon[k,]$exon)
				t.start<-max(start,t.exon[k,]$start)
				t.end<-min(end,t.exon[k,]$end)
				coverage<-append(coverage,(t.end-t.start+1)/(t.exon[k,]$end-t.exon[k,]$start+1))
			}
			filter[j,]$anno<-paste(anno,collapse=",")
			filter[j,]$coverage<-paste(coverage,collapse=",")
		}
		
	}
}
write.table(filter,paste0(outdir,"/",sampleID,".DMD.bin",bin,".threshold",threshold,".txt"),col.names=T,row.names=F,sep="\t",quote=F)
write.table(filter[filter$anno!="",],paste0(outdir,"/",sampleID,".DMD.bin",bin,".threshold",threshold,".EXON.txt"),col.names=T,row.names=F,sep="\t",quote=F)
