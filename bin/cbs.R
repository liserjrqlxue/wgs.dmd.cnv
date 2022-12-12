library(DNAcopy)
args<-commandArgs(T)
input<-args[1]
output<-args[2]

a<-read.table(input,stringsAsFactors=T,sep="\t",fill=T,header=T)
cna<-segment(
    smooth.CNA(
        CNA(
            log2(a$fixRatio+0.00001),
            a$chr,
            a$start,
            data.type = "logratio",
        )
    )
)
cnv<-cna$output
write.table(cnv,output,col.names=T,row.names=F,sep="\t",quote=F)
