#!/usr/bin/Rscript
# Rscript script.R [ref_id_genename_table.txt] [bedfile] [output_name]
# [ref_id_genename] is ../other_annotation/refseq/ref_id_to_genename_table.txt
# [bedfile] is a regular bed file, format: chr \t start \t end \t refseqID
args<-commandArgs(trailingOnly = TRUE)
ref_id_genename=read.table(args[1],header=FALSE,sep="\t",stringsAsFactors=FALSE)
data=read.table(args[2],header=FALSE,sep="\t",stringsAsFactors=FALSE)
colnames(ref_id_genename)=c("refseq_id", "genename")
colnames(data)=c("chr", "start", "end", "refseq_id")
combine=merge(data,ref_id_genename,by.x="refseq_id",by.y="refseq_id")
combine = combine[combine$genename != "n/a",]
write.table(combine[,c("chr","start","end","genename")],args[3],col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
quit(save="no")