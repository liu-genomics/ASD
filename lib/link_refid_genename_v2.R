#!/usr/bin/Rscript
# Rscript link_mutation_with_gene.R ref_id_genename_table.txt genomic_feature_overlap_with_mutation.bed (output from bedtools intersect -a -b -wa -wb) output_prefix
args<-commandArgs(TRUE)
ref_id_genename=read.table(args[1],header=FALSE,sep="\t",stringsAsFactors=FALSE)
data=read.table(args[2],header=FALSE,sep="\t",stringsAsFactors=FALSE)
colnames(ref_id_genename)=as.character(seq(1:2))
colnames(data)=as.character(seq(1:9))
combine=merge(data,ref_id_genename,by.x="8",by.y="1")
combine_temp=combine[,c(2:4,10,5)]
combine_temp=combine_temp[!duplicated(combine_temp),]
combine_temp=combine_temp[combine_temp[,2]!="n/a",]
combine_temp[,5]=paste("E",combine_temp[,5],sep="")
write.table(combine_temp,args[3],col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
quit(save="no")
