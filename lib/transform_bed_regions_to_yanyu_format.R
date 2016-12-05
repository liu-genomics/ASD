#!/usr/bin/Rscript
# Rscript transform_bed_regions_to_yanyu_format.R [ref_id_genename_table.txt] [enhancer_partitions] [output_folder] [mode] [output_label] [output_postfix]
# [ref_id_genename_table.txt] is ../other_annotation/refseq/ref_id_to_genename_table.txt
# [enhancer_partition] is the files generated from pipeline_00_epigenomic_partition_run.sh, could be enhancer partitions or promoter/5'utr/3'utr partitions
# [output_folder_label] is the folder that holds the output files
# [mode] if enhancer partitions mode =1; other wise mode = 0
# [output_label] is "intersect__refseq_id_3utr_pair_table.txt__", 
#                    "intersect__refseq_id_5utr_pair_table.txt__", 
#                    "intersect__refseq_id_promoter_1kb_upstream_pair_table.txt__",
#                     "rest_intersect__0__10000__"
#                     "rest_intersect__10000__20000__"
# [output_postfix] is the epigenomic dataset name appended to the [output_label]
# The aim is to tranform the format of each epigenomic partition file, in the final output file, there are 4 columns : chr \t start \t end \t refseq ID
# the input [enhancer_partition] has gene names as annotations. So essentially, this code translates genenames to refseq IDs. 
# For now, I will transform each gene name to the first refseq ID. 
# the downstream analysis is not performed at transcript level, so for each genename 

args<-commandArgs(trailingOnly = TRUE)
print(args)
reference <- args[1]
epigenomic_folder <-args[2]
epigenomic_partition <- args[3]
output_folder_label <- args[4]
mode <- args[5]
output_label <- args[6]
output_postfix <- args[7]

 

ref_id_genename <- read.table(reference,header=FALSE,sep="\t",stringsAsFactors=FALSE)
data <- read.delim(paste(epigenomic_folder,epigenomic_partition,sep = "/"),header=FALSE,sep="\t",stringsAsFactors=FALSE)

find_1st_refseq_id <- function(genename){
  refseq <- ref_id_genename[ref_id_genename[,2] == genename & !is.na(ref_id_genename[,2]),1]
  refseq[1]
} 

if (mode == 1){
  data$V6 <- mapply(find_1st_refseq_id, data$V4)
  data <- data.frame(data[,1:3], data[,6], data[,4], data[,5])
} else if (mode == 0){
  data$V5 <- mapply(find_1st_refseq_id, data$V4)
  data <- data.frame(data[,1:3], data[,5])
}

output_folder <- paste(epigenomic_folder,output_folder_label,sep = "/")
system(paste("mkdir -p ",output_folder))
output_file_name <- paste(output_folder,"/",output_label, output_postfix, sep = "")
write.table(data, output_file_name, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
quit(save="no")
