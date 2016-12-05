# code to transform the epigenomic partitions from an epigenomic dataset (e.g., Noonan_brain_union.bed)to Yanyu's format 


# Rscript 160927_transform_to_yanyu_format.R [epigenomic_prefix] [reference] [epigenomic_folder] [transform_flag]
# [epigenomic_prefix] is the prefix of a set of epigenomic partition files from one epigenomic dataset, e.g., Noonan_brain_roadmap_union
# [reference] is ../other_annotation/refseq/ref_id_to_genename_table.txt
# [epigenomic_folder] is the fold that has all the epigenomic partition files, it is also the folder that a subfolder of transformed files will be created
# [transform_flag] is a identifier to the transformation method. Could be different based on how enhancers are partitioned.
# e.g., whether I use Brain_roadmap_union.enhancers.10000.bp_within_TSS.bed or Brain_roadmap_union.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed
# 

source("../lib/transform_bed_regions_to_yanyu_format.R")

args<-commandArgs(TRUE, trailingOnly = TRUE)
print(args)
epigenomic_prefix <- args[1] # e.g., Noonan_brain_roadmap_union
reference <- args[2] # 
epigenomic_folder <- args[3]
transform_flag <- args[4]

