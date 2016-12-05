# code to transform the epigenomic partitions from an epigenomic dataset (e.g., Noonan_brain_union.bed)to Yanyu's format 


# sh 160927_transform_to_yanyu_format.sh [epigenomic_prefix] [reference] [epigenomic_folder] [transform_flag]
# [epigenomic_prefix] is the prefix of a set of epigenomic partition files from one epigenomic dataset, e.g., Noonan_brain_roadmap_union
# [reference] is ../other_annotation/refseq/ref_id_to_genename_table.txt
# [epigenomic_folder] is the fold that has all the epigenomic partition files, it is also the folder that a subfolder of transformed files will be created
# [transform_flag] is a identifier to the transformation method. Could be different based on how enhancers are partitioned.
# e.g., whether I use Brain_roadmap_union.enhancers.10000.bp_within_TSS.bed or Brain_roadmap_union.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed
# [output_postfix] is the epigenomic identifer that will be consistent with yanyu's annotations. e.g., brain_roadmap_union.bed

epigenomic_prefix=$1 # e.g., Noonan_brain_roadmap_union
reference=$2 # 
epigenomic_folder=$3
transform_flag=$4
output_postfix=$5

Rscript ../lib/transform_bed_regions_to_yanyu_format.R $reference ${epigenomic_folder} ${epigenomic_prefix}.enhancers.10000.bp_within_TSS.bed ${epigenomic_prefix}_yanyu_format_${transform_flag} 1 rest_intersect__0__10000__ ${output_postfix}
Rscript ../lib/transform_bed_regions_to_yanyu_format.R $reference ${epigenomic_folder} ${epigenomic_prefix}.enhancers.25000.bp_within_TSS.bed ${epigenomic_prefix}_yanyu_format_${transform_flag} 1 rest_intersect__10000__20000__ ${output_postfix}
Rscript ../lib/transform_bed_regions_to_yanyu_format.R $reference ${epigenomic_folder} ${epigenomic_prefix}.promoter_1kb_genename.bed ${epigenomic_prefix}_yanyu_format_${transform_flag} 0 intersect__refseq_id_promoter_1kb_upstream_pair_table.txt__ ${output_postfix}
Rscript ../lib/transform_bed_regions_to_yanyu_format.R $reference ${epigenomic_folder} ${epigenomic_prefix}.5utr_genename.bed ${epigenomic_prefix}_yanyu_format_${transform_flag} 0 intersect__refseq_id_5utr_pair_table.txt__ ${output_postfix}
Rscript ../lib/transform_bed_regions_to_yanyu_format.R $reference ${epigenomic_folder} ${epigenomic_prefix}.3utr_genename.bed ${epigenomic_prefix}_yanyu_format_${transform_flag} 0 intersect__refseq_id_3utr_pair_table.txt__ ${output_postfix}


