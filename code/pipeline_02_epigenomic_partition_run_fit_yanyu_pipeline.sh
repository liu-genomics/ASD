# This script will partition enhancers within (rather than overlap with) certain distancer range of gene TSSs
# also 0-1000bp and 1000-2000 bp (for example) will be disjoint


#whole genome enhancers 10,000bp, 10,000-20,000bp
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/Whole_genome.bed 10000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Whole_genome 20000


# noonan union 10,000bp for enhancers, 10,000-20,000bp
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/Noonan_union.bed 10000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Noonan_brain_union 20000


# Brain roadmap union 10,000bp for enhancers, 10,000-20,000bp
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/brain_roadmap_union.bed 10000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Brain_roadmap_union 20000

# Roadmap DHS union 10,000bp for enhancers, 10,000-20,000bp
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/encode_DHS_union.bed  10000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Encode_DHS_union 20000

# Fantom union (extended 500bp) 10,000bp for enhancers, 10,000-20,000bp
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt_fetal_brain_160721.list_0.01_mode2.bed.extend.500.bed  \
10000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Fantom_fetal_brain_500bp_union 20000

# Noonan_brain_roadmap_union 10,000bp for enhancers, 10,000-20,000bp
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed  \
10000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Noonan_brain_roadmap_union 20000

# Noonan_brain_roadmap_intersect 10,000bp for enhancers, 10,000-20,000bp
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_intersect.bed  \
10000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Noonan_brain_roadmap_intersect 20000

# Psychencode_yale_ASD_DFC_H3K27ac_union.bed 10,000bp for enhancers, 10,000-20,000bp
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/Psychencode_yale_ASD_DFC_H3K27ac_union.bed  \
10000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Yale_ASD_DFC_H3K27ac 20000

# Psychencode_yale_ASD_CBC_H3K27ac_union.bed 10,000bp for enhancers, 10,000-20,000bp
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/Psychencode_yale_ASD_CBC_H3K27ac_union.bed  \
10000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Yale_ASD_CBC_H3K27ac 20000
