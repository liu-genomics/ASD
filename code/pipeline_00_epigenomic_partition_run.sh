#whole genome enhancers 10,000bp
sh ../lib/160920_partition_whole_genome_to_regions_within_distance_to_genes.sh ../other_annotation/epigenomic_annotation/Whole_genome.bed 10000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Whole_genome

# whole genome enhancers 25,000bp
sh ../lib/160920_partition_whole_genome_to_regions_within_distance_to_genes.sh ../other_annotation/epigenomic_annotation/Whole_genome.bed 25000 ../other_annotation/refseq/ ../other_annotation/epigenomic

# whole genome annotation
sh ../lib/160920_partition_promoters_and_utrs.sh ../other_annotation/epigenomic_annotation/Whole_genome.bed ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Whole_genome

# noonan union 10,000bp for enhancers
sh ../lib/160920_partition_enhancers.sh ../other_annotation/epigenomic_annotation/Noonan_union.bed 10000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Noonan_brain_union

# noonan unon 25,000bp for enhancers
sh ../lib/160920_partition_enhancers.sh ../other_annotation/epigenomic_annotation/Noonan_union.bed 25000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Noonan_brain_union

# noonan promoter, 5-utr and 3-utr
sh ../lib/160920_partition_promoters_and_utrs.sh ../other_annotation/epigenomic_annotation/Noonan_union.bed ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Noonan_brain_union

# Brain roadmap union 10,000bp for enhancers
sh ../lib/160920_partition_enhancers.sh ../other_annotation/epigenomic_annotation/brain_roadmap_union.bed 10000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Brain_roadmap_union

# Brain roadmap union 25,000bp for enhancers
sh ../lib/160920_partition_enhancers.sh ../other_annotation/epigenomic_annotation/brain_roadmap_union.bed 25000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Brain_roadmap_union

# Brain roadmap union promoter, 5-utr and 3-utr
sh ../lib/160920_partition_promoters_and_utrs.sh ../other_annotation/epigenomic_annotation/brain_roadmap_union.bed ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Brain_roadmap_union


# Roadmap DHS union 10,000bp for enhancers
sh ../lib/160920_partition_enhancers.sh ../other_annotation/epigenomic_annotation/encode_DHS_union.bed  10000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Encode_DHS_union

# Roadmap DHS union 25,000bp for enhancers
sh ../lib/160920_partition_enhancers.sh ../other_annotation/epigenomic_annotation/encode_DHS_union.bed  25000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Encode_DHS_union

# Roadmap DHS union promoter, 5-utr and 3-utr
sh ../lib/160920_partition_promoters_and_utrs.sh ../other_annotation/epigenomic_annotation/encode_DHS_union.bed  ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Encode_DHS_union


# Fantom union (extended 500bp) 10,000bp for enhancers
sh ../lib/160920_partition_enhancers.sh ../other_annotation/epigenomic_annotation/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt_fetal_brain_160721.list_0.01_mode2.bed.extend.500.bed  10000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Fantom_fetal_brain_500bp_union

# Fantom union (extended 500bp) union 25,000bp for enhancers
sh ../lib/160920_partition_enhancers.sh ../other_annotation/epigenomic_annotation/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt_fetal_brain_160721.list_0.01_mode2.bed.extend.500.bed  25000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Fantom_fetal_brain_500bp_union


# Fantom union (extended 500bp) union promoter, 5-utr and 3-utr
sh ../lib/160920_partition_promoters_and_utrs.sh ../other_annotation/epigenomic_annotation/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt_fetal_brain_160721.list_0.01_mode2.bed.extend.500.bed ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Fantom_fetal_brain_500bp_union


# Noonan brain and roadmap union 10,000bp for enhancers
sh ../lib/160920_partition_enhancers.sh ../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed 10000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Noonan_brain_roadmap_union

# Noonan brain and roadmap union 25,000bp for enhancers
sh ../lib/160920_partition_enhancers.sh ../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed  25000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Noonan_brain_roadmap_union


# Noonan brain and roadmap union promoter, 5-utr and 3-utr
sh ../lib/160920_partition_promoters_and_utrs.sh ../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Noonan_brain_roadmap_union


# Psychencode_yale_ASD_DFC_H3K27ac_union.bed promoter, 5-utr and 3-utr
sh ../lib/160920_partition_promoters_and_utrs.sh ../other_annotation/epigenomic_annotation/Psychencode_yale_ASD_DFC_H3K27ac_union.bed ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Yale_ASD_DFC_H3K27ac


# Psychencode_yale_ASD_DFC_H3K27ac_union.bed promoter, 5-utr and 3-utr
sh ../lib/160920_partition_promoters_and_utrs.sh ../other_annotation/epigenomic_annotation/Psychencode_yale_ASD_CBC_H3K27ac_union.bed ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Yale_ASD_CBC_H3K27ac

#Differentillay regulated fetal brain H3K27ac sites from Sun et al., 2016

sh ../lib/160920_partition_promoters_and_utrs.sh ../other_annotation/epigenomic_annotation/ASD_K27ac_CB_down_Sun_et_al_Cell_2016.bed ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ ASD_K27ac_CB_down_Sun_et_al_Cell_2016
sh ../lib/160920_partition_promoters_and_utrs.sh ../other_annotation/epigenomic_annotation/ASD_K27ac_CB_up_Sun_et_al_Cell_2016.bed ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ ASD_K27ac_CB_up_Sun_et_al_Cell_2016
sh ../lib/160920_partition_promoters_and_utrs.sh ../other_annotation/epigenomic_annotation/ASD_K27ac_PFC_down_Sun_et_al_Cell_2016.bed ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ ASD_K27ac_PFC_down_Sun_et_al_Cell_2016
sh ../lib/160920_partition_promoters_and_utrs.sh ../other_annotation/epigenomic_annotation/ASD_K27ac_PFC_up_Sun_et_al_Cell_2016.bed ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ ASD_K27ac_PFC_up_Sun_et_al_Cell_2016
sh ../lib/160920_partition_promoters_and_utrs.sh ../other_annotation/epigenomic_annotation/ASD_K27ac_TC_down_Sun_et_al_Cell_2016.bed ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ ASD_K27ac_TC_down_Sun_et_al_Cell_2016
sh ../lib/160920_partition_promoters_and_utrs.sh ../other_annotation/epigenomic_annotation/ASD_K27ac_TC_up_Sun_et_al_Cell_2016.bed ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ ASD_K27ac_TC_up_Sun_et_al_Cell_2016
