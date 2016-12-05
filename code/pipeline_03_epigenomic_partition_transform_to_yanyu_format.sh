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


# move to a further distance, 20kb and 20kb-to 100kb 

#whole genome enhancers  20kb and 20kb-to 100kb 
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/Whole_genome.bed 20000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Whole_genome 100000


# noonan union  20kb and 20kb-to 100kb 
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/Noonan_union.bed 20000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Noonan_brain_union 100000


# Brain roadmap union  20kb and 20kb-to 100kb 
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/brain_roadmap_union.bed 20000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Brain_roadmap_union 100000

# Roadmap DHS union  20kb and 20kb-to 100kb 
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/encode_DHS_union.bed  20000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Encode_DHS_union 100000

# Fantom union (extended 500bp)  20kb and 20kb-to 100kb 
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt_fetal_brain_160721.list_0.01_mode2.bed.extend.500.bed  \
20000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Fantom_fetal_brain_500bp_union 100000

# Noonan_brain_roadmap_union  20kb and 20kb-to 100kb 
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed  \
20000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Noonan_brain_roadmap_union 100000

# Noonan_brain_roadmap_intersect  20kb and 20kb-to 100kb 
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_intersect.bed  \
20000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Noonan_brain_roadmap_intersect 100000


#whole genome enhancers  20kb and 20kb-to 500kb 
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/Whole_genome.bed 20000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Whole_genome 500000


# noonan union  20kb and 20kb-to500kb
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/Noonan_union.bed 20000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Noonan_brain_union 500000


# Brain roadmap union  20kb and 20kb-to 500kb 
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/brain_roadmap_union.bed 20000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Brain_roadmap_union 500000

# Roadmap DHS union  20kb and 20kb-to 500kb 
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/encode_DHS_union.bed  20000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Encode_DHS_union 500000

# Fantom union (extended 500bp)  20kb and 20kb-to 500kb 
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt_fetal_brain_160721.list_0.01_mode2.bed.extend.500.bed  \
20000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Fantom_fetal_brain_500bp_union 500000

# Noonan_brain_roadmap_union  20kb and 20kb-to 500kb 
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed  \
20000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Noonan_brain_roadmap_union 500000

# Noonan_brain_roadmap_intersect  20kb and 20kb-to 500kb
sh ../lib/160920_partition_enhancers_for_yanyu_pipeline.sh \
../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_intersect.bed  \
20000 ../other_annotation/refseq/ ../other_annotation/epigenomic_annotation/ Noonan_brain_roadmap_intersect 500000


