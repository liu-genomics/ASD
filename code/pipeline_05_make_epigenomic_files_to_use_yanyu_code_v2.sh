# This script will generate epigenomic partition files in a format that is usable by yanyu's code/.
# flag set to 1001, using yanyu's pipeline to partition enhancers

sh ../lib/160927_transform_to_yanyu_format_v2.sh Brain_roadmap_union ../other_annotation/refseq/ref_id_to_genename_table.txt \
../other_annotation/epigenomic_annotation/ 1001 brain_roadmap_union.bed

sh ../lib/160927_transform_to_yanyu_format_v2.sh Encode_DHS_union ../other_annotation/refseq/ref_id_to_genename_table.txt \
../other_annotation/epigenomic_annotation/ 1001 encode_DHS_union.bed

sh ../lib/160927_transform_to_yanyu_format_v2.sh Fantom_fetal_brain_500bp_union ../other_annotation/refseq/ref_id_to_genename_table.txt \
../other_annotation/epigenomic_annotation/ 1001 fantom_fetal_brain_extend_500.bed

sh ../lib/160927_transform_to_yanyu_format_v2.sh Noonan_brain_roadmap_union ../other_annotation/refseq/ref_id_to_genename_table.txt \
../other_annotation/epigenomic_annotation/ 1001 Noonan_brain_roadmap_union.bed

sh ../lib/160927_transform_to_yanyu_format_v2.sh Noonan_brain_union ../other_annotation/refseq/ref_id_to_genename_table.txt \
../other_annotation/epigenomic_annotation/ 1001 Noonan_union.bed

sh ../lib/160927_transform_to_yanyu_format_v2.sh Whole_genome ../other_annotation/refseq/ref_id_to_genename_table.txt \
../other_annotation/epigenomic_annotation/ 1001 whole_genome.bed

sh ../lib/160927_transform_to_yanyu_format_v2.sh Noonan_brain_roadmap_intersect ../other_annotation/refseq/ref_id_to_genename_table.txt \
../other_annotation/epigenomic_annotation/ 1001 Noonan_brain_roadmap_intersect.bed