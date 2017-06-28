### Extract the mutation rate for each window of an epigenomic annotation and get the sum mutation rate, for partitioning the non-coding risk contributing to ASD
awk {'print $1"\t"$2"\t"$3"\t"NR'} ../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.promoter_1kb_genename.bed > \
../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.promoter_1kb_genename.bed.tmp

/media/yuwen/F/tools/bigWigAverageOverBed ../other_annotation/mutation_rate/Daly_mutrate.bw ../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.promoter_1kb_genename.bed.tmp \
../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.promoter_1kb_genename.bed.mutrate

rm ../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.promoter_1kb_genename.bed.tmp

awk {'print $1"\t"$2"\t"$3"\t"NR'} ../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed > \
../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.tmp

/media/yuwen/F/tools/bigWigAverageOverBed ../other_annotation/mutation_rate/Daly_mutrate.bw ../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.tmp \
../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.mutrate

rm ../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.tmp