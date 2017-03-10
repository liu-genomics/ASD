# will reassign genes to 50bp window. 
# The assginment will be done for each window.
# Previously, the assignment (../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename.bed) is based on bigger intervals
# which were then chopped into 50-bp windows

#remove genes not assembled to 23 chromosomes and transform to BED format
grep -vP "^.*\t.*\t.*_.*\t" ../other_annotation/refseq/tss_position.txt | awk {'print $3"\t"$4"\t"$4"\t"$2'} > ../other_annotation/refseq/tss_position.txt.temp
#adjust distance
bedtools flank -i ../other_annotation/refseq/tss_position.txt.temp -g ../other_annotation/genome_build/hg19.genome -l 10000 -r 10000 > ../other_annotation/refseq/tss_position.txt.temp2
grep -v "n/a" ../other_annotation/refseq/tss_position.txt.temp2 > ../other_annotation/refseq/tss_position.txt.temp3

# assign non-coding (nf) 50bp windows to nearest TSSs
bedtools intersect -a ../other_annotation/epigenomic_annotation/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.50bp_window.bed \
-b ../other_annotation/refseq/tss_position.txt.temp3 -wa -wb | awk {'print $1"\t"$2"\t"$3"\t"$8"\t"$4'} > ../other_annotation/refseq/tss_position.txt.temp4

perl ../lib/160920_find_nearest_gene_for_each_intergenic_enhancer.pl ../other_annotation/refseq/tss_position.txt.temp4 ../other_annotation/refseq/ ../other_annotation/refseq/tss_position.txt.temp5

awk {'print $5"\t"$4'} ../other_annotation/refseq/tss_position.txt.temp5 > \
../other_annotation/epigenomic_annotation/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.50bp_window.bed_genename_corrected.bed # now the distance to TSS is based on each 50bp, rather than a bigger interval

# add up different parts of genename assignment together with the corrected nf gene assignment
cat ../other_annotation/epigenomic_annotation/Whole_genome.promoter_1kb_genename.bed.50bp_window.bed_genename.bed \
../other_annotation/epigenomic_annotation/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.50bp_window.bed_genename_corrected.bed \
../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.50bp_window.bed_genename.bed > \
../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename_corrected.bed
