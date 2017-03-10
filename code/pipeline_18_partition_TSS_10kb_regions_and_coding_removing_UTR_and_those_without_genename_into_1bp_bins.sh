
# What is different here from Pipeline 14 is that for coding I will remove regions that overlap with 5'UTR or 3'UTR and those
# that don't have a refseq gene name (For promoters I always only kept regions with a refseq genename). 

#Use predefined regions (within 10kb, removing utr, coding sequence, and promoter), calculate CG content for each. 
 
# Need to use some of the files generated in Pipeline 14. 


# The aim is to prepare data for using GLM Pisson model to calibrate mutation rates
# The coding region now is all ref-seq exons (including those with a gene name and those without, subtracted any part that overlaps with promoters)
 # bedtools subtract -a ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb.bed -b ../other_annotation/refseq/refseq_id_3_UTR_list.txt > ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR.bed

# bedtools subtract -a ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR.bed -b ../other_annotation/refseq/refseq_id_5utr_pair_table.txt > ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR.bed

#Rscript ../lib/link_refid_genename_for_regular_bed.R ../other_annotation/refseq/ref_id_to_genename_table.txt \ 
# ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR.bed 
#\../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed

sort -k1,1 -k2,2n ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed > \
../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene_temp1.bed

bedtools merge -i ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene_temp1.bed > \
../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene_temp2.bed

#also need to remove sequences from unassembled chromatin structure
bedtools makewindows -w 1 -b ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene_temp2.bed | grep -v _ |awk {'print $1"\t"$2"\t"$3"\tcoding"NR'} > ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.1bp_window.bed

mv ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.1bp_window.bed > /media/yuwen/Elements/ASD_temp_storage/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.1bp_window.bed

nohup awk {'print $1"\t"$2-25"\t"$3+25"\t"$4'} /media/yuwen/Elements/ASD_temp_storage/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.1bp_window.bed > /media/yuwen/Elements/ASD_temp_storage/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.1bp_window.bed_50bp_surrounding.bed &

nohup bedtools getfasta -tab -name -bed \
/media/yuwen/Elements/ASD_temp_storage/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.1bp_window.bed_50bp_surrounding.bed \
-fi ../other_annotation/genome_build/hg19.fasta \
-fo /media/yuwen/Elements/ASD_temp_storage/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.1bp_window.bed_50bp_surrounding.bed.fasta &

nohup Rscript ../lib/161117_calculate_cg.R \
/media/yuwen/Elements/ASD_temp_storage/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.1bp_window.bed_50bp_surrounding.bed.fasta &

cat /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_1kb_genename.bed.1bp_window.bed \
/media/yuwen/Elements/ASD_temp_storage/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.1bp_window.bed \
/media/yuwen/Elements/ASD_temp_storage/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.1bp_window.bed \
> /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed




nohup cat /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_1kb_genename.bed.1bp_window.bed_50bp_surrounding.bed.fasta.cg  \
/media/yuwen/Elements/ASD_temp_storage/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.1bp_window.bed_50bp_surrounding.bed.fasta.cg  \
/media/yuwen/Elements/ASD_temp_storage/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.1bp_window.bed_50bp_surrounding.bed.fasta.cg > \
/media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed_50bp_surrounding.bed.fasta.cg &

#mutation rate , need to split and combine
nohup split -d -l 10000000 /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed.split &

nohup sh -c 'for i in /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed.split*; do ../lib/bigWigAverageOverBed ../other_annotation/mutation_rate/Daly_mutrate.bw $i $i.mutrate; done' &

nohup cat /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed.split*mutrate > /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed.mutrate &

rm /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed.split*mutrate

# only keep two relevant columns
awk {'print $1"\t"$4'} /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed.mutrate > \
/media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed.mutrate.temp

mv /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed.mutrate.temp \
/media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed.mutrate

# get phastcons conservation socres for intervals

nohup sh -c 'for i in /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed.split*; do ../lib/bigWigAverageOverBed ../other_annotation/conservation_annotation/hg19.100way.phastCons.bw $i $i.phastcons100way; done' &

nohup cat /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed.split*.phastcons100way \
> /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed.phastcons100way &

# only keeps two relevant columns
awk {'print $1"\t"$4'} /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed.phastcons100way > \
/media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed.phastcons100way.temp

mv /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed.phastcons100way.temp \
/media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed.phastcons100way

rm /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed.split*phastcons100way


# need to assign 1-bp window to genes

# for 10kb non-coding regions

nohup bedtools intersect -a /media/yuwen/Elements/ASD_temp_storage/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.1bp_window.bed \
 -b ../other_annotation/epigenomic_annotation/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed -wa -wb | awk {'print $4"\t"$8'} > \
 /media/yuwen/Elements/ASD_temp_storage/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.1bp_window.bed_genename.bed & 

# for promoter regions

nohup bedtools intersect -a /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_1kb_genename.bed.1bp_window.bed \
-b ../other_annotation/epigenomic_annotation/Whole_genome.promoter_1kb_genename.bed -wa -wb | awk {'print $4"\t"$8'} > \
/media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_1kb_genename.bed.1bp_window.bed_genename.bed & 


# for coding genes
# There is no overlap in /media/yuwen/Elements/ASD_temp_storage/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.1bp_window.bed
# but there are a lot of overlap regions in /media/yuwen/Elements/ASD_temp_storage/other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed

# so /media/yuwen/Elements/ASD_temp_storage/other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.1bp_window.bed_genename.bed has a lot of redundant lines that need to be collapsed. 
# The merge function in R would still work in this scenario.. 

nohup bedtools intersect -a /media/yuwen/Elements/ASD_temp_storage/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.1bp_window.bed \
-b ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed -wa -wb \
| awk {'print $4"\t"$8'} | sort | uniq > /media/yuwen/Elements/ASD_temp_storage/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.1bp_window.bed_genename.bed & 

# add up different parts of genename assignment together
nohup cat /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_1kb_genename.bed.1bp_window.bed_genename.bed \
/media/yuwen/Elements/ASD_temp_storage/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.1bp_window.bed_genename.bed \
/media/yuwen/Elements/ASD_temp_storage/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.1bp_window.bed_genename.bed > \
 /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.1bp_window.bed_genename.bed


### only save two columns for for .phastcons100way file and 