#Use predefined regions (within 10kb, removing utr, coding sequence, and promoter), calculate CG content for each. 

# nf sequences, within 10kb, non-promoter, non-coding and non-ute
awk {'print $1"\t"$2"\t"$3'} ../other_annotation/epigenomic_annotation/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed | sort -k1,1 -k2,2n > ../other_annotation/epigenomic_annotation/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.temp1

bedtools merge -i ../other_annotation/epigenomic_annotation/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.temp1 > ../other_annotation/epigenomic_annotation/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.temp2

bedtools makewindows -b ../other_annotation/epigenomic_annotation/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.temp2 -w 1 |awk {'print $1"\t"$2"\t"$3"\tnf"NR'} > ../other_annotation/epigenomic_annotation/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.1bp_window.bed

# move files to a bigger repository

mkdir -p /media/yuwen/Elements/ASD_temp_storage/

mv ../other_annotation/epigenomic_annotation/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.1bp_window.bed /media/yuwen/Elements/ASD_temp_storage/

nohup awk {'print $1"\t"$2-25"\t"$3+25"\t"$4'} /media/yuwen/Elements/ASD_temp_storage/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.1bp_window.bed > \ 
/media/yuwen/Elements/ASD_temp_storage/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.1bp_window.bed_50bp_surrounding.bed &


nohup bedtools getfasta -tab -name -fi ../other_annotation/genome_build/hg19.fasta -bed /media/yuwen/Elements/ASD_temp_storage/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.1bp_window.bed_50bp_surrounding.bed -fo /media/yuwen/Elements/ASD_temp_storage/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.1bp_window.bed_50bp_surrounding.bed.fasta &

nohup split -d -l 10000000 /media/yuwen/Elements/ASD_temp_storage/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.1bp_window.bed_50bp_surrounding.bed.fasta Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.1bp_window.bed_50bp_surrounding.bed.fasta.split &

nohup sh -c 'for i in /media/yuwen/Elements/ASD_temp_storage/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.1bp_window.bed_50bp_surrounding.bed.fasta.split*; do Rscript ../lib/161230_calculate_cg_parallele_version.R $i 6; done' &


cat /media/yuwen/Elements/ASD_temp_storage/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.1bp_window.bed_50bp_surrounding.bed.fasta.split*cg \
> /media/yuwen/Elements/ASD_temp_storage/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.1bp_window.bed_50bp_surrounding.bed.fasta.cg

rm /media/yuwen/Elements/ASD_temp_storage/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.1bp_window.bed_50bp_surrounding.bed.fasta.split*


rm ../other_annotation/epigenomic_annotation/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.temp*

# promoter

awk {'print $1"\t"$2"\t"$3'} ../other_annotation/epigenomic_annotation/Whole_genome.promoter_1kb_genename.bed | sort -k1,1 -k2,2n > ../other_annotation/epigenomic_annotation/Whole_genome.promoter_1kb_genename.bed.temp1

bedtools merge -i ../other_annotation/epigenomic_annotation/Whole_genome.promoter_1kb_genename.bed.temp1 > ../other_annotation/epigenomic_annotation/Whole_genome.promoter_1kb_genename.bed.temp2

bedtools makewindows -b ../other_annotation/epigenomic_annotation/Whole_genome.promoter_1kb_genename.bed.temp2 -w 1 |awk {'print $1"\t"$2"\t"$3"\tpromoter"NR'} > ../other_annotation/epigenomic_annotation/Whole_genome.promoter_1kb_genename.bed.1bp_window.bed

mv ../other_annotation/epigenomic_annotation/Whole_genome.promoter_1kb_genename.bed.1bp_window.bed /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_1kb_genename.bed.1bp_window.bed

nohup awk {'print $1"\t"$2-25"\t"$3+25"\t"$4'} /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_1kb_genename.bed.1bp_window.bed > /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_1kb_genename.bed.1bp_window.bed_50bp_surrounding.bed & 



nohup bedtools getfasta -tab -name -fi ../other_annotation/genome_build/hg19.fasta -bed /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_1kb_genename.bed.1bp_window.bed_50bp_surrounding.bed -fo /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_1kb_genename.bed.1bp_window.bed_50bp_surrounding.bed.fasta &

nohup Rscript ../lib/161230_calculate_cg_parallele_version.R /media/yuwen/Elements/ASD_temp_storage/Whole_genome.promoter_1kb_genename.bed.1bp_window.bed_50bp_surrounding.bed.fasta 5 &

rm ../other_annotation/epigenomic_annotation/Whole_genome.promoter_1kb_genename.bed.temp*


# then coding regions. 


# The aim is to prepare data for using GLM Pisson model to calibrate mutation rates
# The coding region now is all ref-seq exons (including those with a gene name and those without, subtracted any part that overlaps with promoters)
bedtools subtract -a ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes.bed -b ../other_annotation/epigenomic_annotation/Whole_genome.promoter_1kb_genename.bed > ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb.bed

sort -k1,1 -k2,2n ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb.bed > ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_temp1.bed

bedtools merge -i ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_temp1.bed > ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_temp2.bed

#also need to remove sequences from unassembled chromatin structure
bedtools makewindows -w 50 -b ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_temp2.bed | grep -v _ |awk {'print $1"\t"$2"\t"$3"\tcoding"NR'} > ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb.bed.50bp_window.bed

bedtools getfasta -tab -name -bed ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb.bed.50bp_window.bed -fi ../other_annotation/genome_build/hg19.fasta -fo ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb.bed.50bp_window.bed.fasta

nohup Rscript ../lib/161117_calculate_cg.R ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb.bed.50bp_window.bed.fasta &
