
# What is different here from Pipeline 14 is that for coding I will remove regions that overlap with 5'UTR or 3'UTR and those
# that don't have a refseq gene name (For promoters I always only kept regions with a refseq genename). 

# I also included intron here, as later on I need to estimate the splicing mutation rates.

#Use predefined regions (within 10kb, removing utr, coding sequence, and promoter), calculate CG content for each. 
 
# Need to use some of the files generated in Pipeline 14. 

# Now need to get intron coordinates


# The intron region now is all ref-seq introns (including those with a gene name and those without, need to subtract any part that overlaps with promoters)
bedtools subtract -a ../other_annotation/refseq/refseq_id_introns_genomic_coordinates.bed -b ../other_annotation/epigenomic_annotation/Whole_genome.promoter_1kb_genename.bed > ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb.bed

sort -k1,1 -k2,2n ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb.bed > ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_temp1.bed

bedtools merge -i ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_temp1.bed > ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_temp2.bed

#also need to remove sequences from unassembled chromatin structure
bedtools makewindows -w 50 -b ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_temp2.bed | grep -v _ |awk {'print $1"\t"$2"\t"$3"\tintron"NR'} > ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb.bed.50bp_window.bed

bedtools getfasta -tab -name -bed ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb.bed.50bp_window.bed -fi ../other_annotation/genome_build/hg19.fasta -fo ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb.bed.50bp_window.bed.fasta

nohup Rscript ../lib/161117_calculate_cg.R ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb.bed.50bp_window.bed.fasta &



# The aim is to prepare data for using GLM Pisson model to calibrate mutation rates
# The coding region now is all ref-seq exons (including those with a gene name and those without, subtracted any part that overlaps with promoters)
# now need to further remove 3UTR and 5UTR for introns, need to remove exons, need to remove non-coding regions that have already been included in "nf", otherwise, regions for each gene won't be disjoint

bedtools subtract -a ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb.bed -b ../other_annotation/refseq/refseq_id_3_UTR_list.txt > ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_3UTR.bed

bedtools subtract -a ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_3UTR.bed -b ../other_annotation/refseq/refseq_id_5utr_pair_table.txt > ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_3UTR_5UTR.bed

bedtools subtract -a ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_3UTR_5UTR.bed \
-b ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed \
> ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_3UTR_5UTR_exon.bed


bedtools subtract -a ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_3UTR_5UTR_exon.bed \
-b ../other_annotation/epigenomic_annotation/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed \
> ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_3UTR_5UTR_exon_nf.bed  # could be thought as intron specific

Rscript ../lib/link_refid_genename_for_regular_bed.R ../other_annotation/refseq/ref_id_to_genename_table.txt \ 
../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_3UTR_5UTR_exon_nf.bed  \
../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_3UTR_5UTR_exon_nf_removing_na_gene.bed

sort -k1,1 -k2,2n ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_3UTR_5UTR_exon_nf_removing_na_gene.bed > \
../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_3UTR_5UTR_exon_nf_removing_na_gene_temp1.bed

bedtools merge -i ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_3UTR_5UTR_exon_nf_removing_na_gene_temp1.bed > \
../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_3UTR_5UTR_exon_nf_removing_na_gene_temp2.bed




#also need to remove sequences from unassembled chromatin structure
bedtools makewindows -w 50 -b ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_3UTR_5UTR_exon_nf_removing_na_gene_temp2.bed | grep -v _ |awk {'print $1"\t"$2"\t"$3"\tintron"NR'} > ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_3UTR_5UTR_exon_nf_removing_na_gene.bed.50bp_window.bed

bedtools getfasta -tab -name -bed \
../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_3UTR_5UTR_exon_nf_removing_na_gene.bed.50bp_window.bed \
-fi ../other_annotation/genome_build/hg19.fasta \
-fo ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_3UTR_5UTR_exon_nf_removing_na_gene.bed.50bp_window.bed.fasta

nohup Rscript ../lib/161117_calculate_cg.R ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_3UTR_5UTR_exon_nf_removing_na_gene.bed.50bp_window.bed.fasta &

# implement up to here

nohup cat ../other_annotation/epigenomic_annotation/Whole_genome.promoter_1kb_genename.bed.50bp_window.bed \
../other_annotation/epigenomic_annotation/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.50bp_window.bed \
../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.50bp_window.bed \
../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_3UTR_5UTR_exon_nf_removing_na_gene.bed.50bp_window.bed \
> ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_introns_no_utr_no_na_genes.50bp_window.bed &

#cat ../other_annotation/epigenomic_annotation/Whole_genome.promoter_1kb_genename.bed.50bp_window.bed.fasta \
#../other_annotation/epigenomic_annotation/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.50bp_window.bed.fasta \
#../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.50bp_window.bed.fasta \
#> ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.fasta

nohup cat ../other_annotation/epigenomic_annotation/Whole_genome.promoter_1kb_genename.bed.50bp_window.bed.fasta.cg \
../other_annotation/epigenomic_annotation/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.50bp_window.bed.fasta.cg \
../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.50bp_window.bed.fasta.cg \
../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_3UTR_5UTR_exon_nf_removing_na_gene.bed.50bp_window.bed.fasta.cg \
> ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_introns_no_utr_no_na_genes.50bp_window.bed.fasta.cg &



#### change the strategy to only include regions that have defined splicing scores, in order to reduce computation burdens. 5UTR and 3UTR were included if there are defined spedix scores
bedtools subtract -a ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb.bed \
-b ../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed \
> ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_exon.bed

bedtools subtract -a ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_exon.bed \
-b ../other_annotation/epigenomic_annotation/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed \
> ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_exon_nf.bed

# only retain the sequences that have precomputed spidex scores
bedtools intersect -a ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_exon_nf.bed \
-b ../other_annotation/spidex_database/spidex_public_noncommercial_v1_0_collapse_unique_base_merged.bed \
> ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_exon_nf_with_spidex.bed


nohup Rscript ../lib/link_refid_genename_for_regular_bed.R ../other_annotation/refseq/ref_id_to_genename_table.txt \ 
../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_exon_nf_with_spidex.bed  \
../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_exon_nf_with_spidex_removing_na_gene.bed &

sort -k1,1 -k2,2n ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_exon_nf_with_spidex_removing_na_gene.bed > \
../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_exon_nf_with_spidex_removing_na_gene_temp1.bed

bedtools merge -i ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_exon_nf_with_spidex_removing_na_gene_temp1.bed > \
../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_exon_nf_with_spidex_removing_na_gene_temp2.bed

#also need to remove sequences from unassembled chromatin structure
bedtools makewindows -w 50 -b ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_exon_nf_with_spidex_removing_na_gene_temp2.bed | grep -v _ |awk {'print $1"\t"$2"\t"$3"\tspdx_intron"NR'} \
> ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_exon_nf_with_spidex_removing_na_gene.bed.50bp_window.bed

bedtools getfasta -tab -name -bed \
../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_exon_nf_with_spidex_removing_na_gene.bed.50bp_window.bed \
-fi ../other_annotation/genome_build/hg19.fasta \
-fo ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_exon_nf_with_spidex_removing_na_gene.bed.50bp_window.bed.fasta

nohup Rscript ../lib/161117_calculate_cg.R ../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_exon_nf_with_spidex_removing_na_gene.bed.50bp_window.bed.fasta &

#combine previous files with the intron files, intron part may have utrs that have spedix scores

cat ../other_annotation/epigenomic_annotation/Whole_genome.promoter_1kb_genename.bed.50bp_window.bed \
../other_annotation/epigenomic_annotation/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.50bp_window.bed \
../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.50bp_window.bed \
../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_exon_nf_with_spidex_removing_na_gene.bed.50bp_window.bed \
> ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_introns_with_spidex_no_na_genes.50bp_window.bed

cat ../other_annotation/epigenomic_annotation/Whole_genome.promoter_1kb_genename.bed.50bp_window.bed.fasta.cg \
../other_annotation/epigenomic_annotation/Whole_genome.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed.50bp_window.bed.fasta.cg \
../other_annotation/refseq/refseq_id_exons_genomic_coordiantes_subtract_promoter_1kb_3UTR_5UTR_removing_na_gene.bed.50bp_window.bed.fasta.cg \
../other_annotation/refseq/refseq_id_introns_genomic_coordinates_subtract_promoter_1kb_exon_nf_with_spidex_removing_na_gene.bed.50bp_window.bed.fasta.cg \
> ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_introns_with_spidex_no_na_genes.50bp_window.bed.fasta.cg


nohup ../lib/bigWigAverageOverBed ../other_annotation/mutation_rate/Daly_mutrate.bw \
../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_introns_with_spidex_no_na_genes.50bp_window.bed \
../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_introns_with_spidex_no_na_genes.50bp_window.bed.mutrate &

# get phastcons conservation socres for intervals
../lib/bigWigAverageOverBed ../other_annotation/conservation_annotation/hg19.100way.phastCons.bw \
../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed \
../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.phastcons100way



