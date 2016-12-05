#!/bin/bash
# This script have enhancer bed files as INPUT (disjoin, nonoverlapped regions), and output mutations that are in the enhancers that are within certain kilobaese of TSSs
# The final output would be a tow-column file, with the 1st column as mutation ID and the second column as associated gene names.

# specifically used for whole_genome.bed

# $1 is the name of the enhancer region bed file (chr \t start \t end \t enhancer ID)
# $2 is the maximal distance that is allowed to assign an enhancer to a TSS. 
# $3 is the folder with refseq annotations (promoter, 5'UTR, 3'UTR and etc.), default is '../other_annotation/refseq/'
# $4 is the output folder, default is '../other_annotation/epigenomic_annotation/'
# $5 is the output prefix, e.g., "Noonan_brain"

# The output of this script should be a set of files (epigenomic partitions), I will prioritize genomic annotations in the following rank
# promoter （will allow multiple gene assignment for now, as it is generally short）
# 5'UTR (will allow multiple gene assignment for now, as it is generally short)
# 3'UTR (will allow multiple gene assignment for now, as it is generally short）
# enhancers (will be assigned to nearest genes)

#remove other genomic features from enhancers, need to include introns at this point
prefix=`date +%s` # prefix for temporary files that will be deleted at the end of the pipeline
bedtools subtract -a $1 -b $3/refseq_id_promoter_1kb_upstream_pair_table.txt > $prefix.tmp
bedtools subtract -a $prefix.tmp -b $3/refseq_id_5utr_pair_table.txt > $prefix.2.tmp
bedtools subtract -a $prefix.2.tmp -b $3/refseq_id_3utr_pair_table.txt > $prefix.3.tmp
bedtools subtract -a $prefix.3.tmp -b $3/refseq_id_exons_genomic_coordiantes.bed > $prefix.4.tmp
#bedtools subtract -a temp2.txt -b /media/yuwen/F/tools/refseq_gene_coversions/table_files/refseq_id_intron_pair_table.txt > $1.temp

#fix the bug that previous enhancer identifier would not be unique anymore, for example E12345, with an exon in the middle would become two seprarated enhancers with a same
#enhancer name. So I need to give each separated enhancer a new name, or some enhancers would be missed
                                                                                                                                                            
awk {'print $1"\t"$2"\t"$3"\t"NR'} $prefix.4.tmp > $prefix.5.tmp

#remove genes not assembled to 23 chromosomes and transform to BED format
grep -vP "^.*\t.*\t.*_.*\t" $3/tss_position.txt | awk {'print $3"\t"$4"\t"$4"\t"$2'} > $prefix.6.tmp
#adjust distance
bedtools flank -i $prefix.6.tmp -g ../other_annotation/genome_build/hg19.genome -l $2 -r $2 > $prefix.7.tmp
grep -v "n/a" $prefix.7.tmp > $prefix.8.tmp

# the code below is different from 160920_partition_enhancers_into_different_genomic_features.sh
# only genomic regions within distance cut-offs will be retained. But for enhancers, anything that has overlap with extended boundaries will be retained
# here it is possible that one mutation will be assigned to multiple genes. 
bedtools intersect -a $prefix.8.tmp -b $prefix.5.tmp > $prefix.9.tmp

# will not use ../lib/160920_find_nearest_gene_for_each_intergenic_enhancer.pl to assign nearest genes, so a whole-genome mutation could be assigned to multiple genes

awk {'print $1"\t"$2"\t"$3"\t"$4'} $prefix.9.tmp > $4/$5.enhancers.$2.bp_within_TSS.bed

rm $prefix.*.tmp



