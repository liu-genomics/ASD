#!/bin/bash
# This script have enhancer bed files as INPUT (disjoin, nonoverlapped regions), and output mutations that are in the enhancers that are within certain kilobaese of TSSs
# The final output would be a tow-column file, with the 1st column as mutation ID and the second column as associated gene names.


# $1 is the name of the enhancer region bed file (chr \t start \t end \t enhancer ID)
# $2 is the folder with refseq annotations (promoter, 5'UTR, 3'UTR and etc.), default is '../other_annotation/refseq/'
# $3 is the output folder, default is '../other_annotation/epigenomic_annotation/'
# $4 is the output prefix, e.g., "Noonan_brain"

# The output of this script should be a set of bed files with epigenomics marks covered promoters, 5'utrs and 3'utrs.

# run promoters
prefix=`date +%s` # prefix for temporary files that will be deleted at the end of the pipeline

bedtools intersect -a $1 -b $2/refseq_id_promoter_1kb_upstream_pair_table.txt -wb> $prefix.tmp
Rscript ../lib/link_refid_genename_v2.R $2/ref_id_to_genename_table.txt $prefix.tmp $prefix.promoter.tmp
awk {'print $1"\t"$2"\t"$3"\t"$5"\t"$4'} $prefix.promoter.tmp | grep -v "n/a" > $prefix.promoter.1.tmp

perl ../lib/160920_merge_regions_for_each_gene.pl $prefix.promoter.1.tmp $prefix.promoter.2.tmp $prefix

awk {'print $1"\t"$2"\t"$3"\t"$4'} $prefix.promoter.2.tmp > $3/$4.promoter_1kb_genename.bed

rm $prefix*tmp

# run 5-utrs
bedtools intersect -a $1 -b $2/refseq_id_5utr_pair_table.txt -wb> $prefix.tmp
Rscript ../lib/link_refid_genename_v2.R $2/ref_id_to_genename_table.txt $prefix.tmp $prefix.utr5.tmp
awk {'print $1"\t"$2"\t"$3"\t"$5"\t"$4'} $prefix.utr5.tmp | grep -v "n/a" > $prefix.utr5.1.tmp

perl ../lib/160920_merge_regions_for_each_gene.pl $prefix.utr5.1.tmp $prefix.utr5.2.tmp $prefix

awk {'print $1"\t"$2"\t"$3"\t"$4'} $prefix.utr5.2.tmp > $3/$4.5utr_genename.bed

rm $prefix*tmp

#run 3-utrs
bedtools intersect -a $1 -b $2/refseq_id_3utr_pair_table.txt -wb> $prefix.tmp
Rscript ../lib/link_refid_genename_v2.R $2/ref_id_to_genename_table.txt $prefix.tmp $prefix.utr3.tmp
awk {'print $1"\t"$2"\t"$3"\t"$5"\t"$4'} $prefix.utr3.tmp | grep -v "n/a" > $prefix.utr3.1.tmp

perl ../lib/160920_merge_regions_for_each_gene.pl $prefix.utr3.1.tmp $prefix.utr3.2.tmp $prefix

awk {'print $1"\t"$2"\t"$3"\t"$4'} $prefix.utr3.2.tmp > $3/$4.3utr_genename.bed

rm $prefix*tmp

