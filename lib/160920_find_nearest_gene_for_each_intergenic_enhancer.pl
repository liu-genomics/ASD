
#first build a hash for to store start position of promoter of each gene, then for each intergenic enhancer find the nearest TSS.
#perl script.pl union_for_k27ac_10kb_for_each_gene_enhancer_pair_not_overlapped_genebody_not_promoter_utr_exon.txt folder_of_tss_posiion.txt target_file_nearest_gene.txt
use POSIX;
use warnings;
use strict;
use List::Util qw( min max );


my $infile1 = shift; #k27ac covered region 
my $infile2 = shift; # the folder of tss_position.txt
$infile2 .= "/tss_position.txt"; # refid to genename transition table



open INPUT, $infile2 or die "can't open file\n";
my %data;
while(<INPUT>){
	chomp;
	my @temp = split;
	my $gene = $temp[1];
        my $refseq = $temp[0];
        my $tss = $temp[3];
        $data{$gene}{$refseq} = $tss;
}
close (INPUT);

open INPUT, $infile1 or die "can't open file\n";
my %enhancer_position;
while(<INPUT>){
        chomp;
	my @temp = split;
	my $enhancer = $temp[4];
        my $gene = $temp[3];
        my $start = $temp[1];
        my $end = $temp[2];
        my $middle = ($end - $start)/2+$start;
        $enhancer_position{$enhancer}{$gene} = $middle;       
}
close (INPUT);

open INPUT, $infile1 or die "can't open file\n";
my %enhancer_info;
while(<INPUT>){
        chomp;
	my @temp = split;
	my $enhancer = $temp[4];
        my $gene = $temp[3];
        my $start = $temp[1];
        my $end = $temp[2];
	my $chr = $temp[0];
        $enhancer_info{$enhancer}{"chr"} = $chr;
	$enhancer_info{$enhancer}{"start"} = $start;
	$enhancer_info{$enhancer}{"end"} = $end;
	$enhancer_info{$enhancer}{"name"} = $enhancer       
}
close (INPUT);

my $infile3 = shift;
open (OUTPUT, ">$infile3") or die "can't open file\n";

#for each intergenic enhancer assign nearest gene
foreach my $enhancer (sort keys %enhancer_position){
	my %gene_distance; # a hash where key is the gene, value is the minimal distance to that gene (gene may have multiple tss) 
	foreach my $gene (sort keys %{$enhancer_position{$enhancer}}){
		if(exists $data{$gene}){		
			my @distance;
			foreach my $refseq (sort keys %{$data{$gene}}){
				my $distance_refseq = abs($enhancer_position{$enhancer}{$gene} - $data{$gene}{$refseq});
				#print OUTPUT "\n"."check".$refseq."\t".$distance_refseq."\n";
				push @distance, $distance_refseq;
}
			my $min_distance_this_gene = min @distance;
			$gene_distance{$gene} = $min_distance_this_gene;
}
		else{
			$gene_distance{$gene} = 100000000 # assign an arbitory large value, so it won't be chosen unless it is the only one gene within certain prevously chosen cutoff distance
}
}
	my @temp =  sort {$gene_distance{$a} <=> $gene_distance{$b}} keys %gene_distance;
	my $nearest_gene = $temp[0];
	
	print OUTPUT $enhancer_info{$enhancer}{"chr"}."\t".$enhancer_info{$enhancer}{"start"}."\t".$enhancer_info{$enhancer}{"end"}."\t".$nearest_gene."\t".$enhancer_info{$enhancer}{"name"}."\n";
} 

close(OUTPUT);



