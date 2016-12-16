#use this script to generate merged bed regions for each gene
#! /usr/bin/perl
#perl script.pl region_file(chr \t start \t end \t refseq_ID \t geneID)final_output.file random_prefix
# for each gene, will merge all bed regions that from multiple lines in the region_file.
use POSIX;
use warnings;
use strict;


my $infile1 = shift;
my $final_output = shift; # the final only output file name
my $tmp = shift; # the random prefix


my %data;

#read in region files
open INPUT, $infile1 or die "can't open file\n";

while(<INPUT>){
        chomp;
        my @temp = split;
        my $gene = $temp[4];
        my $enhancer = $temp[3];
        my $chr = $temp[0];
        my $start = $temp[1];
        my $end = $temp[2];
        my $bed = [$chr, $start, $end];
        if(exists $data{$gene}{$enhancer}){
            push @{$data{$gene}{$enhancer}},$bed;  
}
        else{
            $data{$gene}{$enhancer} = [$bed];
}

}
close (INPUT);

#my $output = "_gene_enhancer_pair_temp.txt"; # will contain gene enhancer pair
my $output = $tmp."_gene_enhancer_pair_temp.txt"; # will contain gene enhancer pair

open (OUTPUT, ">$output");

my $output2 = $tmp."_merged_bed_file_temp.txt"; # wil contain merged regions



#print merged regions
foreach my $gene (sort keys %data){
        my $temp_out = $tmp."temp.out.txt";
        open (TEMPOUT, ">$temp_out");
        foreach my $enhancer (sort keys %{$data{$gene}}){                        
            foreach my $array_ref (@{$data{$gene}{$enhancer}}){
                my @array = @{$array_ref};
                    foreach(@array){
                        print TEMPOUT $_."\t";
}
                print TEMPOUT "\n";
}
}
            close(TEMPOUT); 
            my $command1 = "sort -k1,1 -k2,2n ".$tmp."temp.out.txt > ".$tmp."temp.out.txt.sorted";
            system($command1);
            my $command2 = "bedtools merge -i ".$tmp."temp.out.txt.sorted > ".$tmp."temp.out.txt.sorted.merged";
            system($command2);
            my $command3 = "cat ".$tmp."temp.out.txt.sorted.merged >> $output2";
            system($command3);
            my $command4 = "wc -l ".$tmp."temp.out.txt.sorted.merged";
            my $line_count = `$command4`; 
            my @line = split(/\s/, $line_count);
            
            for(my $i = 0; $i < $line[0]; $i++){
                print OUTPUT $gene."\n";
}
}
close(OUTPUT);
#combine $output1 and $output2 column wise
system("paste $output2 $output | column -s \$'\t' -t > $final_output");
system("rm $output2");
system("rm $output");
my $command5 = "rm ".$tmp."temp.out*";
system($command5);