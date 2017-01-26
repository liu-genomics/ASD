# this code is used to transform yanyu's data structure to mine so that previous burden analysis code could be used

args <- commandArgs(trailingOnly = TRUE)
print(args)
Rdata <- args[1] # prefix of Rdata that needs to be transformed, e.g., 0703_region_list_080216_data_matrix
mutation_type <- args[2] #SNVs or indels
motif_pct_cutoff<- as.numeric(args[3]) # for example 0.9
annovar_folder <- args[4] # for example /media/yuwen/Elements/ANNOVAR/annovar/
spidex_folder <- args[5] # for example /media/yuwen/Elements/spidex_database/, this folder at least has spidex_public_noncommercial_v1_0.tab.gz and spedix_header

#regions_file <- args[2]
#region_out <- basename(regions_file)
#motifs <- args[3]
#rm(args)

# read in data
load(paste(Rdata,".Rdata", sep = ""))
mutation = data.frame(chr = data_matrix$Chrom, start = data_matrix$Start, end = data_matrix$End, index = data_matrix$ID, mut_type = data_matrix$Type, 
                      study = "SCC", 
                      phenotype = data_matrix$Prediction)
mutation$phenotype = as.character(mutation$phenotype)
mutation[mutation$phenotype == "dnv_proband",]$phenotype = "ASD"
mutation[mutation$phenotype == "dnv_sibling",]$phenotype = "control"
# set mutation ID labels for cases and controls

ASD_effective_SNV_ID = mutation[mutation$phenotype == "ASD" & mutation$mut_type == "SNV",]$index
control_effective_SNV_ID = mutation[mutation$phenotype == "control" & mutation$mut_type == "SNV",]$index

# set filtering flags, the Rdata that has been read in
unique_dummy = rep(1, length(mutation[,1]))
individual_dummy_140 = rep(1, length(mutation[,1]))
repeat_dummy = rep(1, length(mutation[,1]))

#unique mutation information
unique_mut_ID = seq(1:length(mutation$mutation[,1]))
unique_mut_ID = unique_mut_ID[mutation$unique_dummy ==1]

#flagged patients information  , remove patients with more than 140 de novo mutations

individual_mut_ID = seq(1:length(mutation$mutation[,1]))
individual_mut_ID = individual_mut_ID[mutation$individual_dummy_140 ==1]

#flag indels that are in repeat regions

nonrepeat_mut_indel_ID = seq(1:length(mutation$mutation[,1]))
nonrepeat_mut_indel_ID = nonrepeat_mut_indel_ID[mutation$repeat_dummy != 1 & mutation$mutation$mut_type != "SNV"]


# set annotations for mutations
CADD_score = data.frame(data_matrix$CADD13_PHRED)
gerp = data.frame(NA,NA,NA, data_matrix$ID, NA, NA, NA, data_matrix$GERP)
mut_eigen = data_matrix$Eigen
phylop = data.frame(data_matrix$ID, NA, NA, NA, NA, data_matrix$phyloP46wayAllElements)
motif_cutoff = quantile(data_matrix$Noonan_brain_motif_top50pct_expressed_motif.txt.ref_alt.3.m.pvalue, c(motif_pct_cutoff))
jaspar_2014_motif_q0.1_dummy = as.numeric(data_matrix$Noonan_brain_motif_top50pct_expressed_motif.txt.ref_alt.3.m.pvalue >= motif_cutoff)

# there are some annotations that haven't been added in Yanyu's data
fitcons = data.frame(data_matrix$ID, NA, NA, NA, NA, NA)

# add a data.frame to have alt_ref_allele information
ref_alt_allele = data_matrix[,4:5]

# use Annovar to annotate coding mutations. 
prefix=system("date +%s", intern = TRUE) # prefix for temporary files that will be deleted at the end of the pipeline
temp_for_annovar = data.frame(mutation[,1:3], ref_alt_allele, mutation$mut_type, mutation$index)
write.table(temp_for_annovar, paste(prefix, "_for_annovar_temp.txt",sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

command = paste("perl ", annovar_folder, "annotate_variation.pl -geneanno -buildver hg19 ", paste(prefix, "_for_annovar_temp.txt",sep = ""), " ", annovar_folder, "humandb/", sep = "")
system(command)
coding_anno = read.delim(paste(paste(prefix, "_for_annovar_temp.txt",sep = ""),".exonic_variant_function", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
coding_gene = unlist(lapply(strsplit(coding_anno[,3],split = ":"),`[[`,1))
exon_mut_type = rep("unknown", length(mutation[,1]))
exon_mut_type[is.element(mutation$index, coding_anno[coding_anno[,2] == "frameshift deletion",10])] = "frameshift deletion"
exon_mut_type[is.element(mutation$index, coding_anno[coding_anno[,2] == "frameshift insertion",10])] = "frameshift insertion"
exon_mut_type[is.element(mutation$index, coding_anno[coding_anno[,2] == "nonframeshift deletion",10])] = "nonframeshift deletion"
exon_mut_type[is.element(mutation$index, coding_anno[coding_anno[,2] == "nonsynonymous SNV",10])] = "nonsynonymous"
exon_mut_type[is.element(mutation$index, coding_anno[coding_anno[,2] == "stopgain",10])] = "stopgain"
exon_mut_type[is.element(mutation$index, coding_anno[coding_anno[,2] == "synonymous SNV",10])] = "synonymous"
system(paste("rm ", prefix, "_for_annovar_temp.txt*",sep = ""))

# here th synonymous mutation number in cases and in controls are totally determined by annotation from ANNOVAR. Previously, I only considered mutations in genes that are in the 
# list of genes with muation rate. 

syn_number_ASD = length(mutation[exon_mut_type == "synonymous" & mutation$phenotype == "ASD",1])
syn_number_control = length(mutation[exon_mut_type == "synonymous" & mutation$phenotype == "control",1])

coding_mut_to_gene = data.frame(coding_anno[,10], coding_gene)


#calculate CG content for 50bp, 100bp 200bp and 500bp window, 
CG_content = data.frame(index = mutation$index, CG_50bp = numeric(nrow(mutation)), CG_100bp = numeric(nrow(mutation)), CG_200bp = numeric(nrow(mutation)), CG_500bp = numeric(nrow(mutation)))
calculate_cg <-function(seq){
  letter = as.data.frame(strsplit(c(seq),split=""))
  GC_count = length(letter[letter[,1] == "C"|letter[,1] == "c" | letter[,1] == "G" | letter[,1] == "g",1])
  GC_pct = GC_count/length(letter[,1]) 
  GC_pct
}
prefix=system("date +%s", intern = TRUE) # prefix for temporary files that will be deleted at the end of the pipeline

write.table(mutation[,1:4], paste(prefix, "mutation_all.bed",sep = "_"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

command = paste("awk {'print $1\"\t\"$2-26\"\t\"$3+25\"\t\"$4'} ", paste(prefix, "mutation_all.bed",sep = "_")," > ", paste(prefix,"mutation_all_50bp_around.bed",sep = ""),sep = "")
system(command)
command = paste("bedtools getfasta -tab -name -fi ../other_annotation/genome_build/hg19.fasta -bed ",paste(prefix,"mutation_all_50bp_around.bed",sep = "")," -fo ", paste(prefix,"mutation_all_50bp_around.bed.fasta",sep = ""),sep = "")
system(command)
temp = read.delim(paste(prefix,"mutation_all_50bp_around.bed.fasta",sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
cg_pct = data.frame(temp[,1],mapply(calculate_cg,temp[,2]))
CG_content$CG_50bp = cg_pct[,2]

command = paste("awk {'print $1\"\t\"$2-51\"\t\"$3+50\"\t\"$4'} ", paste(prefix, "mutation_all.bed",sep = "_")," > ", paste(prefix,"mutation_all_100bp_around.bed",sep = ""),sep = "")
system(command)
command = paste("bedtools getfasta -tab -name -fi ../other_annotation/genome_build/hg19.fasta -bed ",paste(prefix,"mutation_all_100bp_around.bed",sep = "")," -fo ", paste(prefix,"mutation_all_100bp_around.bed.fasta",sep = ""),sep = "")
system(command)
temp = read.delim(paste(prefix,"mutation_all_100bp_around.bed.fasta",sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
cg_pct = data.frame(temp[,1],mapply(calculate_cg,temp[,2]))
CG_content$CG_100bp = cg_pct[,2]

command = paste("awk {'print $1\"\t\"$2-101\"\t\"$3+100\"\t\"$4'} ", paste(prefix, "mutation_all.bed",sep = "_")," > ", paste(prefix,"mutation_all_200bp_around.bed",sep = ""),sep = "")
system(command)
command = paste("bedtools getfasta -tab -name -fi ../other_annotation/genome_build/hg19.fasta -bed ",paste(prefix,"mutation_all_200bp_around.bed",sep = "")," -fo ", paste(prefix,"mutation_all_200bp_around.bed.fasta",sep = ""),sep = "")
system(command)
temp = read.delim(paste(prefix,"mutation_all_200bp_around.bed.fasta",sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
cg_pct = data.frame(temp[,1],mapply(calculate_cg,temp[,2]))
CG_content$CG_200bp = cg_pct[,2]

command = paste("awk {'print $1\"\t\"$2-251\"\t\"$3+250\"\t\"$4'} ", paste(prefix, "mutation_all.bed",sep = "_")," > ", paste(prefix,"mutation_all_500bp_around.bed",sep = ""),sep = "")
system(command)
command = paste("bedtools getfasta -tab -name -fi ../other_annotation/genome_build/hg19.fasta -bed ",paste(prefix,"mutation_all_500bp_around.bed",sep = "")," -fo ", paste(prefix,"mutation_all_500bp_around.bed.fasta",sep = ""),sep = "")
system(command)
temp = read.delim(paste(prefix,"mutation_all_500bp_around.bed.fasta",sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
cg_pct = data.frame(temp[,1],mapply(calculate_cg,temp[,2]))
CG_content$CG_500bp = cg_pct[,2]
CG_content$index = cg_pct[,1]

system(paste("rm ", prefix, "*mutation*bed*",sep = ""))

### get the spidex output
temp_for_spedix = data.frame(mutation[,1:2], mutation$index, ref_alt_allele)
write.table(temp_for_spedix, paste(prefix, "_for_spedix_temp.txt",sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
command <- paste("python ../lib/yanyu_CRI_query_tabix_for_python_3.py ", paste(prefix, "_for_spedix_temp.txt ",sep = ""), paste(spidex_folder, "spidex_public_noncommercial_v1_0.tab.gz ",sep = "/"), paste(spidex_folder, "spidex_header", sep = "/"), " > ", paste(prefix,"_for_spedix_temp.output",sep = ""), sep = "")
system(command)
mut_spedix <- read.delim(paste(prefix,"_for_spedix_temp.output",sep = ""), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
system(paste("rm ", prefix, "_for_spidex_temp.txt",sep = ""))
system(paste("rm ", prefix, "_for_spidex_temp.output",sep = ""))

### save to output
output = paste(Rdata, "transformed_for_old_code.Rdata", sep = "")
save.image(output)
q(save = "no")


