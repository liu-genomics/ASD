library(data.table)
library(parallel)
# This is a function to run glm for 50bp windows in order to adjust for batch-effects on observed mutation rates, then will output a Rmd object that has base level mutation rates 
adjust_mutation_rate_window_to_base <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks = "no", overlap = 0.5, rm_nonsyn = FALSE, annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", annovar_input = "no.txt", 
                                                gene_assign_file, gene_prior_file, report_proportion, node_n =6){
  # [mut_file] is the mutation location file from one dataset. e.g.,../analysis/160229_data_for_analysis.bed, chr \t start \t end \t index,  1-based start and end
  # [window_file] is the file with genome partitioned into windows. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed
  # [cg_file] is the file that has CG-content calculated, matching [window_file]. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.cg
  # [mutrate_file] is the file that has mutrate data. e.g., "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate"
  # [sample_size] is the number of individuals
  # [epitenomic_marks] if it is "no", then glm won't use if a mutation is in h3k27ac as a categorical variable, if it is a epigenomic annotation name.
  # for example, "../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed". It is fine to directly use original epigenomics bed files
  # instead of the partitioned files. Because the window_file already has annotations for each window in terms of whetehr they are promoter or coding. 
  # [overlap] the fraction of overlap (the paramter f in bedtools intersect), default is 0.5, For the ease of TADA modeling, it is better to assign it to be 1E-9(1bp)
  # [rm_nonsyn] TRUE or FALSE. Indicating whether nonysnonymous mutations should be removed when fitting the model. Use TRUE when fitting model using ASD mutations. 
  # [annovar_folder] is where the annovar program is, e.g., . e.g., /media/yuwen/Elements/ANNOVAR/annovar/
  # [annovar_input] is a input file generated from mutations that will be annotated by annovar. The format is chr \t start \t end \t ref_allele \t alt_allele \t index 1-based for starts and ends
  # [gene_assign_file] is the file that assigns each window to a gene. for example ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename.bed
  # [gene_prior_file], a file that has prior (derived from posterior and prior)for a gene as a risk gene. genename \t probability e.g. ../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.94.txt
  # [report_proportion] Choose the top X% TADA genes to report. 
  # [node_n] is the number of nodes used to run , default is 6
  prefix=system("date +%s", intern = TRUE) # prefix for temporary files that will be deleted at the end of the pipeline
  
  mut = read.delim(mut_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  if(rm_nonsyn){
    # command to transform the start base to 0-based in order to use bedtools to do overlap
    command = paste("awk {'print $1\"\t\"$2-1\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$7'} ", annovar_input, " > ", paste(prefix, "_annovar_for_bedtools.bed",sep = ""),sep = "")
    system(command)
    command = paste("grep coding ", window_file, " > ", paste(prefix, "_coding_window_temp.bed", sep = ""), sep = "")
    system(command)
    command = paste("bedtools intersect -a ",paste(prefix, "_annovar_for_bedtools.bed",sep = ""), " -b ", paste(prefix, "_coding_window_temp.bed", sep = ""), " -wa | sort | uniq | awk {'print $1\"\t\"$2+1\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$7'} > ",paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""), sep = "")
    system(command)
    command = paste("perl ", annovar_folder, "annotate_variation.pl -geneanno -buildver hg19 ", paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""), " ", annovar_folder, "humandb/", sep = "")
    system(command)
    coding_anno = read.delim(paste(paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""),".exonic_variant_function", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    removed_id = coding_anno[coding_anno[,2] != "synonymous SNV",10]
    # now remove non-synonymous mutations (every coding mutations that are not synonymous mutations)
    mut = mut[!is.element(mut[,4], removed_id),]
  }
  mut[,2] = mut[,2] - 1  # change 1-based start to 0-based start. 
  write.table(mut, paste(prefix,"_temp.bed", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  command = paste("bedtools coverage -a ", paste(prefix,"_temp.bed", sep = ""), " -b ", window_file, " > ", paste(prefix,"_temp_coverage.bed", sep = ""), sep = "")
  system(command)
  coverage = read.delim(paste(prefix,"_temp_coverage.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  coverage = coverage[,1:5]
  colnames(coverage) = c("chr","start","end","site_index","mut_count")
  coverage$window_size = coverage$end - coverage$start # some regions don't have full length of, say, 50bp
  window = read.delim(window_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  cg = read.delim(cg_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(cg) = c("site_index", "cg")
  coverage = merge(coverage, cg, by.x = "site_index", by.y ="site_index")
  coverage$coding = as.numeric(mapply(function(x,y) grepl(y,x),coverage$site_index, MoreArgs = list(y="coding")))
  coverage$promoter = as.numeric(mapply(function(x,y) grepl(y,x),coverage$site_index, MoreArgs = list(y="promoter")))
  coverage$nf = as.numeric(mapply(function(x,y) grepl(y,x),coverage$site_index, MoreArgs = list(y="nf")))
  mutrate = read.delim(mutrate_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  mutrate = data.frame(site_index = mutrate[,1], mutrate = mutrate[,4])
  coverage = merge(coverage, mutrate, by.x = "site_index", by.y = "site_index")
  coverage = coverage[coverage$mutrate !=0,]
  #fit mutation rate model using a fixed set of features not including histone modification marks. 
  out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
  # get the adjusted mutation rate per base per individual
  coverage = data.frame(coverage, adjusted_mutrate = out.offset$fitted.values/(2*sample_size))
  # assign gene name to windows
  gene_assign = read.delim(gene_assign_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_assign) = c("site_index", "genename")
  coverage = merge(coverage, gene_assign, by.x = "site_index", by.y = "site_index")
  # get the piror probability of genes.
  gene_prior = read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_prior) = c("genename", "prior")
  # change prior to the 
  gene_prior$prior = 1- gene_prior$prior
  
  if(report_proportion !=1){
    genes_for_report = gene_prior[order(gene_prior[,2],decreasing = TRUE),1]
    genes_for_report = genes_for_report[1:floor(length(genes_for_report)*report_proportion)]
    coverage = coverage[is.element(coverage$genename, genes_for_report),]
  }
  
  # get the adjustion effect for each window by dividing adjusted_mutrate by mutrate
  coverage$adjusting_effect = coverage$adjusted_mutrate/coverage$mutrate
  
  # now need to extropolate window mutation file to base level mutation file, now removes coding region, and don't differntiate between promoter and nf
  coverage_noncoding <-coverage[coverage$coding != 1,c("chr","start","end","adjusting_effect","genename")]
  coverage_noncoding$ID <- paste(coverage_noncoding$genename, coverage_noncoding$start, sep = "_")
  
  #funtion to expand windows to base level
  window_expansion <- function(table_row){
    start <- seq(as.numeric(table_row[2]),as.numeric(table_row[3])-1)
    data.frame(table_row[1], start, start+1, paste(table_row[5],start,sep = "_"), table_row[6])
  }
  
  options(warn=-1)
  #options(scipen=999)
  cl <- makeCluster(node_n)
  clusterExport(cl, "window_expansion")
  # use parallel computing and rbindlist to increase the speed by thousands of times. 
  coverage_noncoding_for_base_mutrate <- rbindlist(parApply(cl, coverage_noncoding, 1, window_expansion))
  stopCluster(cl)
  options(warn = 0)
 #options(scipen=0)
  colnames(coverage_noncoding_for_base_mutrate) <- c("chr","start","end","base_ID","ID")
  coverage_noncoding_for_base_mutrate$start <- as.integer(coverage_noncoding_for_base_mutrate$start)
  coverage_noncoding_for_base_mutrate$end <- as.integer(coverage_noncoding_for_base_mutrate$end)
  
  # write out a bed file to get base-level mutation rates
  fwrite(coverage_noncoding_for_base_mutrate[,1:4],paste(prefix, "_temp_for_mutrate.bed", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  command <- paste("../lib/bigWigAverageOverBed ", mutrate_ref_file, " ", paste(prefix, "_temp_for_mutrate.bed", sep = ""), " ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = "" ), sep = "")
  system(command)
  command <- paste("awk {'print $1\"\t\"$4'} ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = ""), " > ", paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), sep = "")
  system(command)
  coverage_noncoding_base_mutrate <-fread(paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  colnames(coverage_noncoding_base_mutrate) <- c("base_ID","base_mutrate")
  a_temp <- coverage_noncoding_base_mutrate[coverage_noncoding_for_base_mutrate, on = "base_ID"]
  a_temp <- a_temp[as.data.table(coverage_noncoding), on = "ID"]
  coverage_noncoding_mutrate_adjusted <- a_temp[,c("base_ID","chr","start","end","genename")]
  coverage_noncoding_mutrate_adjusted$adjusted_base_mutrate <- a_temp[,c("base_mutrate")]*a_temp[,c("adjusting_effect")]
  rm(a_temp)
  
  # get the mutation count for each base
  command = paste("bedtools coverage -a ", paste(prefix,"_temp.bed", sep = ""), " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""), " > ", paste(prefix,"_temp_base_level_coverage.bed", sep = ""), sep = "")
  system(command)
  
  # merge table
  coverage_noncoding_mutrate_adjusted <- coverage_noncoding_mutrate_adjusted[mutation_overlap_base, on = "base_ID"]
  mutation_overlap_base <-fread(paste(prefix,"_temp_base_level_coverage.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  mutation_overlap_base <- mutation_overlap_base[,c("V4","V5")]
  colnames(mutation_overlap_base) <- c("base_ID","mut_count")
  coverage_noncoding_mutrate_adjusted <- coverage_noncoding_mutrate_adjusted[mutation_overlap_base, on = "base_ID"]
  rm(mutation_overlap_base)
  
  # overlap with epi feature.
  command <- paste("bedtools intersect -f ", overlap, " -a ", paste(prefix, "_temp_for_mutrate.bed", sep = "") , " -b ", epigenomic_marks, " -wa > ", paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""),sep = "")
  system(command)
  base_in_epi <- fread(paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  base_in_epi <- data.table(base_ID = base_in_epi[,4], epi = 1)
  colnames(base_in_epi) <- c("base_ID", "epi")
  coverage_noncoding_mutrate_adjusted <- base_in_epi[coverage_noncoding_mutrate_adjusted, on = "base_ID"]
  coverage_noncoding_mutrate_adjusted[is.na(get("epi")), ("epi"):=0]
  fwrite(coverage_noncoding_mutrate_adjusted, "170126_single_base_mutation_table.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  
  system(paste("rm ", prefix, "_temp*", sep = ""))
  
  out.offset
}



new_control = new.env()
load("../data/debug_region_list_073116_using_258_control_data_matrix.Rdata", new_control)

prefix <-system("date +%s", intern = TRUE) # prefix for temporary files that will be deleted at the end of the pipeline

write.table(new_control$data_matrix[new_control$data_matrix$Prediction == "dnv_proband" & new_control$data_matrix$Type == "SNV",c(1,2,3,7)], paste(prefix,"temp.bed",sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(new_control$data_matrix[new_control$data_matrix$Prediction == "dnv_proband" ,c(1:5,8,7)], paste(prefix,"annovar_input_temp.bed",sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)


mut_file = paste(prefix,"temp.bed",sep = "")
window_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed"
cg_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.fasta.cg"
mutrate_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.mutrate"
sample_size = 314
epigenomic_marks = "../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed"
overlap = 0.5
rm_nonsyn = FALSE
annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/"
annovar_input = "no.txt"
gene_assign_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename.bed"
gene_prior_file = "../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.94.txt"
report_proportion = 0.06
mutrate_ref_file = "/media/yuwen/Elements/mutation_rate_raw/MarkDaly/Daly_mutrate.bw"
node_n = 6
