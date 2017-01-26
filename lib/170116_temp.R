mut_file = paste(prefix,"temp.bed",sep = "")
window_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed"
cg_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.fasta.cg"
mutrate_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.mutrate"
sample_size = 314
epigenomic_marks_list = c("../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed")
overlap = 1e-9
gene_assign_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename.bed"
rm_nonsyn = FALSE
annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/"
annovar_input = paste(prefix,"annovar_input_temp.bed",sep = "")
gene_prior_file = "../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.94.txt"
sequence_annotation_list = "no"
mut_spidex = proband_SNV_mut_spedix
rr = c(3,2.1)
coding_bayes_table_file = "../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_BF_for_de_novo_mutations.txt" 
TADA_p0 = 0.94



# The Version 2 of the function used all splicing mutations, including those mutations that are not overlapped by the window_file.
bayes_factor_for_each_gene_v2 <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks_list = "no", overlap = 0.5, gene_assign_file,rm_nonsyn = FALSE, annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", annovar_input = "no.txt", gene_prior_file, 
                                       sequence_annotation_list = "no", mut_spidex, rr, coding_bayes_table_file, TADA_p0 = 0.94){
  # [mut_file] is the mutation location file from one dataset. e.g.,../analysis/160229_data_for_analysis.bed, chr \t start \t end \t index,  1-based start and end
  # [window_file] is the file with genome partitioned into windows. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed
  # [cg_file] is the file that has CG-content calculated, matching [window_file]. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.cg
  # [mutrate_file] is the file that has mutrate data. e.g., "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate"
  # [sample_size] is the number of individuals
  # [epitenomic_marks_list] if it is "no", then glm won't use if a mutation is in h3k27ac as a categorical variable, if it is a vector of epigenomic annotation names, epigenomic mark will be added one by one.
  # for example, c("../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed","../other_annotation/epigenomic_annotation/Epigenome_E081_E082_intersection.bed","../other_annotation/epigenomic_annotation/encode_DHS_union.bed"). It is fine to directly use original epigenomics bed files
  # instead of the partitioned files. Because the window_file already has annotations for each window in terms of whetehr they are promoter or coding. 
  # [overlap] the fraction of overlap (the paramter f in bedtools intersect), default is 0.5, For the ease of TADA modeling, it is better to assign it to be 1E-9(1bp)
  # [gene_assign_file] is the file that assigns each window to a gene. for example ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename.bed
  # [rm_nonsyn] TRUE or FALSE. Indicating whether nonysnonymous mutations should be removed when fitting the model. Use TRUE when fitting model using ASD mutations. 
  # [annovar_folder] is where the annovar program is, e.g., . e.g., /media/yuwen/Elements/ANNOVAR/annovar/
  # [annovar_input] is a input file generated from mutations that will be annotated by annovar. The format is chr \t start \t end \t ref_allele \t alt_allele \t index 1-based for starts and ends
  # [gene_prior_file], a file that has prior (derived from posterior and prior)for a gene as a risk gene. genename \t probability e.g. ../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.05.txt
  # [sequence_annotation_list] a list of annotation files for each window in the window_file, e.g., after extracting conservation score for each window. c("../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.phastcons100way")
  # [mut_spidex] is proband_SNV_mut_spedix in the R object. e.g., "../data/debug_region_list_073116_data_matrixtransformed_for_old_code.Rdata"
  # [rr] is the estimated risk. it is vector with relative risk estimates in order from epigenomic_annotations, sequence_annotations to splicing mutations.
  # [coding_bayes_table_file], a file that has Bayes factors calculated from previous coding studies. e.g., "../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_BF_for_de_novo_mutations.txt"
  # [TADA_prior] the prior proportion of genes being non-risk genes. Will used in TADA, default = 0.94
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
  if(epigenomic_marks_list == "no"){
    out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
  }
  else{
    epi_index = 1 # start from the first epi marker, will add 1 to itself after incorporate each epigenomics mark. 
    num_addition_par = 0
    for(epigenomic_marks in epigenomic_marks_list){
      command = paste("bedtools intersect -f ", overlap, " -a ", window_file , " -b ", epigenomic_marks, " -wa > ", paste(prefix,"_temp_overlap_epi.bed", sep = ""),sep = "")
      system(command)
      window_in_epi = read.delim(paste(prefix,"_temp_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      window_in_epi = data.frame(site_index = window_in_epi[,4])
      window_in_epi[,paste("epi", epi_index, sep="")] = 1
      coverage = merge(coverage, window_in_epi, by.x = "site_index", by.y = "site_index", all.x = TRUE)
      coverage[is.na(coverage[,paste("epi", epi_index, sep="")]),paste("epi", epi_index, sep="")] = 0
      epi_index = epi_index + 1
    }
    num_addition_par = length(epigenomic_marks_list) # number of additional feature parameters
    if(sequence_annotation_list != "no"){
      annotation_index = 1
      for(seq_annotation in sequence_annotation_list){
        annotation_raw = read.delim(seq_annotation, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        annotation = data.frame(site_index = annotation_raw[,1])
        annotation[,paste("annotation",annotation_index,sep = "")] = annotation_raw[,4]
        coverage = merge(coverage, annotation, by.x = "site_index", by.y = "site_index")
        annotation_index = annotation_index + 1
      }
      num_addition_par = length(epigenomic_marks_list) + length(sequence_annotation_list) # number of additional feature parameters
    }
    # only use a set of fixed predictors to run glm in order to adjust for mutation rates. 
    out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
    coverage = data.frame(coverage, adjusted_mutrate = out.offset$fitted.values)
    
    # get the number of synonymous mutations that are captured in the windows, which will be used to scale factor for splicing mutations, could be all mutations if rm_nonsy is set to be FALSE
    syn_num = sum(coverage[coverage$coding == 1 & coverage$mut_count ==1,]$coding)
    
    # get mutations that have dpsi_zscore under some certain cutoff, and these mutations will be defined as splicing mutations
    splicing_threshold = quantile(mut_spidex$dpsi_zscore,seq(0,1,0.1),na.rm = TRUE)[2] # lower 10% 
    splicing_mut = mut_spidex[mut_spidex$dpsi_zscore <= splicing_threshold & !is.na(mut_spidex$dpsi_zscore),]
    splicing_mut = as.data.frame(table(splicing_mut$gene))
    colnames(splicing_mut) = c("genename","raw_splicing_mut")
    
    #get splicing mutation number
   # splicing_num = sum(splicing_mut$raw_splicing_mut) # only get 51 for ASD_manuscript_cases
    
    
    gene_assign = read.delim(gene_assign_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(gene_assign) = c("site_index", "genename")
    
    coverage = merge(coverage, gene_assign, by.x = "site_index", by.y = "site_index")
    # now I need to use the prior probability of each gene to generate risk gene indicator and mutation numbers
    gene_prior = read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    #test_count = mapply(function(x){rpois(1,x)},mutrate_temp$mutrate2)
    
    coding_bayes_table = read.delim(coding_bayes_table_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    coding_bayes_table[,2] = log(coding_bayes_table[,2])
    colnames(coding_bayes_table) = c("genename", "logBF_coding")
    
    logBF_noncoding = by(coverage, coverage[,"genename"],
                         function(x) sum(x$mut_count*(x$promoter+x$nf)*(as.matrix(x[,12:(12+num_addition_par-1)])%*%rr[1:num_addition_par]) - x$adjusted_mutrate*(exp((x$promoter+x$nf)*(as.matrix(x[,12:(12+num_addition_par-1)])%*%rr[1:num_addition_par]))-1)))
    
    splicing_rr =rr[num_addition_par+1]
    
    splicing_mutrate = by(coverage, coverage[,"genename"],
                          function(x) sum(x[x$coding == 1,]$adjusted_mutrate))
    splicing_mutrate = data.frame(genename = names(splicing_mutrate), adjusted_mutrate = as.vector(splicing_mutrate))
    splicing_data = merge(splicing_mutrate, splicing_mut, by = "genename", all.x = TRUE)
    splicing_data[is.na(splicing_data$raw_splicing_mut),]$raw_splicing_mut = 0
    
    #scale to get the splicing mutation rate for each gene
    splicing_data$adjusted_mutrate = splicing_data$adjusted_mutrate*sum(splicing_data$raw_splicing_mut)/syn_num
    
    logBF_splicing = by(splicing_data, splicing_data[,"genename"],
                        function(x) (sum(x$raw_splicing_mut)*splicing_rr - sum(x$adjusted_mutrate)*(exp(splicing_rr)-1)))  
    logBF_splicing = data.frame(genename = names(logBF_splicing), logBF_splicing = as.vector(logBF_splicing))
    
    gene_BF_table = data.frame(genename = names(logBF_noncoding), logBF_noncoding = as.vector(logBF_noncoding))
    gene_BF_table = merge(gene_BF_table, logBF_splicing, by = "genename")                           

    gene_BF_table = merge(gene_BF_table, coding_bayes_table, by = "genename")
    gene_BF_table$logBF_all = gene_BF_table$logBF_noncoding + gene_BF_table$logBF_splicing + gene_BF_table$logBF_coding
    gene_BF_table[,c("BF_noncoding","BF_splicing","BF_coding","BF_all")] = exp(gene_BF_table[,c("logBF_noncoding","logBF_splicing","logBF_coding","logBF_all")])
    gene_BF_table = gene_BF_table[order(gene_BF_table$BF_coding, decreasing = TRUE),]
    gene_BF_table$FDR_coding = Bayesian.FDR(gene_BF_table$BF_coding, TADA_p0)$FDR
    gene_BF_table = gene_BF_table[order(gene_BF_table$BF_all, decreasing = TRUE),]
    gene_BF_table$FDR_all = Bayesian.FDR(gene_BF_table$BF_all, TADA_p0)$FDR
  }
  system(paste("rm ", prefix, "_temp*", sep = ""))
  gene_BF_table
}