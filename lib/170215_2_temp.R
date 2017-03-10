# The Version 2 of the function used all splicing mutations, including those mutations that are not overlapped by the window_file.
# This function calculates BF for splicing mutations, and will combine the BF with non-coding BF from the base-level model (output from [bayes_factor_for_each_gene_base_level_from_collapsed_data]$gene_BF_table)
# For now will have to refit the mutation model with the same parameter as used before in order to adjust for the mutation rate for the base-level model. But it doesn't take too much time.
bayes_factor_for_splicing_mut_for_each_gene <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, gene_assign_file, mutrate_baseline = "coding", annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", 
                                                        annovar_input = "no.txt", gene_prior_file, mut_spidex, rr, coding_bayes_table_file){
  # [mut_file] is the mutation location file from one dataset. e.g.,../analysis/160229_data_for_analysis.bed, chr \t start \t end \t index,  1-based start and end
  # [window_file] is the file with genome partitioned into windows. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed
  # [cg_file] is the file that has CG-content calculated, matching [window_file]. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.cg
  # [mutrate_file] is the file that has mutrate data. e.g., "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate"
  # [sample_size] is the number of individuals
  # [gene_assign_file] is the file that assigns each window to a gene. for example ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename.bed
  # [mutrate_baseline] Use what types of mutations to scale for the splicing mutation rate. If "syn", then only synonymous mutations will be used to get the scaling factor for infering splicing mutations.
  # [annovar_folder] is where the annovar program is, e.g., . e.g., /media/yuwen/Elements/ANNOVAR/annovar/
  # [annovar_input] is a input file generated from mutations that will be annotated by annovar. The format is chr \t start \t end \t ref_allele \t alt_allele \t index 1-based for starts and ends
  # [gene_prior_file], a file that has prior (derived from posterior and prior)for a gene as a risk gene. genename \t probability e.g. ../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.05.txt
  # [mut_spidex] is the proband_SNV_mut_spedix in the R object. e.g., "../data/debug_region_list_073116_data_matrixtransformed_for_old_code.Rdata"
  # [rr] is the estimated risk. it is vector with relative risk estimates in order from epigenomic_annotations, sequence_annotations to splicing mutations. Remeber to use only the mutations that are either in cases or controls depending on the needs. 
  # The papramter is in the log scale. For example, if the real rr is 3, the parameter that needs to be set here is log(3).
  # [coding_bayes_table_file], a file that has Bayes factors calculated from previous coding studies. e.g., "../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_BF_for_de_novo_mutations.txt"
  prefix <- system("date +%s", intern = TRUE) # prefix for temporary files that will be deleted at the end of the pipeline
  
  mut <- fread(mut_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
#   # command to transform the start base to 0-based in order to use bedtools to do overlap
#   command <- paste("awk {'print $1\"\t\"$2-1\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$7'} ", annovar_input, " > ", paste(prefix, "_annovar_for_bedtools.bed",sep = ""),sep = "")
#   system(command)
#   command <- paste("grep coding ", window_file, " > ", paste(prefix, "_coding_window_temp.bed", sep = ""), sep = "")
#   system(command)
#   command <- paste("bedtools intersect -a ",paste(prefix, "_annovar_for_bedtools.bed",sep = ""), " -b ", paste(prefix, "_coding_window_temp.bed", sep = ""), " -wa | sort | uniq | awk {'print $1\"\t\"$2+1\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$7'} > ",paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""), sep = "")
#   system(command)
#   command <- paste("perl ", annovar_folder, "annotate_variation.pl -geneanno -buildver hg19 ", paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""), " ", annovar_folder, "humandb/", sep = "")
#   system(command)
#   coding_anno <- fread(paste(paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""),".exonic_variant_function", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  mut[,2] <- mut[,2] - 1  # change 1-based start to 0-based start. 
  fwrite(mut, paste(prefix,"_temp.bed", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  command <- paste("bedtools coverage -a ", paste(prefix,"_temp.bed", sep = ""), " -b ", window_file, " > ", paste(prefix,"_temp_coverage.bed", sep = ""), sep = "")
  system(command)
  coverage <- fread(paste(prefix,"_temp_coverage.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  coverage <- coverage[,1:5]
  colnames(coverage) <- c("chr","start","end","site_index","mut_count")
  coverage$window_size <- coverage$end - coverage$start # some regions don't have full length of, say, 50bp
  window <- fread(window_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  cg <- fread(cg_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(cg) <- c("site_index", "cg")
  coverage <- coverage[cg, on="site_index"]
  
  binerize_adjust_feature <- function(x) {return(list(coding = as.numeric(grepl("coding",x)),
                                                      promoter = as.numeric(grepl("promoter",x)),
                                                      nf = as.numeric(grepl("nf",x))))
  }
  
  coverage[,c("coding","promoter","nf") := binerize_adjust_feature(site_index), with = FALSE]
  
  mutrate <- fread(mutrate_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  mutrate<-mutrate[,c("V1","V4"),with = FALSE]
  colnames(mutrate) <- c("site_index", "mutrate")
  coverage <-coverage[mutrate, on = "site_index"]
  coverage <- coverage[mutrate !=0]
  #fit mutation rate model using a fixed set of features not including histone modification marks. 
  out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
  
  coverage$adjusted_mutrate = out.offset$fitted.values
    
  # get the number of coding mutations, and use the number of these mutations as a reference to infer splicing mutation rate. 
  if(mutrate_baseline == "coding"){
    baseline_mut_num = sum(coverage[coverage$coding == 1 & coverage$mut_count ==1,]$coding)
  }
  # get mutations that have dpsi_zscore under some certain cutoff, and these mutations will be defined as splicing mutations
  splicing_threshold = quantile(mut_spidex$dpsi_zscore,seq(0,1,0.1),na.rm = TRUE)[2] # lower 10% 
  splicing_mut = mut_spidex[mut_spidex$dpsi_zscore <= splicing_threshold & !is.na(mut_spidex$dpsi_zscore),]
  splicing_mut = as.data.table(table(splicing_mut$gene))
  colnames(splicing_mut) = c("genename","raw_splicing_mut")
    
  gene_assign <- fread(gene_assign_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_assign) = c("site_index", "genename")
  # add gene name information to 50-bp windows
  coverage <- gene_assign[coverage, on = "site_index"]
  # now get the reference mutation rate, this depends on which number we used as a reference to infer the splicing mutation rates. If "coding" is chosen, then will sum all coding mutation rates for every gene here. 
  splicing_mutrate <- by(coverage, coverage[,"genename"],function(x) sum(x[x$coding == 1,]$adjusted_mutrate))
  splicing_mutrate = data.frame(genename = names(splicing_mutrate), adjusted_mutrate = as.vector(splicing_mutrate))
  splicing_data = merge(splicing_mutrate, splicing_mut, by = "genename", all.x = TRUE)
  splicing_data[is.na(splicing_data$raw_splicing_mut),]$raw_splicing_mut = 0
  
  
  
  
  
  
  # now I need to use the prior probability of each gene to generate risk gene indicator and mutation numbers
    gene_prior = read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    #test_count = mapply(function(x){rpois(1,x)},mutrate_temp$mutrate2)
    
    coding_bayes_table = read.delim(coding_bayes_table_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    coding_bayes_table[,2] = log(coding_bayes_table[,2])
    colnames(coding_bayes_table) = c("genename", "logBF_coding")
    
    logBF_noncoding = by(coverage, coverage[,"genename"],
                         function(x) sum(x$mut_count*(x$promoter+x$nf)*(as.matrix(x[,12:(12+num_addition_par-1)])%*%rr[1:num_addition_par]) - x$adjusted_mutrate*(exp((x$promoter+x$nf)*(as.matrix(x[,12:(12+num_addition_par-1)])%*%rr[1:num_addition_par]))-1)))
    gene_with_epi = by(coverage, coverage[,"genename"],
                       function(x) sum(apply(x[,c(9,10)]*x[,12:(12+length(epigenomic_marks_list)-1)],1,sum)))
    
    genes_with_no_epi = names(gene_with_epi[gene_with_epi == 0])
    
    splicing_rr =rr[num_addition_par+1]
    
    
    
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
  list(gene_BF_table = gene_BF_table, genes_with_no_epi = genes_with_no_epi)
}