# mut_file = paste(prefix,"temp.bed",sep = "")
# window_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_introns_with_spidex_no_na_genes.50bp_window.bed"
# cg_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_introns_with_spidex_no_na_genes.50bp_window.bed.fasta.cg"
# mutrate_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_introns_with_spidex_no_na_genes.50bp_window.bed.mutrate"
# sample_size = 162
# overlap = 0.5
# rm_nonsyn = FALSE
# annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/"
# annovar_input = "no.txt"
# gene_prior_file = "../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.94.txt"
# report_proportion = 1000/18665
# spidex_mutrate_prefix = "../other_annotation/spidex_database/mutrate_spedix_score_by_chr/lower_10pct_-1.416/spidex_public_noncommercial_v1_0_collapse_unique_base.extended_3bp.bed.fasta.trinuc.mutrate_with_spedix_gene.lower_10pct"
# splicing_cutoff = -1.416
# proband_SNV_mut_spidex = proband_SNV_mut_spidex





#This patch version corrects the [splicing_info] output from [adjust_mutation_rate_window_to_base_compact_more_feature_splicing_hybrid] 
# [adjust_mutation_rate_window_to_base_compact_more_feature_splicing_hybrid] outputs a [splicing_info] that doesn't have the correct information
adjust_mutation_rate_window_to_base_compact_more_feature_splicing_hybrid_patch <- function(mut_file, window_file, cg_file, mutrate_file, sample_size,
                                                                                     overlap = 0.5, rm_nonsyn = FALSE, annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", annovar_input = "no.txt", 
                                                                                     gene_prior_file, report_proportion, 
                                                                                     spidex_mutrate_prefix, splicing_cutoff, proband_SNV_mut_spidex){
  # [mut_file] is the mutation location file from one dataset. e.g.,../analysis/160229_data_for_analysis.bed, chr \t start \t end \t index,  1-based start and end
  # [window_file] is the file with genome partitioned into windows. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed
  # [cg_file] is the file that has CG-content calculated, matching [window_file]. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.cg
  # [mutrate_file] is the file that has mutrate data. e.g., "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate"
  # [sample_size] is the number of individuals
 
  # [overlap] the fraction of overlap (the paramter f in bedtools intersect), default is 0.5, For the ease of TADA modeling, it is better to assign it to be 1E-9(1bp)
  # [rm_nonsyn] TRUE or FALSE. Indicating whether nonysnonymous mutations should be removed when fitting the model. Use TRUE when fitting model using ASD mutations. 
  # [annovar_folder] is where the annovar program is, e.g., . e.g., /media/yuwen/Elements/ANNOVAR/annovar/
  # [annovar_input] is a input file generated from mutations that will be annotated by annovar. The format is chr \t start \t end \t ref_allele \t alt_allele \t index 1-based for starts and ends
  # [gene_prior_file], a file that has prior (derived from posterior and prior)for a gene as a risk gene. genename \t probability e.g. ../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.94.txt
  # [report_proportion] Choose the top X% TADA genes to report. 
  # [spidex_mutrate_prefix] the prefix of a series files that have spidex scores and locations below certain cutoffs. Default is below 10%
  # [splicing_cutoff] cutoff numeric value consistent with [spidex_mutrate_prefix], default is -1.416, which is the lower ten percentile of all the bases with a spedix score
  # [proband_SNV_mut_spidex], spidex annotation of mutations, e.g., new_control$mut_spedix[is.element(new_control$mut_spedix$index, SNV_ID),], which is the output of spidex score
  
  prefix <- system("date +%s", intern = TRUE) # prefix for temporary files that will be deleted at the end of the pipeline
  
  
  
  mut <- fread(mut_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  if(rm_nonsyn){
    # command to transform the start base to 0-based in order to use bedtools to do overlap
    command <- paste("awk {'print $1\"\t\"$2-1\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$7'} ", annovar_input, " > ", paste(prefix, "_annovar_for_bedtools.bed",sep = ""),sep = "")
    system(command)
    command <- paste("grep coding ", window_file, " > ", paste(prefix, "_coding_window_temp.bed", sep = ""), sep = "")
    system(command)
    command <- paste("bedtools intersect -a ",paste(prefix, "_annovar_for_bedtools.bed",sep = ""), " -b ", paste(prefix, "_coding_window_temp.bed", sep = ""), " -wa | sort | uniq | awk {'print $1\"\t\"$2+1\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$7'} > ",paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""), sep = "")
    system(command)
    command <- paste("perl ", annovar_folder, "annotate_variation.pl -geneanno -buildver hg19 ", paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""), " ", annovar_folder, "humandb/", sep = "")
    system(command)
    coding_anno <- fread(paste(paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""),".exonic_variant_function", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    removed_id <- coding_anno[coding_anno[,2] != "synonymous SNV",10]
    # now remove non-synonymous mutations (every coding mutations that are not synonymous mutations)
    mut <- mut[!is.element(mut[,4], removed_id),]
  }
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
                                                      nf = as.numeric(grepl("nf",x)),
                                                      spidex_intron = as.numeric(grepl("spdx_intron",x))))
  }
  
  coverage[,c("coding","promoter","nf","spdx_intron") := binerize_adjust_feature(site_index), with = FALSE]
  
  mutrate <- fread(mutrate_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  mutrate<-mutrate[,c("V1","V4"),with = FALSE]
  colnames(mutrate) <- c("site_index", "mutrate")
  coverage <-coverage[mutrate, on = "site_index"]
  coverage <- coverage[mutrate !=0]
  #fit mutation rate model using a fixed set of features not including histone modification marks. 
  out.offset <- glm(mut_count ~ coding+promoter+spdx_intron+cg+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
  # get the adjusted mutation rate per base per individual
  coverage$adjusted_mutrate = out.offset$fitted.values/(2*sample_size)
  

  # get the piror probability of genes.
  gene_prior <- read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_prior) <- c("genename", "prior")
  # change prior to the posterior probability of alternative (risk) hypothesis
  gene_prior$prior <- 1- gene_prior$prior
 
  
  # remove genes based on TADA prior probability
  if(report_proportion !=1){
    genes_for_report <- gene_prior[order(gene_prior[,2],decreasing = TRUE),1]
    genes_for_report <- genes_for_report[1:floor(length(genes_for_report)*report_proportion)]
  }else{
    genes_for_report  <- gene_prior[,1] # choose all the genes in TADA coding table 
  }
  
 
  
  # Now will parse splicing information in a compact form, e.g., data_partition[[1]][[1]]
  # feature_vector would have two different choices, either 0 or 1. 
  # get the adjustion effect for each window by dividing adjusted_mutrate by mutrate
  coverage$adjusting_effect = coverage$adjusted_mutrate/coverage$mutrate
  coverage_by_chr <- split(coverage, coverage$chr)
  
  # now get the splicing mutation coordinates
  # now merge chromosome location with all splicing data
  proband_SNV_mut_spidex_with_location <- merge(proband_SNV_mut_spidex, mut, by.x = "index", by.y = "V4")
  proband_SNV_mut_spidex_with_location <- proband_SNV_mut_spidex_with_location[proband_SNV_mut_spidex_with_location$dpsi_zscore < splicing_cutoff & !is.na(proband_SNV_mut_spidex_with_location$dpsi_zscore),c("V1","V2","V3","index")]
  proband_SNV_mut_spidex_with_location <- data.table(proband_SNV_mut_spidex_with_location)
  colnames(proband_SNV_mut_spidex_with_location) <- c("chr","start","end","index")
  
  splicing_data_partition <-list()
  
  for(i in 1:length(coverage_by_chr)){
    #write to a bed file, didn't remove mutations that have other annotations for now. 
    fwrite(coverage_by_chr[[i]][,c("chr","start","end","adjusting_effect")], paste(prefix, "temp_for_spidex_rate.bed",sep = "_"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    spidex_mutrate_file <- paste(spidex_mutrate_prefix,coverage_by_chr[[i]][1,c("chr")], sep = ".")
    command <- paste("bedtools intersect -a ", spidex_mutrate_file, " -b ", paste(prefix, "temp_for_spidex_rate.bed",sep = "_"), " -wa -wb > ", paste(prefix, "temp_scaling_factor_spidex_rate.bed",sep = "_"), sep = "")
    system(command)
    mutrate_scaling_spidex_rate <- fread(paste(prefix, "temp_scaling_factor_spidex_rate.bed",sep = "_"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    mutrate_scaling_spidex_rate <- mutrate_scaling_spidex_rate[,c(1,2,3,4,6,10)]
    colnames(mutrate_scaling_spidex_rate) <- c("chr","start","end","mutrate","genename","scaling_factor")
    mutrate_scaling_spidex_rate$adjusted_mutrate <- mutrate_scaling_spidex_rate$mutrate * mutrate_scaling_spidex_rate$scaling_factor*2*sample_size # remember to add sample size in. 
    mutrate_scaling_spidex_rate <- mutrate_scaling_spidex_rate[,sum(.SD$adjusted_mutrate),by = c("chr","start","end","genename"), .SDcols = c("adjusted_mutrate")]
    colnames(mutrate_scaling_spidex_rate)[5] <- "adjusted_base_mutrate" # this is the splicing mutation rate summing up allele-specific mutation rates at each base
    
    # get splicing mutation coverage at the base level
    fwrite(mutrate_scaling_spidex_rate, paste(prefix,"_temp_splicing_for_getting_mutation_coverage.bed", sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    proband_SNV_mut_spidex_with_location_chr <- proband_SNV_mut_spidex_with_location # will calculate the coverage of one chromosome using all mutations anyway
    # previously I only use mutations that are on the current chromosome, but there are cases where no mutations are on this chromosome.
    # To account for this, I skip the steps from here in the current iterations, which turn out to be totally wrong.
    # This leas to overestimate the RR of splicing mutations, because whenever there is no splicing mutations, I removed the whole chromosome in the RR estimation.
    
    fwrite(proband_SNV_mut_spidex_with_location_chr, paste(prefix, "_temp_splicing_mutation.bed", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    command <- paste("bedtools coverage -a ", paste(prefix, "_temp_splicing_mutation.bed", sep = ""), " -b ", paste(prefix,"_temp_splicing_for_getting_mutation_coverage.bed", sep = ""), " > ", paste(prefix,"_temp_splicing_mut_coverage.bed", sep = ""), sep = "")
    system(command)
    mutrate_scaling_spidex_rate <- fread(paste(prefix,"_temp_splicing_mut_coverage.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    mutrate_scaling_spidex_rate <- mutrate_scaling_spidex_rate[,c(1,2,3,4,5,6)]
    colnames(mutrate_scaling_spidex_rate) <- c("chr","start","end","genename","adjusted_base_mutrate","mut_count")
    
    mutrate_scaling_spidex_rate$splicing_feature = 1 # add a feature to everybase, because here, every base does have at least one splicing mutation if mutated
      #only retain genes that are considered in [genes_for_report], which is defined earlier by [report_proportion]
    mutrate_scaling_spidex_rate <- mutrate_scaling_spidex_rate[is.element(mutrate_scaling_spidex_rate$genename, genes_for_report),]
    if(nrow(mutrate_scaling_spidex_rate) == 0) next # if no genes (in genes for report) on this current chromosome have splicing mutations, exit the current iteration without adding anything to [splicing_data_partition] and go to the next one
    
    partition_splicing_feature <- function(pbg){
      # input is one element of the list of partition_by_gene
      pbg_split <- split(pbg, pbg$splicing_feature,drop = TRUE)
      feature_combination_number <- length(pbg_split)
      # this function below is different from the function used in dealing with dataset without reading by chunk. Here, prior is not incoporated at this step.
      info_for_each_feature <- function(feature_set){
        list(feature_vector = c(1), sum_mut_rate_count = sum(feature_set$mut_count*log(feature_set$adjusted_base_mutrate)), sum_mut_rate = sum(feature_set$adjusted_base_mutrate), sum_mut_count = sum(feature_set$mut_count), log_fcount = sum(log(factorial(feature_set$mut_count))))
      }
      sapply(pbg_split, info_for_each_feature,simplify = FALSE)
    }
    
    #then partition by gene for the current chromosome
    mutrate_scaling_spidex_rate <- split(mutrate_scaling_spidex_rate, mutrate_scaling_spidex_rate$genename)
    # then partition by feature configuration for each gene in the current chunk
    mutrate_scaling_spidex_rate <- sapply(mutrate_scaling_spidex_rate, partition_splicing_feature, simplify = FALSE)
    # add into [splicing_data_partition]
    splicing_data_partition <- append(splicing_data_partition, mutrate_scaling_spidex_rate)
    
  }
  system(paste("rm ", prefix, "_temp*", sep = ""))
  list(splicing_info = splicing_data_partition)
}
