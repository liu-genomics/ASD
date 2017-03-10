new_control = new.env()
load("../data/debug_region_list_073116_data_matrixtransformed_for_old_code.Rdata", new_control)

prefix <-system("date +%s", intern = TRUE) # prefix for temporary files that will be deleted at the end of the pipeline

write.table(new_control$data_matrix[new_control$data_matrix$Prediction == "dnv_proband" & new_control$data_matrix$Type == "SNV",c(1,2,3,7)], paste(prefix,"temp.bed",sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
#write.table(new_control$data_matrix[new_control$data_matrix$Prediction == "dnv_proband" ,c(1:5,8,7)], paste(prefix,"annovar_input_temp.bed",sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

# SNV_ID for cases or controls. 
SNV_ID = new_control$data_matrix[new_control$data_matrix$Prediction == "dnv_proband" & new_control$data_matrix$Type == "SNV",c(1,2,3,7)]$ID
proband_SNV_mut_spidex <- new_control$mut_spedix[is.element(new_control$mut_spedix$index, SNV_ID),] # this will be paseed as an argument to the function below




mut_file = paste(prefix,"temp.bed",sep = "")
window_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_introns_with_spidex_no_na_genes.50bp_window.bed"
cg_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_introns_with_spidex_no_na_genes.50bp_window.bed.fasta.cg"
sample_size = 314
rm_nonsyn = FALSE
mutrate_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_introns_with_spidex_no_na_genes.50bp_window.bed.mutrate"
spidex_mutrate_prefix = "../other_annotation/spidex_database/mutrate_spedix_score_by_chr/lower_10pct_-1.416/spidex_public_noncommercial_v1_0_collapse_unique_base.extended_3bp.bed.fasta.trinuc.mutrate_with_spedix_gene.lower_10pct"
splicing_cutoff = -1.416
log_splicing_rr = log(2.1)

# This function has windows including intron regions, adjust for mutation rates
generate_BF_for_splicing_mutations <- function(mut_file, 
                                              window_file, 
                                              cg_file, 
                                              mutrate_file, 
                                              sample_size,
                                              rm_nonsyn = FALSE,
                                              annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", 
                                              annovar_input = "no.txt", 
                                              spidex_mutrate_prefix = "../other_annotation/spidex_database/mutrate_spedix_score_by_chr/lower_10pct_-1.416/spidex_public_noncommercial_v1_0_collapse_unique_base.extended_3bp.bed.fasta.trinuc.mutrate_with_spedix_gene.lower_10pct",
                                              splicing_cutoff = -1.416,
                                              log_splicing_rr = log(2.1),
                                              proband_SNV_mut_spidex){
  # [mut_file] is the mutation location file from one dataset. e.g.,../analysis/160229_data_for_analysis.bed, chr \t start \t end \t index,  1-based start and end
  # [window_file] is the file with genome partitioned into windows. e.g., "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_introns_with_spidex_no_na_genes.50bp_window.bed"
  # [cg_file] is the file that has CG-content calculated, matching [window_file]. e.g., "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_introns_with_spidex_no_na_genes.50bp_window.bed.fasta.cg"
  # [mutrate_file] is the file that has mutrate data. e.g., "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_introns_with_spidex_no_na_genes.50bp_window.bed.mutrate"
  # [sample_size] is the number of individuals
  # [rm_nonsyn] TRUE or FALSE. Indicating whether nonysnonymous mutations should be removed when fitting the model. Use TRUE when fitting model using ASD mutations. 
  # [annovar_folder] is where the annovar program is, e.g., . e.g., /media/yuwen/Elements/ANNOVAR/annovar/
  # [annovar_input] is a input file generated from mutations that will be annotated by annovar. The format is chr \t start \t end \t ref_allele \t alt_allele \t index 1-based for starts and ends
  # [spidex_mutrate_prefix] the prefix of a series files that have spidex scores and locations below certain cutoffs. Default is below 10%
  # [splicing_cutoff] cutoff numeric value consistent with [spidex_mutrate_prefix], default is -1.416, which is the lower ten percentile of all the bases with a spedix score
  # [log_splicing_rr] the relative risk of splicing mutations on the log scale.
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

  # now get the splicing mutation coordinates
  # now merge chromosome location with all splicing data
  proband_SNV_mut_spidex_with_location <- merge(proband_SNV_mut_spidex, mut, by.x = "index", by.y = "V4")
  proband_SNV_mut_spidex_with_location <- proband_SNV_mut_spidex_with_location[proband_SNV_mut_spidex_with_location$dpsi_zscore < splicing_cutoff & !is.na(proband_SNV_mut_spidex_with_location$dpsi_zscore),c("V1","V2","V3","index")]
  proband_SNV_mut_spidex_with_location <- data.table(proband_SNV_mut_spidex_with_location)
  colnames(proband_SNV_mut_spidex_with_location) <- c("chr","start","end","index")
  
  # get the adjustion effect for each window by dividing adjusted_mutrate by mutrate
  coverage$adjusting_effect = coverage$adjusted_mutrate/coverage$mutrate
  coverage_by_chr <- split(coverage, coverage$chr)
  output_splicing_logBF <- data.table(genename = character(), logBF_splicing = numeric())
  for(i in 1:length(coverage_by_chr)){
    #write to a bef file
    fwrite(coverage_by_chr[[i]][,c("chr","start","end","adjusting_effect")], paste(prefix, "temp_for_spidex_rate.bed",sep = "_"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    spidex_mutrate_file <- paste(spidex_mutrate_prefix,coverage_by_chr[[i]][1,c("chr")], sep = ".")
    command <- paste("bedtools intersect -a ", spidex_mutrate_file, " -b ", paste(prefix, "temp_for_spidex_rate.bed",sep = "_"), " -wa -wb > ", paste(prefix, "temp_scaling_factor_spidex_rate.bed",sep = "_"), sep = "")
    system(command)
    mutrate_scaling_spidex_rate <- fread(paste(prefix, "temp_scaling_factor_spidex_rate.bed",sep = "_"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    mutrate_scaling_spidex_rate <- mutrate_scaling_spidex_rate[,c(1,2,3,4,6,10)]
    colnames(mutrate_scaling_spidex_rate) <- c("chr","start","end","mutrate","genename","scaling_factor")
    mutrate_scaling_spidex_rate$adjusted_mutrate <- mutrate_scaling_spidex_rate$mutrate * mutrate_scaling_spidex_rate$scaling_factor
    mutrate_scaling_spidex_rate <- mutrate_scaling_spidex_rate[,sum(.SD$adjusted_mutrate),by = c("chr","start","end","genename"), .SDcols = c("adjusted_mutrate")]
    colnames(mutrate_scaling_spidex_rate)[5] <- "collapsed_adjusted_mutrate" # this is the splicing mutation rate summing up allele-specific mutation rates at each base
    # write to a bed file
    fwrite(mutrate_scaling_spidex_rate[,c("chr","start","end")], paste(prefix, "_temp_scalfold_by_chr.bed", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    proband_SNV_mut_spidex_with_location_chr <- proband_SNV_mut_spidex_with_location[chr == coverage_by_chr[[i]][1,c("chr")]] # find the mutations that are in the current chromosome
    # make sure there is indeed at least one splicing mutatin here in this current chromosome.
    if(nrow(proband_SNV_mut_spidex_with_location_chr) != 0){
      # write to a bed file
      fwrite(proband_SNV_mut_spidex_with_location_chr, paste(prefix, "_temp_splicing_mutation.bed", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    
      #calculate splicing mutations coverage on scalfold bed file
      command <- paste("bedtools coverage -a ", paste(prefix, "_temp_splicing_mutation.bed", sep = ""), " -b ", paste(prefix, "_temp_scalfold_by_chr.bed", sep = ""), " > ", paste(prefix, "_temp_scalfold_by_chr.bed.coverageBed", sep = ""), sep = "")
      system(command)
      mutrate_scaling_spidex_coverage <- fread(paste(prefix, "_temp_scalfold_by_chr.bed.coverageBed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      mutrate_scaling_spidex_rate$mut_count <- mutrate_scaling_spidex_coverage$V4
    } else{ # if there is no splicing mutation in this chromosome
      mutrate_scaling_spidex_rate$mut_count <- 0
    }
    # define a function to calculate splicing logBF for each gene
    splicing_logBF_per_gene <- function(splicing_table){
      logBF <- sum(splicing_table$mut_count) * 1 * log_splicing_rr - 2 * sample_size * sum(splicing_table$collapsed_adjusted_mutrate) * (exp(1 * log_splicing_rr) - 1)
    }
    splicing_logBF_by_chr <- mutrate_scaling_spidex_rate[,splicing_logBF_per_gene(.SD), by = "genename", .SDcols = c("mut_count", "collapsed_adjusted_mutrate")]
    colnames(splicing_logBF_by_chr)[2] <- "logBF_splicing"
    output_splicing_logBF <- rbind(output_splicing_logBF, splicing_logBF_by_chr)
  }
  system(paste("rm ", prefix, "_temp*", sep = "")) # remove temp files
  list(splicing_BF_table = output_splicing_logBF)
}
