mut_collapsed_data_list = list(adjusted_base_level_rr_1$base_info, adjusted_base_level_rr_2$base_info, adjusted_base_level_rr_3$base_info, adjusted_base_level_rr_4$base_info, adjusted_base_level_rr_5$base_info)
features = c(1)
rr = 0.6944592
coding_bayes_table_file = "../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_BF_for_de_novo_mutations.txt"
TADA_p0 = 0.94


# all the features that need to be used are already encoded in a collapsed form the in the input data (output from [adjust_mutation_rate_window_to_base_compact_more_feature_v2]$info_base).
# so if splicing info hasn't been included in the data in collapsed form, splicing mutations won't be carried in the analysis.
# Takes categorized data from multiple studies.
bayes_factor_for_each_gene_base_level_from_collapsed_data_multiple_studies <- function(mut_collapsed_data_list,features, rr, coding_bayes_table_file, TADA_p0 = 0.94){
  # [mut_collapsed_data_list] is a list of (output from [adjust_mutation_rate_window_to_base_compact_more_feature_splicing_hybrid]$info_base).
  # [features] determines which features in the mut_collapsed_data are going to be used for calculating bayes factors. e.g., if [feature_vector] is c(1,3). Then the first and third element of feature_vector in mut_collapsed_data will be selected
  # [rr] is the estimated risk of selected features. Estiamted from [estimate_effect_size_for_simulation_data_mixture_with_categorization_for_compact_data] using adjust_mutation_rate_window_to_base_compact_more_feature_v2]$info_base
  # rr is on log scale
  # [coding_bayes_table_file], a file that has Bayes factors calculated from previous coding studies. e.g., "../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_BF_for_de_novo_mutations.txt"
  # [TADA_p0] the prior proportion of genes being non-risk genes. Will used in TADA, default = 0.94
  
  coding_bayes_table <- fread(coding_bayes_table_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  coding_bayes_table[,2] = log(coding_bayes_table[,2])
  colnames(coding_bayes_table) = c("genename", "logBF_coding")
  
  cal_logP_Zg1 <- function(data_partition_element){
    cal_logP_Zg1_level2 <-function(data_partition_element_level2){
      data_partition_element_level2[[2]]+data_partition_element_level2[[4]]*data_partition_element_level2[[1]][features]%*%rr-data_partition_element_level2[[3]]*exp(data_partition_element_level2[[1]][features]%*%rr)-data_partition_element_level2[[5]]
    }
    sum(sapply(data_partition_element, cal_logP_Zg1_level2))
  }
  
  cal_logP_Zg0 <- function(data_partition_element){
    cal_logP_Zg0_level2 <-function(data_partition_element_level2){
      data_partition_element_level2[[2]]-data_partition_element_level2[[3]]-data_partition_element_level2[[5]]
    }
    sum(sapply(data_partition_element, cal_logP_Zg0_level2))
  }
  
  
  logP_Zg1 <- sapply(mut_collapsed_data, cal_logP_Zg1)
  logP_Zg0 <- sapply(mut_collapsed_data, cal_logP_Zg0)
  
  # define a function to flag genes that don't have any bases with any features that are considered to affect relative risk in the model. 
  # for example, if we only have Brain K27ac in the model, then we will flag genes that don't have any non-coding bases covered by this histone marks.
  # And the non-coding bayes factor for these genes would be always 1. Because the null model and the alternative model would be exactly the same. 
  no_feature_flag <- function(data_partition_element){
    no_feature_flag_level2 <-function(data_partition_element_level2){
      sum(data_partition_element_level2[[1]][features]^2)
    }
    as.numeric(sum(sapply(data_partition_element, no_feature_flag_level2)) >0) # if there are any bases that have any features that are considered to have rr in the model
  }
  
  noncoding_rr_feature_flag <- sapply(mut_collapsed_data, no_feature_flag)
  # get the names of genes without bases that have any rr features
  genes_without_rr_feature <- names(noncoding_rr_feature_flag[noncoding_rr_feature_flag == 0])
  
  
  logBF_noncoding <-data.table(genename = names(logP_Zg1), logBF_noncoding = logP_Zg1 - logP_Zg0)
  # for each gene sum the logBF_noncoding, I did this because there is possibility that a same gene was partitioned before collapsing feature data. 
  temp <- by(logBF_noncoding, logBF_noncoding$genename, function(x) {sum(x$logBF_noncoding)})
  logBF_noncoding <- data.table(genename = names(temp), logBF_noncoding = as.vector(temp))
  
  # generate a gene_BF_table to hold up everything, including bayse factors for coding and non-coding mutations
  gene_BF_table <- logBF_noncoding
  # add in coding_BF info
  gene_BF_table <- gene_BF_table[coding_bayes_table, on = "genename"]
  
  #get genes that appear in the TADA table, but somehow don't have non-coding bases(no matter if have rr features or not)
  
  genes_without_noncoding_bases <- gene_BF_table[is.na(logBF_noncoding)]$genename
  
  # change NA to 0, NA would be generated if no non-coding annotations/mutations of a gene were recorded in  mut_collapsed_data.
  # However, in mut_collapsed_data, a gene would still be included if it has non-coding bases included even in the situation where non of the bases overlap with any epigenomic marks considered or are marked by any other base-level annotation as deletrious mutations
  # For such a gene, need to consider if we want to call it as novel gene, if the combined FDR becomes a little bit less than a cutoff and its coding FDR is a little greater than this cutoff.
  gene_BF_table[is.na(logBF_noncoding)]$logBF_noncoding = 0
  # add coding and non-coding logBF
  gene_BF_table$logBF_all <- gene_BF_table$logBF_noncoding + gene_BF_table$logBF_coding
  gene_BF_table[,c("BF_noncoding","BF_coding","BF_all")] = exp(gene_BF_table[,c("logBF_noncoding","logBF_coding","logBF_all")])
  gene_BF_table = gene_BF_table[order(gene_BF_table$BF_coding, decreasing = TRUE),]
  gene_BF_table$FDR_coding = Bayesian.FDR(gene_BF_table$BF_coding, TADA_p0)$FDR
  gene_BF_table = gene_BF_table[order(gene_BF_table$BF_all, decreasing = TRUE),]
  gene_BF_table$FDR_all = Bayesian.FDR(gene_BF_table$BF_all, TADA_p0)$FDR
  
  
  list(gene_BF_table = gene_BF_table, genes_with_no_epi = genes_without_rr_feature)
}

