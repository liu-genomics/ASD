#function to calculate burden after correcting CG content
burden_CG_corrected <-function(mutation_table, data_matrix, gene_mut_table, CG_content_table, geneset, annotation_list=c("NA"), annotation_cutoff=0,
                               CG_window_size = "CG_50bp", mut_type = "SNV", ref = "self", all_genes = "no", annotation = "yes"){
  #[mutation_table] is the mutation data.frame in mutation.Rdata
  #[data_matrix] is the data_matrix data.frame in mutation.Rdata
  #[gene_mut_table] is the output of {generate_mut_enhancer_pair_for_burden_analysis}, the 1st column is mutation ID, second column is gene name
  #[CG_content_table] is the CG_content data.frame in the mutation.Rdata
  #[geneset] is the geneset that is going to be focused
  #[annotation_list] is a vector giving the annotation that needs to be used for example c("GERP", "Eigen", "CADD13_PHRED", "phyloP46wayAllElements")
  #[annotation_cutoff] is a vecotr giving the cutoff of each chosen annotation, needs to have the same annotation sequence as [annotation_list]
  #[CG_window_size] is the window size to calculate CG content around mutations. four different choices "CG_50bp", "CG_100bp", "CG_200bp", "CG_500bp"
  #[mut_type], currently it only works for SNV. 
  #[ref] "self" or "all_mutations", if self, only use all the mutations that are in epigenomic marks as baseline, suitable for 160229 data and Scherer data
  #[all_genes] is "yes" or "no", when is set to "no", only genes in the geneset will be considered. otherwise, geneset information is 
  #[annotation] is "yes" if using base level annotation, otherwise, 
  
  # first only take SNV mutations here, temp_data is the all mutations that are not SNVs for here
  temp_data = mutation_table[mutation_table$mut_type == mut_type,]
  data = data.frame(y = temp_data$phenotype, index = temp_data$index, annotation = 0 , CG = CG_content_table[mutation_table$mut_type == mut_type,CG_window_size] )
  
  # then get burden analysis done
  if(all_genes == "no" & ref == "self" & annotation == "yes"){ # with base-level annotation and geneset, compared to epigenomic baseline
    data = data[is.element(data$index, gene_mut_table[,1]),] # only include mutations that are in epigenomic marks
    ID_in_epi_and_gene = gene_mut_table[is.element(gene_mut_table[,2], geneset),1]
    annotation_index = ID_in_epi_and_gene
    for(i in 1:length(annotation_list)){
      annotation_index = intersect(annotation_index, data_matrix[data_matrix[,annotation_list[i]] >= annotation_cutoff[i] & !is.na(data_matrix[,annotation_list[i]]),]$ID)
    }
    data[is.element(data$index, annotation_index),]$annotation = 1
  }
  
  if(all_genes == "no" & ref == "total_ratio" & annotation == "yes"){ # with base-level annotation and geneset, compared to all SNV baseline
    #data = data[is.element(data$index, gene_mut_table[,1]),] # only include mutations that are in epigenomic marks
    ID_in_epi_and_gene = gene_mut_table[is.element(gene_mut_table[,2], geneset),1]
    annotation_index = ID_in_epi_and_gene
    for(i in 1:length(annotation_list)){
      annotation_index = intersect(annotation_index, data_matrix[data_matrix[,annotation_list[i]] >= annotation_cutoff[i] & !is.na(data_matrix[,annotation_list[i]]),]$ID)
    }
    data[is.element(data$index, annotation_index),]$annotation = 1
  }
  
  else if(all_genes == "no" & ref == "self" & annotation == "no"){ #without base-level annotation with geneset, compared to epigenome baseline
    data = data[is.element(data$index, gene_mut_table[,1]),] # only include mutations that are in epigenomic marks
    ID_in_epi_and_gene = gene_mut_table[is.element(gene_mut_table[,2], geneset),1]
    annotation_index = ID_in_epi_and_gene
    data[is.element(data$index, annotation_index),]$annotation = 1
  }
  
  else if(all_genes == "no" & ref == "total_ratio" & annotation == "no"){ #without base-level annotation with geneset, compared to all SNVs baseline
    #data = data[is.element(data$index, gene_mut_table[,1]),] # only include mutations that are in epigenomic marks
    ID_in_epi_and_gene = gene_mut_table[is.element(gene_mut_table[,2], geneset),1]
    annotation_index = ID_in_epi_and_gene
    data[is.element(data$index, annotation_index),]$annotation = 1
  }
  
  else if(all_genes == "yes" & ref == "self" & annotation == "yes"){ # for burden with annoation but no geneset , compared to epigenomic baseline 
    data = data[is.element(data$index, gene_mut_table[,1]),] # only include mutations that are in epigenomic marks
    ID_in_epi_and_gene = data$index
    annotation_index = ID_in_epi_and_gene
    for(i in 1:length(annotation_list)){
      annotation_index = intersect(annotation_index, data_matrix[data_matrix[,annotation_list[i]] >= annotation_cutoff[i] & !is.na(data_matrix[,annotation_list[i]]),]$ID)
    }
    data[is.element(data$index, annotation_index),]$annotation = 1
  }
  
  else if(all_genes == "yes" & ref == "total_ratio" & annotation == "yes"){ # for burden with annoation but no geneset , compared to SNV baseline 
    #data = data[is.element(data$index, gene_mut_table[,1]),] # only include mutations that are in epigenomic marks
    ID_in_epi_and_gene = intersect(data$index,gene_mut_table[,1]) 
    annotation_index = ID_in_epi_and_gene
    for(i in 1:length(annotation_list)){
      annotation_index = intersect(annotation_index, data_matrix[data_matrix[,annotation_list[i]] >= annotation_cutoff[i] & !is.na(data_matrix[,annotation_list[i]]),]$ID)
    }
    data[is.element(data$index, annotation_index),]$annotation = 1
  }
  
  else if(all_genes == "yes" & ref == "self" & annotation == "no"){ #without base-level annotation without geneset, compared to epigenome baseline
    data = data[is.element(data$index, gene_mut_table[,1]),] # only include mutations that are in epigenomic marks
    #ID_in_epi_and_gene = gene_mut_table[is.element(gene_mut_table[,2], geneset),1]
    annotation_index = data$index
    data[is.element(data$index, annotation_index),]$annotation = 1
  }
  
  else if(all_genes == "yes" & ref == "total_ratio" & annotation == "no"){ #without base-level annotation without geneset, compared to all SNVs baseline
    #data = data[is.element(data$index, gene_mut_table[,1]),] # only include mutations that are in epigenomic marks
    ID_in_epi_and_gene = intersect(data$index,gene_mut_table[,1]) 
    annotation_index = ID_in_epi_and_gene
    data[is.element(data$index, annotation_index),]$annotation = 1
  }
  
  data$y = factor(data$y, levels = c("control","ASD"))
  model = glm(y~annotation+CG, family = binomial(link='logit'), data = data)
  odds_ratio = as.numeric(model$coefficients["annotation"])
  if(is.na(model$coefficients["annotation"])){# in some situations will be NA, if all the data has annotation as 1 or other situations
    odds_ratio = 0
    pvalue = 1
    sd = 0
  }
  else {
    pvalue = summary(model)$coeff["annotation","Pr(>|z|)"] # p-value here is the two-sided pvalue
    sd = summary(model)$coeff["annotation","Std. Error"]
  }
  c(nrow(data[data$annotation ==1 & data$y == "ASD",]), nrow(data[data$annotation ==1 & data$y == "control",]), odds_ratio, exp(odds_ratio), exp(odds_ratio), pvalue,
    exp(odds_ratio-1.96*sd), exp(odds_ratio-1.96*sd))
}

show_burden_for_selected_features_with_geneset_CG_corrected <-function(mut_file_name, input_type = "file", mut_type = "noncoding",ref = "self", syn_A = length(mutation$ASD_effective_SNV_ID), syn_C = length(mutation$control_effective_SNV_ID), geneset){
  #input_type could be file or a table matching regions with mutations. 
  # mut_type could be "noncoding" or "coding
  # ref = "total_ratio" means using total number of SNVs as background, and use fisher.test to 
  # get pvalues
  # syn_A is the number of SNVs in ASD
  # syn_C is the number of SNVs in control
  # geneset is the geneset that will be tested, predefined in upper section of this code. 
  if(input_type == "file"){
    mut = read.delim(mut_file_name, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  }
  else if (input_type == "table"){
    mut = mut_file_name
  }
  if(mut_type == "noncoding"){
    
    output = data.frame(ASD = numeric(0), control = numeric(0), frequency_ratio_to_background = numeric(0), burden_to_baseline = numeric(0), odds_ratio = numeric(0), pvalue = numeric(0), lowerbound = numeric(0), upperbound = numeric(0))
    output['baseline',] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                              geneset = c("all_genes"), annotation_list = c("NA"), annotation_cutoff = 0, CG_window_size = "CG_50bp", mut_type = "SNV", 
                                              ref = "total_ratio", all_genes = "yes", annotation = "no") # here for baseline, ref will always be "total ratio" or NA will be generate
    
    output[deparse(substitute(geneset)),] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                                                geneset = geneset, annotation_list = c("NA"), annotation_cutoff = 0, CG_window_size = "CG_50bp", mut_type = "SNV", 
                                                                ref = ref, all_genes = "no", annotation = "no")
    output['+motif',] = c(0,0,0,0) # haven't updated motif burden here
    
    output['baseline+CADD_95pct',] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                                         geneset = c("all_genes"), annotation_list = c("CADD13_PHRED"), annotation_cutoff = as.numeric(quantile(mutation$CADD_score[,1], c(0.95),na.rm = TRUE)),
                                                         CG_window_size = "CG_50bp", mut_type = "SNV", ref = ref, all_genes = "yes", annotation = "yes")
    output['baseline+CADD_90pct',] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                                         geneset = c("all_genes"), annotation_list = c("CADD13_PHRED"), annotation_cutoff = as.numeric(quantile(mutation$CADD_score[,1], c(0.9),na.rm = TRUE)),
                                                         CG_window_size = "CG_50bp", mut_type = "SNV", ref = ref, all_genes = "yes", annotation = "yes")
    output['baseline+CADD_gt10',] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                                        geneset = c("all_genes"), annotation_list = c("CADD13_PHRED"), annotation_cutoff = 10,
                                                        CG_window_size = "CG_50bp", mut_type = "SNV", ref = ref, all_genes = "yes", annotation = "yes")
    output[paste(deparse(substitute(geneset)),"+CADD_95pct", sep = ""),] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                                                                               geneset = geneset, annotation_list = c("CADD13_PHRED"), annotation_cutoff = as.numeric(quantile(mutation$CADD_score[,1], c(0.95),na.rm = TRUE)),
                                                                                               CG_window_size = "CG_50bp", mut_type = "SNV", ref = ref, all_genes = "no", annotation = "yes")
    output[paste(deparse(substitute(geneset)),"+CADD_90pct", sep = ""),] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                                                                               geneset = geneset, annotation_list = c("CADD13_PHRED"), annotation_cutoff = as.numeric(quantile(mutation$CADD_score[,1], c(0.9),na.rm = TRUE)),
                                                                                               CG_window_size = "CG_50bp", mut_type = "SNV", ref = ref, all_genes = "no", annotation = "yes")
    output[paste(deparse(substitute(geneset)),"+CADD_gt10", sep = ""),] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                                                                              geneset = geneset, annotation_list = c("CADD13_PHRED"), annotation_cutoff = 10,
                                                                                              CG_window_size = "CG_50bp", mut_type = "SNV", ref = ref, all_genes = "no", annotation = "yes")

    output['baseline+Eigen_95pct',] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                                          geneset = c("all_genes"), annotation_list = c("Eigen"), annotation_cutoff = as.numeric(quantile(mutation$mut_eigen, c(0.95),na.rm = TRUE)),
                                                          CG_window_size = "CG_50bp", mut_type = "SNV", ref = ref, all_genes = "yes", annotation = "yes")
    output['baseline+Eigen_90pct',] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                                          geneset = c("all_genes"), annotation_list = c("Eigen"), annotation_cutoff = as.numeric(quantile(mutation$mut_eigen, c(0.9),na.rm = TRUE)),
                                                          CG_window_size = "CG_50bp", mut_type = "SNV", ref = ref, all_genes = "yes", annotation = "yes")
    output[paste(deparse(substitute(geneset)),"+Eigen_95pct", sep = ""),] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                                                                                geneset = geneset, annotation_list = c("Eigen"), annotation_cutoff = as.numeric(quantile(mutation$mut_eigen, c(0.95),na.rm = TRUE)),
                                                                                                CG_window_size = "CG_50bp", mut_type = "SNV", ref = ref, all_genes = "no", annotation = "yes")
    output[paste(deparse(substitute(geneset)),"+Eigen_90pct", sep = ""),] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                                                                                geneset = geneset, annotation_list = c("Eigen"), annotation_cutoff = as.numeric(quantile(mutation$mut_eigen, c(0.9),na.rm = TRUE)),
                                                                                                CG_window_size = "CG_50bp", mut_type = "SNV", ref = ref, all_genes = "no", annotation = "yes")
    
    output['baseline+Phylop_95pct',] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                                           geneset = c("all_genes"), annotation_list = c("phyloP46wayAllElements"), annotation_cutoff = as.numeric(quantile(mutation$phylop[,6], c(0.95),na.rm = TRUE)),
                                                           CG_window_size = "CG_50bp", mut_type = "SNV", ref = ref, all_genes = "yes", annotation = "yes")
    output['baseline+Phylop_90pct',] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                                           geneset = c("all_genes"), annotation_list = c("phyloP46wayAllElements"), annotation_cutoff = as.numeric(quantile(mutation$phylop[,6], c(0.9),na.rm = TRUE)),
                                                           CG_window_size = "CG_50bp", mut_type = "SNV", ref = ref, all_genes = "yes", annotation = "yes")
    output[paste(deparse(substitute(geneset)),"+Phylop_95pct", sep = ""),] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                                                                                 geneset = geneset, annotation_list = c("phyloP46wayAllElements"), annotation_cutoff = as.numeric(quantile(mutation$phylop[,6], c(0.95),na.rm = TRUE)),
                                                                                                 CG_window_size = "CG_50bp", mut_type = "SNV", ref = ref, all_genes = "no", annotation = "yes")
    output[paste(deparse(substitute(geneset)),"+Phylop_90pct", sep = ""),] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                                                                                 geneset = geneset, annotation_list = c("phyloP46wayAllElements"), annotation_cutoff = as.numeric(quantile(mutation$phylop[,6], c(0.9),na.rm = TRUE)),
                                                                                                 CG_window_size = "CG_50bp", mut_type = "SNV", ref = ref, all_genes = "no", annotation = "yes")
    
    output['baseline+GERP_95pct',] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                                         geneset = c("all_genes"), annotation_list = c("GERP"), annotation_cutoff = as.numeric(quantile(mutation$gerp[,8], c(0.95),na.rm = TRUE)),
                                                         CG_window_size = "CG_50bp", mut_type = "SNV", ref = ref, all_genes = "yes", annotation = "yes")
    output['baseline+GERP_90pct',] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                                         geneset = c("all_genes"), annotation_list = c("GERP"), annotation_cutoff = as.numeric(quantile(mutation$gerp[,8], c(0.9),na.rm = TRUE)),
                                                         CG_window_size = "CG_50bp", mut_type = "SNV", ref = ref, all_genes = "yes", annotation = "yes")
    output['baseline+GERP_gt2',] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                                       geneset = c("all_genes"), annotation_list = c("GERP"), annotation_cutoff = 2,
                                                       CG_window_size = "CG_50bp", mut_type = "SNV", ref = ref, all_genes = "yes", annotation = "yes")
    output[paste(deparse(substitute(geneset)),"+GERP_95pct", sep = ""),] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                          geneset = geneset, annotation_list = c("GERP"), annotation_cutoff = as.numeric(quantile(mutation$gerp[,8], c(0.95),na.rm = TRUE)),
                          CG_window_size = "CG_50bp", mut_type = "SNV", ref = ref, all_genes = "no", annotation = "yes")
    output[paste(deparse(substitute(geneset)),"+GERP_90pct", sep = ""),] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                                                                               geneset = geneset, annotation_list = c("GERP"), annotation_cutoff = as.numeric(quantile(mutation$gerp[,8], c(0.9),na.rm = TRUE)),
                                                                                               CG_window_size = "CG_50bp", mut_type = "SNV", ref = ref, all_genes = "no", annotation = "yes")
    output[paste(deparse(substitute(geneset)),"+GERP_gt2", sep = ""),] = burden_CG_corrected(mutation_table = mutation$mutation, data_matrix = mutation$data_matrix, gene_mut_table = mut, CG_content_table = mutation$CG_content, 
                                                                                              geneset = geneset, annotation_list = c("GERP"), annotation_cutoff = 2,
                                                                                              CG_window_size = "CG_50bp", mut_type = "SNV", ref = ref, all_genes = "no", annotation = "yes")
  }
  else if (mut_type == "coding"){
    output = show_burden_for_coding_selected_features(mut, mutation$exon_mut_type)
  }
  output 
}

