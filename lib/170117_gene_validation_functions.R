# First load gene list definitions.
source("../lib/170117_gene_list_definition.R")




# function to count the number of novel genes and lost genes after incoporating non-coding mutations.
# The input is the output from bayes_factor_for_each_gene_v2
performance_table <- function(tada_table){
  snv_novel_q0.1 <- tada_table[tada_table$FDR_coding >= 0.1 & tada_table$FDR_all < 0.1,]
  snv_novel_q0.2 <- tada_table[tada_table$FDR_coding >= 0.2 & tada_table$FDR_all < 0.2,]
  snv_novel_q0.3 <- tada_table[tada_table$FDR_coding >= 0.3 & tada_table$FDR_all < 0.3,]
  snv_lost_q0.1 <- tada_table[tada_table$FDR_coding < 0.1 & tada_table$FDR_all >= 0.1,]
  snv_lost_q0.2 <- tada_table[tada_table$FDR_coding < 0.2 & tada_table$FDR_all >= 0.2,]
  snv_lost_q0.3 <- tada_table[tada_table$FDR_coding < 0.3 & tada_table$FDR_all >= 0.3,]
  snv_consistent_q0.1 <- tada_table[tada_table$FDR_coding < 0.1 & tada_table$FDR_all < 0.1,]
  snv_consistent_q0.2 <- tada_table[tada_table$FDR_coding < 0.2 & tada_table$FDR_all < 0.2,]
  snv_consistent_q0.3 <- tada_table[tada_table$FDR_coding < 0.3 & tada_table$FDR_all < 0.3,]
  snv_performance = matrix(c(dim(snv_novel_q0.1)[1],dim(snv_novel_q0.2)[1],dim(snv_novel_q0.3)[1],
                           dim(snv_consistent_q0.1)[1],dim(snv_consistent_q0.2)[1],dim(snv_consistent_q0.3)[1],
                           dim(snv_lost_q0.1)[1],dim(snv_lost_q0.2)[1],dim(snv_lost_q0.3)[1]),3,3)
rownames(snv_performance) = c("q0.1","q0.2","q0.3")
colnames(snv_performance) = c("novel", "consistent", "lost")
snv_performance
}

#### function to calculate the enrichment of genes in a given gene list
foe <- function(novel_gene_list,compare_gene_list, search_space){
#[search_space] is the total set of genes that have undergone TADA analysis. Used to be tada_coding_denovo$TadaName
  a = length(intersect(novel_gene_list, compare_gene_list))
  b = length(intersect(search_space,compare_gene_list))
  fc = (a/length(novel_gene_list))/(b/length(search_space))
  pvalue = phyper(a-1, b, length(search_space)-b,length(novel_gene_list),lower.tail = FALSE)
  c(fc, pvalue)
}

foe_table <-function(novel_gene_list, search_space){
  table = matrix(0,2,8)
  table[,1] = foe(novel_gene_list,union(Petrovski_haploinsufficiency_genelist[,1],Huang_Hpi_score[Huang_Hpi_score[,2] >= 0.95,1]),search_space)
  table[,2] = foe(novel_gene_list, FMRP_gene[,1],search_space)
  table[,3] = foe(novel_gene_list, Petrovski_RVIS_full_table[Petrovski_RVIS_full_table[,3]<25 & !is.na(Petrovski_RVIS_full_table[,3]),1],search_space)
  table[,4] = foe(novel_gene_list, gene_with_exp_mean[gene_with_exp_mean[,3]<0.25,1], search_space)
  table[,5] = foe(novel_gene_list, Brain_GO_gene[,1], search_space)
  table[,6] = foe(novel_gene_list, synaptomeDB_presynaptics[,1], search_space)
  table[,7] = foe(novel_gene_list, synaptomeDB_postsynaptics[,1], search_space)
  table[,8] = foe(novel_gene_list, rownames(ExAC_gene[ExAC_gene$z_lof_pct < 0.25,]), search_space)
  rownames(table) = c("Fold of enrichment", "pvalue")
  colnames(table) = c("HI", "FMRP", "RVIS top 25%", "Brainspan top 25%","Brain GO","presynaptic", "postsynaptic","ExAC top 25%")
  table
}


### get the functional charcteriziation of a gene list
get_function_table_for_each_gene_2 <- function(novel_gene_list){
  output = data.frame(genename = novel_gene_list)
  output = data.frame(output,
                      HI = as.numeric(is.element(output$genename, union(Petrovski_haploinsufficiency_genelist[,1],Huang_Hpi_score[Huang_Hpi_score[,2] >= 0.95,1]))),
                      FMRP = as.numeric(is.element(output$genename, FMRP_gene[,1])),
                      RVIS = NA,
                      ExAC_mis = NA,
                      ExAC_LoF = NA,
                      Brainspan = NA,
                      BrainGO = as.numeric(is.element(output$genename, Brain_GO_gene[,1])),
                      SFARI_strong_candidate = as.numeric(is.element(output$genename, SFRAI_gene_strong_candidate[,1])),
                      SFARI_high_confidence = as.numeric(is.element(output$genename, SFRAI_gene_high_confidence[,1])),
                      SFARI_suggestive = as.numeric(is.element(output$genename, SFRAI_gene_suggestive_evidence[,1])),
                      SFARI_syndromic = as.numeric(is.element(output$genename, SFRAI_gene_syndromic[,1])),
                      presynaptic = as.numeric(is.element(output$genename, synaptomeDB_presynaptics[,1])),
                      postsynaptic = as.numeric(is.element(output$genename, synaptomeDB_postsynaptics[,1])),
                      presynaptic_activezone = as.numeric(is.element(output$genename, synaptomeDB_presynaptics_activezone[,1])),
                      vesicles = as.numeric(is.element(output$genename, synaptomeDB_vesicles[,1]))
  )
  for(i in 1:length(output[,1])){
    if(is.element(output[i,]$genename, Petrovski_RVIS_full_table[,1])){
      output[i,]$RVIS = Petrovski_RVIS_full_table[Petrovski_RVIS_full_table[,1] == output[i,]$genename,3]
    }
  }
  for(i in 1:length(output[,1])){
    if(is.element(output[i,]$genename, rownames(ExAC_gene))){
      output[i,]$ExAC_mis = ExAC_gene[rownames(ExAC_gene) == output[i,]$genename,]$z_mis_pct*100
      output[i,]$ExAC_LoF = ExAC_gene[rownames(ExAC_gene) == output[i,]$genename,]$z_lof_pct*100
    }
  }
  for(i in 1:length(output[,1])){
    if(is.element(output[i,]$genename, gene_with_exp_mean[,1])){
      temp = gene_with_exp_mean[gene_with_exp_mean[,1] == as.character(output[i,]$genename),3]
      temp = mean(temp)
      output[i,]$Brainspan = temp*100
    }
  }
  output
}