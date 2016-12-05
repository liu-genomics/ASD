# code to simulate the mutation count for genes, and randomly assign windows to have a feature that will increase mutation rate for cases, if the gene is a risk gene
mutation_simulator <- function(gene_number = 20000, window_number=100, feature_prob = 0.2, base_rate = 2e-6 ,sample_size = 314, effect_size = 1.5, p_risk = 1){
  #[gene_number] the number of genes 
  #[window_number] the number of windows for each gene
  #[feature_prob] the probability of a window having a feature that is risk-related. 
  #[base_rate] the baseline rate of each window, assuming now each window has a fixed mutation rate
  #[sample_size] the number of individuals. 
  #[effect_size] is the increase of poisson rate due to having a relavant feature indicator
  #[p_risk] is the proportion of risk genes. default is 1 assuming all genes are risk genes
  base_mut = matrix(base_rate,gene_number, window_number)
  feature_indicator = matrix(rbinom(gene_number*window_number, 1, feature_prob), gene_number, window_number)
  risk_gene_index = sample(seq(1,gene_number),floor(gene_number*p_risk))
  risk_base_mut = exp(log(base_mut[risk_gene_index,])+log(effect_size)*risk_gene_index[risk_gene_index,])
  
  
}


gene_number = 20000;window_number=100;feature_prob = 0.2;base_rate = 2e-6;sample_size = 314
