# code to simulate the mutation count for genes, and randomly assign windows to have a feature that will increase mutation rate for cases, if the gene is a risk gene.
# The output will be in long format. within each window is a line, and each line has information of adjusted mutation rate, and if it is in a risk gene, how many mutations are there
library(reshape)
mutation_simulator <- function(gene_number = 20000, window_number=100, feature_prob = 0.2, base_rate = 2e-6 ,sample_size = 314, effect_size = 1.5, p_risk = 1){
  #[gene_number] the number of genes 
  #[window_number] the number of windows for each gene
  #[feature_prob] the probability of a window having a feature that is risk-related. 
  #[base_rate] the baseline rate of each window, assuming now each window has a fixed mutation rate
  #[sample_size] the number of individuals. 
  #[effect_size] is the increase of poisson rate due to having a relavant feature indicator
  #[p_risk] is the proportion of risk genes. default is 1 assuming all genes are risk genes
  base_mut = matrix(base_rate,gene_number, window_number)
  rownames(base_mut) = seq(1,nrow(base_mut))
  feature_indicator = matrix(rbinom(gene_number*window_number, 1, feature_prob), gene_number, window_number)
  risk_gene_index = sample(seq(1,gene_number),floor(gene_number*p_risk))
  risk_base_mut = exp(log(base_mut[risk_gene_index,])+log(effect_size)*feature_indicator[risk_gene_index,])
  if(p_risk == 1){
    mut_number = apply(risk_base_mut, c(1,2), function(x){rpois(1,2*sample_size*x)})
  }
  else {
    nonrisk_base_mut = base_mut[!is.element(rownames(base_mut), risk_gene_index),]
    risk_mut_number = apply(risk_base_mut, c(1,2), function(x){rpois(1,2*sample_size*x)})
    nonrisk_mut_number = apply(nonrisk_base_mut, c(1,2), function(x){rpois(1,2*sample_size*x)})
    mut_number = rbind(risk_base_mut, nonrisk_mut_number)
  }
  
  mut_number = melt(mut_number)
  base_mut = melt(base_mut)
  feature_indicator = melt(feature_indicator)
  data = merge(mut_number, base_mut, by.x = c("X1","X2"), by.y = c("X1","X2"))
  data = merge(data, feature_indicator, by.x = c("X1","X2"), by.y = c("X1","X2"))
  colnames(data) = c("gene","window","mut_count","base_mut_rate","feature_indicator")
  
  list(data = data, sample_size = sample_size)
}



mutation_count = mutation_simulator(gene_number = 20000,window_number=100,feature_prob = 0.2,base_rate = 2e-6,sample_size = 314, p_risk = 1, effect_size = 1.5)

estimate_effect_size_for_simulation_data<-function(data){
  #[data] is the returned value from function{mutation_simulator}
  fr<-function(x){
    effect = x
    logP = by(data$data,data$data$gene, function(x) sum(data$data$mut_count*(log(data$data$base_mut_rate)+data$data$feature_indicator*effect+log(data$sample_size))-data$data$base_mut_rate*exp(data$data$feature_indicator*effect)*data$sample_size-log(factorial(data$data$mut_count))))
    sum(logP)
  }
  mle = optim(1, fr)
  mle
}
estimate = estimate_effect_size_for_simulation_data(mutation_count)
save.image("../analysis/161207_mutation_simulation.Rdata")
#gene_number = 20000;window_number=100;feature_prob = 0.2;base_rate = 2e-6;sample_size = 314; p_risk = 1;effect_size = 1.5
