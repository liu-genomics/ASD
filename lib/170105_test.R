# 3 binary and 1 continuous variables as features, use categorization to do optimization for mutation rates could be ajusted
adjust_mutation_rate_categarization <- function(sample_size_per_gene = 100, gene_prior_file = "../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.94.txt", relative_risk = c(1,1.5,2,0.3)){
  #[sample_size_per_gene] is the number of bases under consideraton for each gene, 15000 for each gene if the summation over all genes expected to be about 3e8
  #[gene_prior_file] = "../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.94.txt
  #[relative_risk] the increase of poisson rate if the gene is a risk gene
  # predictors are 3 binary variables and 1 continous variables
  gene_prior <- read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  gene_num <-nrow(gene_prior)
  sample_size <- gene_num*sample_size_per_gene
  predictors <- data.frame(ep1 <-rbinom(sample_size,1,rep(0.05,sample_size)),
                           ep2 <-rbinom(sample_size,1,rep(0.05,sample_size)),
                           ep3 <-rbinom(sample_size,1,rep(0.05,sample_size)),
                           annotation1 <- rnorm(sample_size)
                           
  )
  colnames(predictors) <- c("ep1","ep2","ep3","annotation1")
  predictors$genename <- rep(gene_prior[,1],each = sample_size_per_gene)
  
  gene_risk_assign <- mapply(function(x){rbinom(1,1,1-x)}, gene_prior[,2])
  gene_prior <- data.frame(gene_prior, risk_indicator = gene_risk_assign)
  predictors_temp_risk <- predictors[is.element(predictors$genename, gene_prior[gene_prior$risk_indicator == 1,1]),]
  predictors_temp_nonrisk <- predictors[!is.element(predictors$genename, gene_prior[gene_prior$risk_indicator == 1,1]),]
  
  active_effect <- as.matrix(predictors_temp_risk[,1:4])%*%relative_risk
  predictors_temp_risk <- data.frame(predictors_temp_risk, adjusted_mutrate = 1e-6*300, sim_mut_count = rpois(nrow(predictors_temp_risk),1e-6*300*exp(active_effect)))
  predictors_temp_nonrisk <- data.frame(predictors_temp_nonrisk, adjusted_mutrate = 1e-6*300, sim_mut_count = rpois(nrow(predictors_temp_nonrisk),1e-6*300))
  predictors <- rbind(predictors_temp_risk,predictors_temp_nonrisk)
  rm(predictors_temp_risk)
  rm(predictors_temp_nonrisk)
  
  if(optimization_prop !=1){
    genes_for_optimization <- gene_prior[order(gene_prior[,2],decreasing = FALSE),1]
    genes_for_optimization <- genes_for_optimization[1:floor(length(genes_for_optimization)*optimization_prop)]
    predictors <- predictors[is.element(predictors$genename, genes_for_optimization),]
  }
  #test_count = mapply(function(x){rpois(1,x)},mutrate_temp$mutrate2)
  fr<-function(x){ # the function that will be optimized later for promoter_beta, enhancer_beta
    # the order of relative risk in x is epigenomic_relative_risk followed by anno_relative_risk
    all_rr = x
    logP_Zg1 = by(predictors, predictors[,"genename"],  
                  function(x) sum(x$sim_mut_count*(log(x$adjusted_mutrate)+(as.matrix(x[,1:4])%*%all_rr))-x$adjusted_mutrate*exp((as.matrix(x[,1:4])%*%all_rr))-log(factorial(x$sim_mut_count))))
    logP_Zg0 = by(predictors, predictors[,"genename"],
                  function(x) sum(x$sim_mut_count*log(x$adjusted_mutrate)-x$adjusted_mutrate-log(factorial(x$sim_mut_count))))
    gene_prob_table = data.frame(genename = names(logP_Zg1), Zg1 = as.vector(logP_Zg1), Zg0 = as.vector(logP_Zg0))
    gene_prior = read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    gene_prior = data.frame(genename = gene_prior[,1], prob = 1-gene_prior[,2])
    gene_mle_table = merge(gene_prob_table, gene_prior, by.x = "genename", by.y = "genename")
    sum(log((gene_mle_table$prob*exp(gene_mle_table$Zg1)+(1-gene_mle_table$prob)*exp(gene_mle_table$Zg0)))) # minimization
  }
  run_time <- system.time(mle <- optim(rep(0.1,4), fr, control=list("fnscale"=-1), hessian = TRUE))
  
  fisher_info<-solve(-mle$hessian)
  prop_sigma<-sqrt(diag(fisher_info))
  upper<-mle$par+1.96*prop_sigma
  lower<-mle$par-1.96*prop_sigma
  
  list(mle <- mle, par_interval <- data.frame(upper_bound = upper, lower_bound = lower), run_time <- run_time)
}


# build a function to do categrization of a data.frame
categorization_data <- function(data){
  #[data] is a data.
} 

test <- predictors
test$annotation1 <- cut(test$annotation1, breaks = quantile(test$annotation1, seq(0,1,0.25)),labels = c("1","2","3","4"))
# better if levels of factors are predefined here. 
test.matrix <- model.matrix(sim_mut_count ~ ep1+ep2+ep3+annotation1-1, data = test)

system.time(test.split <- split(as.data.frame(test.matrix[1:100,]), test.matrix[,1:7]))


## a function to generate a table to do glm regression
sim_table <- function(sample_size = 2e8){
  # predictors are 5 binary variables and 2 continous variables
  predictors <- data.frame(ep1 <-rbinom(sample_size,1,rep(0.05,sample_size)),
                           ep2 <-rbinom(sample_size,1,rep(0.05,sample_size)),
                           ep3 <-rbinom(sample_size,1,rep(0.05,sample_size)),
                           ep4 <-rbinom(sample_size,1,rep(0.05,sample_size)),
                           ep5 <-rbinom(sample_size,1,rep(0.05,sample_size)),
                           annotation1 <- rnorm(sample_size),
                           annotation2 <- rnorm(sample_size)
  )
  colnames(predictors) <- c("ep1","ep2","ep3","ep4","ep5","annotation1","annotation2")
  predictors$count <- rpois(sample_size,exp(as.matrix(predictors)%*%c(2,0.5,0.7,1.4,0,0.2,0.5)))
  predictors
}

sample_size <- 1e6
data <- sim_table(sample_size)

running_time1 <- system.time(data.fit1 <- glm(count~ ep1 + ep2 + ep3 + ep4 + ep5 + annotation1 + annotation2, data = data, family = c("poisson")))

adjust_muatation_rate_optimization_test <- function(data){
  data$annotation1 <- cut(data$annotation1, breaks = quantile(data$annotation1, seq(0,1,0.2)),labels = c("1","2","3","4","5"), include.lowest = TRUE)
  data$annotation2 <- cut(data$annotation2, breaks = quantile(data$annotation2, seq(0,1,0.2)),labels = c("1","2","3","4","5"), include.lowest = TRUE)
  response <- data$count
  data <- model.matrix(count ~ ep1+ep2+ep3+ep4+ep5+annotation1+annotation2, data = data, contrasts.arg = lapply(data[,6:7], contrasts, contrasts=FALSE))[,-1]
  data <- data.frame(response, response = data)
  data.split <-split(as.data.frame(data[,-1]), as.data.frame(data[,-1]), drop = TRUE)
  
  fr<-function(x){ # the function that will be optimized later for promoter_beta, enhancer_beta
    # the order of relative risk in x is epigenomic_relative_risk followed by anno_relative_risk
    # annotation variables have been categorized, so the effect will take values like 1*beta, 2*beat, 3*beta, etc..,,
    all_rr <- c(x[1],x[2],x[3],x[4],x[5], 1*x[6],2*x[6],3*x[6],4*x[6],5*x[6],1*x[7],2*x[7],3*x[7],4*x[7],5*x[7])
    
    sum_cate <-function(data_subset){
      # assuming the base level rate is 1 for every base
      exp(as.matrix(data_subset)[1,]%*%all_rr)*nrow(data_subset)
    }
    # summation of part1 of Equation 4
    sum_logP_1 <- sum(data$response*(log(1)+as.matrix(data)[,-1]%*%all_rr))
    # summation of part2 of Equation 4
    sum_logP_2 <- sum(sapply(data.split, sum_cate))
    sum_logP_1 - sum_logP_2 
  }
  running_time2 <- system.time(mle <- optim(rep(0.1, 7), fr,control=list("fnscale"=-1, maxit = 5000), hessian = TRUE))  
  list(running_time = running_time2, mle = mle)
}

adjust_muatation_rate_optimization_test(data)



