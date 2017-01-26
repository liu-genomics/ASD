sim_table_simple <- function(sample_size = 2e8, epi_proportion = 0.1, base_mutrate = 1e-4, effect = 0.5){
  # only has one binary predictor
  # effect is on the log scale
  predictors <- data.frame(ep1 = rbinom(sample_size,1, epi_proportion))
  colnames(predictors) <- c("ep1")
  predictors$count <- rpois(sample_size,base_mutrate*exp(as.matrix(predictors)%*%c(effect)))
  predictors$base_mutrate <- base_mutrate
  predictors
}

# only process generated from sim_table_simple, use categorization to fit the parameter. 
adjust_muatation_rate_optimization_test_simple_categorization <- function(data){
  
  response <- data$count
  data_matrix <- model.matrix(count ~ ep1, data = data)
  data_matrix <- data.frame(response = data$count, data_matrix)
  split_time <- system.time(data.split <-split(as.data.frame(data_matrix[,-1]), as.data.frame(data_matrix[,-1]), drop = TRUE))
  data_with_response <- data_matrix[data_matrix$response !=0,]
  fr<-function(x){ # the function that will be optimized later for promoter_beta, enhancer_beta
    # the order of relative risk in x is epigenomic_relative_risk followed by anno_relative_risk
    # annotation variables have been categorized, so the effect will take values like 1*beta, 2*beat, 3*beta, etc..,,
    all_rr <- c(log(as.numeric(data$base_mutrate[1])),x)
    
    sum_cate <-function(data_subset){
      exp(as.matrix(data_subset)[1,]%*%all_rr)*nrow(data_subset)
    }
    # summation of part1 of Equation 4
    
    sum_logP_1 <- sum(data_with_response$response*(as.matrix(data_with_response)[,-1]%*%all_rr))
    # summation of part2 of Equation 4
    sum_logP_2 <- sum(sapply(data.split, sum_cate))
    sum_logP_1 - sum_logP_2 
  }
  optimization_time <- system.time(mle <- optim(c(0.1), method = "Brent", fr,lower = 0, upper = 10, control=list("fnscale"=-1, maxit = 5000), hessian = TRUE))  
  list(split_time = split_time,optimization_time = optimization_time, mle = mle)
}



# only process generated from sim_table_simple, use categorization to fit the parameter. 
# Version 2 could deal with situations where the base_mutrate is different for each base
# It is also faster than the 1st version because it didn't recalculate the sum of base_mutrate for each partition of dataset, everytime in every optimization iteration. 
adjust_muatation_rate_optimization_test_simple_categorization_v2 <- function(data){
  
  response <- data$count
  data_matrix <- model.matrix(count ~ ep1, data = data)
  data_matrix <- data.frame(response = data$count, base_mutrate = data$base_mutrate, data_matrix)
  split_time <- system.time(data.split <-split(as.data.frame(data_matrix[,-1]), as.data.frame(data_matrix[,-1:-2]), drop = TRUE))
  data_with_response <- data_matrix[data_matrix$response !=0,]
  mutrate_sum_per_group <- sapply(data.split, function(x) sum(x$base_mutrate))
  fr<-function(x){ # the function that will be optimized later for promoter_beta, enhancer_beta
    # the order of relative risk in x is epigenomic_relative_risk followed by anno_relative_risk
    # annotation variables have been categorized, so the effect will take values like 1*beta, 2*beat, 3*beta, etc..,,
    all_rr <- c(log(data$base_mutrate[1]),x)
    
    sum_cate <-function(data_subset, mutrate_sum){
      exp(as.matrix(data_subset)[1,3]*x)*mutrate_sum
    }
    # summation of part1 of Equation 4
    sum_logP_1 <- sum(data_with_response$response*(log(data_with_response$base_mutrate)+as.matrix(data_with_response)[,4]*x))
    # summation of part2 of Equation 4
    sum_logP_2 <- sum(mapply(sum_cate, data.split, mutrate_sum_per_group))
    sum_logP_1 - sum_logP_2 
  }
  optimization_time <- system.time(mle <- optim(c(0.1), method = "Brent", fr,lower = 0, upper = 10, control=list("fnscale"=-1, maxit = 5000), hessian = TRUE))    
  list(split_time = split_time,optimization_time = optimization_time, mle = mle)
}

# only process generated from sim_table_simple, use categorization to fit the parameter. 
# Version 2 could deal with situations where the base_mutrate is different for each base
# It is also faster than the 1st version because it didn't recalculate the sum of base_mutrate for each partition of dataset, everytime in every optimization iteration. 
adjust_muatation_rate_optimization_test_simple_categorization_v2 <- function(data){
  
  response <- data$count
  data_matrix <- model.matrix(count ~ ep1, data = data)
  data_matrix <- data.frame(response = data$count, base_mutrate = data$base_mutrate, data_matrix)
  split_time <- system.time(data.split <-split(as.data.frame(data_matrix[,-1]), as.data.frame(data_matrix[,-1:-2]), drop = TRUE))
  data_with_response <- data_matrix[data_matrix$response !=0,]
  mutrate_sum_per_group <- sapply(data.split, function(x) sum(x$base_mutrate))
  fr<-function(x){ # the function that will be optimized later for promoter_beta, enhancer_beta
    # the order of relative risk in x is epigenomic_relative_risk followed by anno_relative_risk
    # annotation variables have been categorized, so the effect will take values like 1*beta, 2*beat, 3*beta, etc..,,
    all_rr <- c(log(data$base_mutrate[1]),x)
    
    sum_cate <-function(data_subset, mutrate_sum){
      exp(as.matrix(data_subset)[1,3]*x)*mutrate_sum
    }
    # summation of part1 of Equation 4
    sum_logP_1 <- sum(data_with_response$response*(log(data_with_response$base_mutrate)+as.matrix(data_with_response)[,4]*x))
    # summation of part2 of Equation 4
    sum_logP_2 <- sum(mapply(sum_cate, data.split, mutrate_sum_per_group))
    sum_logP_1 - sum_logP_2 
  }
  optimization_time <- system.time(mle <- optim(c(0.1), method = "Brent", fr,lower = 0, upper = 10, control=list("fnscale"=-1, maxit = 5000), hessian = TRUE))    
  list(split_time = split_time,optimization_time = optimization_time, mle = mle)
}

# only process generated from sim_table_simple, use categorization to fit the parameter. 
# Version 2 could deal with situations where the base_mutrate is different for each base
# It is also faster than the 1st version because it didn't recalculate the sum of base_mutrate for each partition of dataset, everytime in every optimization iteration. 
adjust_muatation_rate_optimization_test_simple_categorization_v2_temp <- function(data){
  
  response <- data$count
  data_matrix <- model.matrix(count ~ ep1, data = data)
  data_matrix <- data.frame(response = data$count, base_mutrate = data$base_mutrate, data_matrix)
  split_time <- system.time(data.split <-split(as.data.frame(data_matrix[,-1]), as.data.frame(data_matrix[,-1:-2]), drop = TRUE))
  data_with_response <- data_matrix[data_matrix$response !=0,]
  mutrate_sum_per_group <- sapply(data.split, function(x) sum(x$base_mutrate))
  fr<-function(x){ # the function that will be optimized later for promoter_beta, enhancer_beta
    # the order of relative risk in x is epigenomic_relative_risk followed by anno_relative_risk
    # annotation variables have been categorized, so the effect will take values like 1*beta, 2*beat, 3*beta, etc..,,
    all_rr <- c(log(data$base_mutrate[1]),x)
    
    sum_cate_1 <-function(data_subset){
      exp(as.matrix(data_subset)[1,3]*x)*mutrate_sum_per_group[[1]]
    }
    sum_cate_2 <-function(data_subset){
      exp(as.matrix(data_subset)[1,3]*x)*mutrate_sum_per_group[[2]]
    }
    # summation of part1 of Equation 4
    sum_logP_1 <- sum(data_with_response$response*(log(data_with_response$base_mutrate)+as.matrix(data_with_response)[,4]*x))
    # summation of part2 of Equation 4
    sum_logP_2_1 <- sum(sapply(data.split[1], sum_cate_1))
    sum_logP_2_2 <- sum(sapply(data.split[2], sum_cate_2))
    
    sum_logP_2 <- sum_logP_2_1 + sum_logP_2_2
    sum_logP_1 - sum_logP_2 
  }
  optimization_time <- system.time(mle <- optim(c(0.1), method = "Brent", fr,lower = 0, upper = 10, control=list("fnscale"=-1, maxit = 5000), hessian = TRUE))    
  list(split_time = split_time,optimization_time = optimization_time, mle = mle)
}



# only process generated from sim_table_simple, use categorization to fit the parameter. 
# Version 3 could deal with situations where the base_mutrate is different for each base
# Compared to V2, V3 tries to change the code in a way that avoids using mapply. which turns out to be very slow. 
adjust_muatation_rate_optimization_test_simple_categorization_v3 <- function(data){
  
  response <- data$count
  data_matrix <- model.matrix(count ~ ep1, data = data)
  data_matrix <- data.frame(response = data$count, base_mutrate = data$base_mutrate, data_matrix)
  split_time <- system.time(data.split <-split(as.data.frame(data_matrix[,-1]), as.data.frame(data_matrix[,-1:-2]), drop = TRUE))
  data_with_response <- data_matrix[data_matrix$response !=0,]
  
  partition_simplify <- function(x){
    c(x[1,3],sum(x[,1]))
  }
  
  simplified_partition <- sapply(data.split,partition_simplify, simplify = FALSE)
  
  
  fr<-function(x){ # the function that will be optimized later for promoter_beta, enhancer_beta
    # the order of relative risk in x is epigenomic_relative_risk followed by anno_relative_risk
    # annotation variables have been categorized, so the effect will take values like 1*beta, 2*beat, 3*beta, etc..,,
    all_rr <- c(log(data$base_mutrate[1]),x)
    
    sum_cate <-function(data_subset){
      exp(data_subset[1]*x)*data_subset[2]
    }
    
    # summation of part1 of Equation 4
    sum_logP_1 <- sum(data_with_response$response*(log(data_with_response$base_mutrate)+as.matrix(data_with_response)[,4]*x))
    # summation of part2 of Equation 4
    sum_logP_2 <- sum(sapply(simplified_partition, sum_cate))
    sum_logP_1 - sum_logP_2 
  }
  optimization_time <- system.time(mle <- optim(c(0.1), method = "Brent", fr,lower = 0, upper = 10, control=list("fnscale"=-1, maxit = 5000), hessian = TRUE))    
  list(split_time = split_time,optimization_time = optimization_time, mle = mle)
}


#### function to study the optimization time without partition
adjust_muatation_rate_optimization_test_simple_without_categorization <- function(data){
  response <- data$count
  data_matrix <- model.matrix(count ~ ep1, data = data)
  data_matrix <- data.frame(response = data$count, base_mutrate = data$base_mutrate, data_matrix)
  
  fr<-function(x){ # the function that will be optimized
    sum_logP <- sum(data_matrix$response*(log(data_matrix$base_mutrate)+as.matrix(data_matrix)[,4]*x))-sum(data_matrix$base_mutrate*exp(as.matrix(data_matrix)[,4]*x))
  }
  optimization_time <- system.time(mle <- optim(c(0.1), method = "Brent", fr,lower = 0, upper = 10, control=list("fnscale"=-1, maxit = 5000), hessian = TRUE))    
  list(optimization_time = optimization_time, mle = mle)
}


# only process generated from sim_table_simple, use glm to fit the parameters
adjust_muatation_rate_optimization_test_simple_glm <- function(data){
  optimization_time <- system.time(out.offset <- glm(count ~ ep1+offset(log(data$base_mutrate)), family = poisson, data = data))

  list(optimization_time = optimization_time, mle = out.offset$coefficients)  
}


#### functions to simulate genes and mutations at base level. 
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
    mut_number = rbind(risk_mut_number, nonrisk_mut_number)
  }
  
  mut_number = melt(mut_number)
  base_mut = melt(base_mut)
  feature_indicator = melt(feature_indicator)
  data = merge(mut_number, base_mut, by.x = c("X1","X2"), by.y = c("X1","X2"))
  data = merge(data, feature_indicator, by.x = c("X1","X2"), by.y = c("X1","X2"))
  colnames(data) = c("gene","window","mut_count","base_mut_rate","feature_indicator")
  
  list(data = data, sample_size = sample_size)
}

# mutation_simulator_v2 will generate two features.
# The purpose is to make the estimate_effect_size_for_simulation_data_mixture_poisson more generate to deal with more than 1 features.
# for simplicity these two features will have the same relative risk and proportion
mutation_simulator_v2 <- function(gene_number = 200, window_number=100, feature_prob = 0.2, base_rate = 2e-6 ,sample_size = 314, effect_size = 1.5, p_risk = 1){
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
  feature_indicator_2 = matrix(rbinom(gene_number*window_number, 1, feature_prob), gene_number, window_number)
  risk_gene_index = sample(seq(1,gene_number),floor(gene_number*p_risk))
  risk_base_mut = exp(log(base_mut[risk_gene_index,])+log(effect_size)*feature_indicator[risk_gene_index,]+log(effect_size)*feature_indicator_2[risk_gene_index,])
  if(p_risk == 1){
    mut_number = apply(risk_base_mut, c(1,2), function(x){rpois(1,2*sample_size*x)})
  }
  else {
    nonrisk_base_mut = base_mut[!is.element(rownames(base_mut), risk_gene_index),]
    risk_mut_number = apply(risk_base_mut, c(1,2), function(x){rpois(1,2*sample_size*x)})
    nonrisk_mut_number = apply(nonrisk_base_mut, c(1,2), function(x){rpois(1,2*sample_size*x)})
    mut_number = rbind(risk_mut_number, nonrisk_mut_number)
  }
  
  mut_number = melt(mut_number)
  base_mut = melt(base_mut)
  feature_indicator = melt(feature_indicator)
  feature_indicator_2 = melt(feature_indicator_2)
  
  data = merge(mut_number, base_mut, by.x = c("X1","X2"), by.y = c("X1","X2"))
  data = merge(data, feature_indicator, by.x = c("X1","X2"), by.y = c("X1","X2"))
  data = merge(data, feature_indicator_2, by.x = c("X1","X2"), by.y = c("X1","X2"))
  colnames(data) = c("gene","window","mut_count","base_mut_rate","feature_indicator", "feature_indicator_2")
  
  list(data = data, sample_size = sample_size)
}


estimate_effect_size_for_simulation_data_mixture_with_categorization<-function(data,feature_start = 5, feature_end = 6, feature_number = 2, prior_prob = 0.3){
  #[data] is the returned value from function{mutation_simulator}
  #[feature_start] is the column number of the first feature in data$data
  #[feature_end] is the column number of the last feature in data$data
  #[feature_number] is the number of features
  #[prior_prob] is the prior probability of each gene being a risk gene. For simplicity, set to be equal for all genes, and is eual to the prior probability used in mutation_simulator_v2
  
  # first get the real mutation rate by multiplying 2 and sample size.
  data$data$base_mut_rate = data$data$base_mut_rate * data$sample_size *2
  
  partition_time1 <- system.time(partition_by_gene <- split(data$data, data$data$gene))
  # function to get effective information of each element of partition_by_gene
  # These information are those necessary to compute log-likelihood in the optimization function
  partition_feature <- function(pbg){
    # input is one element of the list of partition_by_gene
    pbg_split <- split(pbg, pbg[,feature_start:feature_end],drop = TRUE)
    info_for_each_feature <- function(feature_set){
      list(feature_vector = as.numeric(feature_set[1,feature_start:feature_end]), sum_mut_rate_count = sum(feature_set[,3]*log(feature_set[,4])), sum_mut_rate = sum(feature_set[,4]), sum_mut_count = sum(feature_set[,3]), log_fcount = sum(log(factorial(feature_set[,3]))))
    }
    sapply(pbg_split, info_for_each_feature,simplify = FALSE)
  }
  partition_time2 <- system.time(data_partition <- sapply(partition_by_gene, partition_feature, simplify = FALSE))
  
  fr<-function(x){
    all_rr = x
    cal_logP_Zg1 <- function(data_partition_element){
      cal_logP_Zg1_level2 <-function(data_partition_element_level2){
        data_partition_element_level2[[2]]+data_partition_element_level2[[4]]*data_partition_element_level2[[1]]%*%all_rr-data_partition_element_level2[[3]]*exp(data_partition_element_level2[[1]]%*%all_rr)-data_partition_element_level2[[5]]
      }
      sum(sapply(data_partition_element, cal_logP_Zg1_level2))
    }
    
    cal_logP_Zg0 <- function(data_partition_element){
      cal_logP_Zg0_level2 <-function(data_partition_element_level2){
        data_partition_element_level2[[2]]-data_partition_element_level2[[3]]-data_partition_element_level2[[5]]
      }
      sum(sapply(data_partition_element, cal_logP_Zg0_level2))
    }
    
    logP_Zg1 = sapply(data_partition, cal_logP_Zg1)
    logP_Zg0 = sapply(data_partition, cal_logP_Zg0)
    sum(log((prior_prob*exp(logP_Zg1)+(1-prior_prob)*exp(logP_Zg0))))
  }
  optimization_time <- system.time(mle <- optim(rep(0.1, feature_number), fr,control=list("fnscale"=-1), hessian = TRUE))
  list(partition_time1 = partition_time1, partition_time2 = partition_time2, optimization_time = optimization_time,mle = mle)
}

# This function estimates the relative risk from simulated data generated by mutation_simulator_v2

estimate_effect_size_for_simulation_data_mixture_without_categorization<-function(data, feature_start = 5, feature_end = 6, feature_number = 2, prior_prob = 0.3){
  #[data] is the returned value from function{mutation_simulator}
  #[feature_start] is the column number of the first feature in data$data
  #[feature_end] is the column number of the last feature in data$data
  #[feature_number] is the number of features
  #[prior_prob] is the prior probability of each gene being a risk gene. For simplicity, set to be equal for all genes, and is eual to the prior probability used in mutation_simulator_v2
  # first get the real mutation rate by multiplying 2 and sample size.
  data$data$base_mut_rate = data$data$base_mut_rate * data$sample_size *2
  fr<-function(x){
    all_rr <- x
    logP_Zg1 <- by(data$data, data$data[,"gene"],  
                   function(x) sum(x$mut_count*(log(x$base_mut_rate)+(as.matrix(x[,feature_start:feature_end])%*%all_rr))-x$base_mut_rate*exp((as.matrix(x[,feature_start:feature_end])%*%all_rr))-log(factorial(x$mut_count))))
    logP_Zg0 <- by(data$data, data$data[,"gene"],
                   function(x) sum(x$mut_count*log(x$base_mut_rate)-x$base_mut_rate-log(factorial(x$mut_count))))
    sum(log((prior_prob*exp(logP_Zg1)+(1-prior_prob)*exp(logP_Zg0)))) # minimization
  }
  optimization_time <- system.time(mle <- optim(rep(0.1, feature_number), fr,control=list("fnscale"=-1), hessian = TRUE))
  list(optimization_time = optimization_time,mle = mle)
}