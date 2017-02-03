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

# version 2 allows for different prior_prob for each gene
estimate_effect_size_for_simulation_data_mixture_with_categorization_v2<-function(data,feature_start = 6, feature_end = 6, feature_number = 1, gene_prior_file, sample_size = 314){
  #[data] a data.table
  #[feature_start] is the column number of the first feature in data$data
  #[feature_end] is the column number of the last feature in data$data
  #[feature_number] is the number of features
  #[gene_prior_file], a file that has prior (derived from posterior and prior)for a gene as a risk gene. genename \t probability e.g. ../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.94.txt
  
  # first get the real mutation rate by multiplying 2 and sample size.
  data$adjusted_base_mutrate <- data$adjusted_base_mutrate * sample_size *2
  
  # remove bases that don't have mutation rates in the reference at the base level (the output would be 0 from bigwigaverageoverbed). 
  data <- data[adjusted_base_mutrate != 0]
  
  # get the piror probability of genes.
  gene_prior = fread(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_prior) = c("genename", "prior")
  # change prior to the 
  gene_prior$prior = 1- gene_prior$prior
  
  data <- gene_prior[data, on = "genename"]
  
  # remove genes that don't have a prior probability
  data <- data[complete.cases(data)]
  
  partition_time1 <- system.time(partition_by_gene <- split(data, data$genename))
  # function to get effective information of each element of partition_by_gene
  # These information are those necessary to compute log-likelihood in the optimization function
  partition_feature <- function(pbg){
    # input is one element of the list of partition_by_gene
    pbg_split <- split(pbg, pbg[,feature_start:feature_end],drop = TRUE)
    feature_combination_number <- length(pbg_split)
    info_for_each_feature <- function(feature_set){
      list(feature_vector = c(as.numeric(feature_set[1,feature_start:feature_end])), sum_mut_rate_count = sum(feature_set$mut_count*log(feature_set$adjusted_base_mutrate)), sum_mut_rate = sum(feature_set$adjusted_base_mutrate), sum_mut_count = sum(feature_set$mut_count), log_fcount = sum(log(factorial(feature_set$mut_count))), Zg1_ave_log_prior = log(feature_set[1,]$prior)/feature_combination_number, Zg0_ave_log_prior = log(1-feature_set[1,]$prior)/feature_combination_number)
    }
    sapply(pbg_split, info_for_each_feature,simplify = FALSE)
  }
  partition_time2 <- system.time(data_partition <- sapply(partition_by_gene, partition_feature, simplify = FALSE))
  
  fr<-function(x){
    all_rr = x
    cal_logP_Zg1 <- function(data_partition_element){
      cal_logP_Zg1_level2 <-function(data_partition_element_level2){
        data_partition_element_level2[[6]]+data_partition_element_level2[[2]]+data_partition_element_level2[[4]]*data_partition_element_level2[[1]]%*%all_rr-data_partition_element_level2[[3]]*exp(data_partition_element_level2[[1]]%*%all_rr)-data_partition_element_level2[[5]]
      }
      sum(sapply(data_partition_element, cal_logP_Zg1_level2))
    }
    
    cal_logP_Zg0 <- function(data_partition_element){
      cal_logP_Zg0_level2 <-function(data_partition_element_level2){
        data_partition_element_level2[[7]]+data_partition_element_level2[[2]]-data_partition_element_level2[[3]]-data_partition_element_level2[[5]]
      }
      sum(sapply(data_partition_element, cal_logP_Zg0_level2))
    }
    
    logP_Zg1 = sapply(data_partition, cal_logP_Zg1)
    logP_Zg0 = sapply(data_partition, cal_logP_Zg0)
    sum(log(exp(logP_Zg1)+exp(logP_Zg0)))
  }
  if (feature_number == 1){
    optimization_time <- system.time(mle <- optim(rep(0.1, feature_number), fr,method = "Brent", lower = -1, upper = 5,control=list("fnscale"=-1), hessian = TRUE))
  }
  else{
    optimization_time <- system.time(mle <- optim(rep(0.1, feature_number), fr,control=list("fnscale"=-1), hessian = TRUE))
  }
  list(partition_time1 = partition_time1, partition_time2 = partition_time2, optimization_time = optimization_time,mle = mle)
}




# This version allows for dealing with really big data, e.g., > 2 billion bases
# will read in [data_file] by chunk. 
estimate_effect_size_for_simulation_data_mixture_with_categorization_for_big_data<-function(data_file,feature_start = 6, feature_end = 6, feature_number = 1, gene_prior_file, sample_size = 314){
  #[data_file] is the file that has base level information (e.g., 170130_base_level_two_feature_mixture_poisson_rr_estimate_all_gene_test_in_RCC.info)
  #[feature_start] is the column number of the first feature in data$data
  #[feature_end] is the column number of the last feature in data$data
  #[feature_number] is the number of features
  #[gene_prior_file], a file that has prior (derived from posterior and prior)for a gene as a risk gene. genename \t probability e.g. ../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.94.txt
  
  
  base_num <- system(paste("wc -l ", data_file, sep = ""), intern = TRUE)
  base_num <- as.numeric(strsplit(base_num, split = " ")[[1]][1])
  # read in data in 10 chunks one by one.
  chunk_num <- 11
  base_per_chunk <- ceiling(base_num/chunk_num)
  
  # get the piror probability of genes.
  gene_prior = fread(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_prior) = c("genename", "prior")
  # change prior to the 
  gene_prior$prior = 1- gene_prior$prior
  gene_prior$genename = as.character(gene_prior$genename)
  
  data_col_names <- c()
  partition_time1 <-c()
  partition_time2 <-c()
  data_partition <-list()
  for(i in 1:chunk_num){
    
    if(i == 1){
      data <-fread(data_file, skip = (i-1)*base_per_chunk, nrows = base_per_chunk, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      data_col_names <- colnames(data)
    }
    else{
      data <-fread(data_file, skip = (i-1)*base_per_chunk, nrows = base_per_chunk, header = FALSE, sep = "\t", stringsAsFactors = FALSE) 
      colnames(data) <- data_col_names
    }
    data<- data[,c("genename","base_ID","mut_count","adjusted_base_mutrate","epi")]
    data$genename <- as.character(data$genename)
    # remove bases that don't have mutation rates in the reference at the base level (the output would be 0 from bigwigaverageoverbed). 
    data <- data[adjusted_base_mutrate != 0]
    
    # first get the real mutation rate by multiplying 2 and sample size.
    data$adjusted_base_mutrate <- data$adjusted_base_mutrate * sample_size *2
    
    data <- gene_prior[data, on = "genename"]
    # remove genes that don't have a prior probability
    data <- data[complete.cases(data)]
    partition_time1[i] <- system.time(partition_by_gene <- split(data, data$genename))
    # function to get effective information of each element of partition_by_gene
    # These information are those necessary to compute log-likelihood in the optimization function
    partition_feature <- function(pbg){
      # input is one element of the list of partition_by_gene
      pbg_split <- split(pbg, pbg[,feature_start:feature_end],drop = TRUE)
      feature_combination_number <- length(pbg_split)
      # this function below is different from the function used in dealing with dataset without reading by chunk. Here, prior is not incoporated at this step.
      info_for_each_feature <- function(feature_set){
        list(feature_vector = c(as.numeric(feature_set[1,feature_start:feature_end])), sum_mut_rate_count = sum(feature_set$mut_count*log(feature_set$adjusted_base_mutrate)), sum_mut_rate = sum(feature_set$adjusted_base_mutrate), sum_mut_count = sum(feature_set$mut_count), log_fcount = sum(log(factorial(feature_set$mut_count))))
      }
      sapply(pbg_split, info_for_each_feature,simplify = FALSE)
    }
    partition_time2[i] <- system.time(data_partition_per_chunk <- sapply(partition_by_gene, partition_feature, simplify = FALSE))
    data_partition <- append(data_partition, data_partition_per_chunk)
  }
  
  # notice the fr function is different from the function that deals with dataset without partition
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
    
    logP_table<-data.table(logP_Zg1 = logP_Zg1, logP_Zg0 = logP_Zg0, genename = names(logP_Zg1))
    
    logP_table <- logP_table[gene_prior, on = "genename"]
    logP_table <- logP_table[complete.cases(logP_table)]
    
    sum(by(logP_table, logP_table$genename, function(x){log(exp(log(x[1,]$prior)+sum(x$logP_Zg1))+exp(log(1-x[1,]$prior)+sum(x$logP_Zg0)))}))
    
  }
  if (feature_number == 1){
    optimization_time <- system.time(mle <- optim(rep(0.1, feature_number), fr,method = "Brent", lower = -1, upper = 5,control=list("fnscale"=-1), hessian = TRUE))
  }
  else{
    optimization_time <- system.time(mle <- optim(rep(0.1, feature_number), fr,control=list("fnscale"=-1), hessian = TRUE))
  }
  list(partition_time1 = partition_time1, partition_time2 = partition_time2, optimization_time = optimization_time,mle = mle)
}


# This version allows for dealing with really big data, e.g., > 2 billion bases
# will read in [data_file] by chunk. 
estimate_effect_size_for_simulation_data_mixture_with_categorization_for_compact_data<-function(data,feature_start = 6, feature_end = 6, feature_number = 1, gene_prior_file, sample_size = 314){
  #[data] is the object returned from [adjust_mutation_rate_window_to_base_compact], which has compact base-level information. 
  #[feature_start] is the column number of the first feature in data$data
  #[feature_end] is the column number of the last feature in data$data
  #[feature_number] is the number of features
  #[gene_prior_file], a file that has prior (derived from posterior and prior)for a gene as a risk gene. genename \t probability e.g. ../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.94.txt
  
  
  # get the piror probability of genes.
  gene_prior = fread(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_prior) = c("genename", "prior")
  # change prior to the 
  gene_prior$prior = 1- gene_prior$prior
  gene_prior$genename = as.character(gene_prior$genename)
  # notice the fr function is different from the function that deals with dataset without partition
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
    
    logP_Zg1 = sapply(data, cal_logP_Zg1)
    logP_Zg0 = sapply(data, cal_logP_Zg0)
    
    logP_table<-data.table(logP_Zg1 = logP_Zg1, logP_Zg0 = logP_Zg0, genename = names(logP_Zg1))
    
    logP_table <- logP_table[gene_prior, on = "genename"]
    logP_table <- logP_table[complete.cases(logP_table)]
    
    sum(by(logP_table, logP_table$genename, function(x){log(exp(log(x[1,]$prior)+sum(x$logP_Zg1))+exp(log(1-x[1,]$prior)+sum(x$logP_Zg0)))}))
    
  }
  if (feature_number == 1){
    optimization_time <- system.time(mle <- optim(rep(0.1, feature_number), fr,method = "Brent", lower = -1, upper = 5,control=list("fnscale"=-1), hessian = TRUE))
  }
  else{
    optimization_time <- system.time(mle <- optim(rep(0.1, feature_number), fr,control=list("fnscale"=-1), hessian = TRUE))
  }
  list(optimization_time = optimization_time, mle = mle)
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

library(data.table)
library(parallel)
# This is a function to run glm for 50bp windows in order to adjust for batch-effects on observed mutation rates, then will output a Rmd object that has base level mutation rates 
adjust_mutation_rate_window_to_base <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks = "no", overlap = 0.5, rm_nonsyn = FALSE, annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", annovar_input = "no.txt", 
                                                gene_assign_file, gene_prior_file, report_proportion,  mutrate_ref_file, node_n =6){
  # [mut_file] is the mutation location file from one dataset. e.g.,../analysis/160229_data_for_analysis.bed, chr \t start \t end \t index,  1-based start and end
  # [window_file] is the file with genome partitioned into windows. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed
  # [cg_file] is the file that has CG-content calculated, matching [window_file]. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.cg
  # [mutrate_file] is the file that has mutrate data. e.g., "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate"
  # [sample_size] is the number of individuals
  # [epitenomic_marks] if it is "no", then glm won't use if a mutation is in h3k27ac as a categorical variable, if it is a epigenomic annotation name.
  # for example, "../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed". It is fine to directly use original epigenomics bed files
  # instead of the partitioned files. Because the window_file already has annotations for each window in terms of whetehr they are promoter or coding. 
  # [overlap] the fraction of overlap (the paramter f in bedtools intersect), default is 0.5, For the ease of TADA modeling, it is better to assign it to be 1E-9(1bp)
  # [rm_nonsyn] TRUE or FALSE. Indicating whether nonysnonymous mutations should be removed when fitting the model. Use TRUE when fitting model using ASD mutations. 
  # [annovar_folder] is where the annovar program is, e.g., . e.g., /media/yuwen/Elements/ANNOVAR/annovar/
  # [annovar_input] is a input file generated from mutations that will be annotated by annovar. The format is chr \t start \t end \t ref_allele \t alt_allele \t index 1-based for starts and ends
  # [gene_assign_file] is the file that assigns each window to a gene. for example ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename.bed
  # [gene_prior_file], a file that has prior (derived from posterior and prior)for a gene as a risk gene. genename \t probability e.g. ../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.94.txt
  # [report_proportion] Choose the top X% TADA genes to report. 
  # [node_n] is the number of nodes used to run , default is 6
  # [mutrate_ref_file] the file with base level mutation rate reference. e.g., "/media/yuwen/Elements/mutation_rate_raw/MarkDaly/Daly_mutrate.bw"
  prefix=system("date +%s", intern = TRUE) # prefix for temporary files that will be deleted at the end of the pipeline
  
  mut = read.delim(mut_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  if(rm_nonsyn){
    # command to transform the start base to 0-based in order to use bedtools to do overlap
    command = paste("awk {'print $1\"\t\"$2-1\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$7'} ", annovar_input, " > ", paste(prefix, "_annovar_for_bedtools.bed",sep = ""),sep = "")
    system(command)
    command = paste("grep coding ", window_file, " > ", paste(prefix, "_coding_window_temp.bed", sep = ""), sep = "")
    system(command)
    command = paste("bedtools intersect -a ",paste(prefix, "_annovar_for_bedtools.bed",sep = ""), " -b ", paste(prefix, "_coding_window_temp.bed", sep = ""), " -wa | sort | uniq | awk {'print $1\"\t\"$2+1\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$7'} > ",paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""), sep = "")
    system(command)
    command = paste("perl ", annovar_folder, "annotate_variation.pl -geneanno -buildver hg19 ", paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""), " ", annovar_folder, "humandb/", sep = "")
    system(command)
    coding_anno = read.delim(paste(paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""),".exonic_variant_function", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    removed_id = coding_anno[coding_anno[,2] != "synonymous SNV",10]
    # now remove non-synonymous mutations (every coding mutations that are not synonymous mutations)
    mut = mut[!is.element(mut[,4], removed_id),]
  }
  mut[,2] = mut[,2] - 1  # change 1-based start to 0-based start. 
  write.table(mut, paste(prefix,"_temp.bed", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  command = paste("bedtools coverage -a ", paste(prefix,"_temp.bed", sep = ""), " -b ", window_file, " > ", paste(prefix,"_temp_coverage.bed", sep = ""), sep = "")
  system(command)
  coverage = read.delim(paste(prefix,"_temp_coverage.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  coverage = coverage[,1:5]
  colnames(coverage) = c("chr","start","end","site_index","mut_count")
  coverage$window_size = coverage$end - coverage$start # some regions don't have full length of, say, 50bp
  window = read.delim(window_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  cg = read.delim(cg_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(cg) = c("site_index", "cg")
  coverage = merge(coverage, cg, by.x = "site_index", by.y ="site_index")
  coverage$coding = as.numeric(mapply(function(x,y) grepl(y,x),coverage$site_index, MoreArgs = list(y="coding")))
  coverage$promoter = as.numeric(mapply(function(x,y) grepl(y,x),coverage$site_index, MoreArgs = list(y="promoter")))
  coverage$nf = as.numeric(mapply(function(x,y) grepl(y,x),coverage$site_index, MoreArgs = list(y="nf")))
  mutrate = read.delim(mutrate_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  mutrate = data.frame(site_index = mutrate[,1], mutrate = mutrate[,4])
  coverage = merge(coverage, mutrate, by.x = "site_index", by.y = "site_index")
  coverage = coverage[coverage$mutrate !=0,]
  #fit mutation rate model using a fixed set of features not including histone modification marks. 
  out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
  # get the adjusted mutation rate per base per individual
  coverage = data.frame(coverage, adjusted_mutrate = out.offset$fitted.values/(2*sample_size))
  # assign gene name to windows
  gene_assign = read.delim(gene_assign_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_assign) = c("site_index", "genename")
  coverage = merge(coverage, gene_assign, by.x = "site_index", by.y = "site_index")
  # get the piror probability of genes.
  gene_prior = read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_prior) = c("genename", "prior")
  # change prior to the 
  gene_prior$prior = 1- gene_prior$prior
  
  if(report_proportion !=1){
    genes_for_report = gene_prior[order(gene_prior[,2],decreasing = TRUE),1]
    genes_for_report = genes_for_report[1:floor(length(genes_for_report)*report_proportion)]
    coverage = coverage[is.element(coverage$genename, genes_for_report),]
  }
  
  # get the adjustion effect for each window by dividing adjusted_mutrate by mutrate
  coverage$adjusting_effect = coverage$adjusted_mutrate/coverage$mutrate
  
  # now need to extropolate window mutation file to base level mutation file, now removes coding region, and don't differntiate between promoter and nf
  coverage_noncoding <-coverage[coverage$coding != 1,c("chr","start","end","adjusting_effect","genename")]
  coverage_noncoding$ID <- paste(coverage_noncoding$genename, coverage_noncoding$start, sep = "_")
  
  #funtion to expand windows to base level
  window_expansion <- function(table_row){
    start <- seq(as.numeric(table_row[2]),as.numeric(table_row[3])-1)
    data.frame(table_row[1], start, start+1, paste(table_row[5],start,sep = "_"), table_row[6])
  }
  
  options(warn=-1)
  #options(scipen=999)
  cl <- makeCluster(node_n)
  #clusterExport(cl, "coverage_noncoding",  environment()) # need to tell clusterExport to look for window_expansion in the current environment inside of the function
  # use parallel computing and rbindlist to increase the speed by thousands of times. 
  environment(window_expansion) <- .GlobalEnv
  coverage_noncoding_for_base_mutrate <- rbindlist(parApply(cl, coverage_noncoding, 1, window_expansion))
  stopCluster(cl)
  options(warn = 0)
  #options(scipen=0)
  colnames(coverage_noncoding_for_base_mutrate) <- c("chr","start","end","base_ID","ID")
  coverage_noncoding_for_base_mutrate$start <- as.integer(coverage_noncoding_for_base_mutrate$start)
  coverage_noncoding_for_base_mutrate$end <- as.integer(coverage_noncoding_for_base_mutrate$end)
  
  # write out a bed file to get base-level mutation rates
  fwrite(coverage_noncoding_for_base_mutrate[,1:4],paste(prefix, "_temp_for_mutrate.bed", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  command <- paste("../lib/bigWigAverageOverBed ", mutrate_ref_file, " ", paste(prefix, "_temp_for_mutrate.bed", sep = ""), " ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = "" ), sep = "")
  system(command)
  command <- paste("awk {'print $1\"\t\"$4'} ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = ""), " > ", paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), sep = "")
  system(command)
  coverage_noncoding_base_mutrate <-fread(paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  colnames(coverage_noncoding_base_mutrate) <- c("base_ID","base_mutrate")
  a_temp <- coverage_noncoding_base_mutrate[coverage_noncoding_for_base_mutrate, on = "base_ID"]
  a_temp <- a_temp[as.data.table(coverage_noncoding), on = "ID"]
  coverage_noncoding_mutrate_adjusted <- a_temp[,c("base_ID","chr","start","end","genename")]
  coverage_noncoding_mutrate_adjusted$adjusted_base_mutrate <- a_temp[,c("base_mutrate")]*a_temp[,c("adjusting_effect")]
  rm(a_temp)
  
  # get the mutation count for each base
  command = paste("bedtools coverage -a ", paste(prefix,"_temp.bed", sep = ""), " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""), " > ", paste(prefix,"_temp_base_level_coverage.bed", sep = ""), sep = "")
  system(command)
  
  # merge table
  mutation_overlap_base <-fread(paste(prefix,"_temp_base_level_coverage.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  mutation_overlap_base <- mutation_overlap_base[,c("V4","V5")]
  colnames(mutation_overlap_base) <- c("base_ID","mut_count")
  coverage_noncoding_mutrate_adjusted <- coverage_noncoding_mutrate_adjusted[mutation_overlap_base, on = "base_ID"]
  rm(mutation_overlap_base)
  
  # overlap with epi feature.
  command <- paste("bedtools intersect -f ", overlap, " -a ", paste(prefix, "_temp_for_mutrate.bed", sep = "") , " -b ", epigenomic_marks, " -wa > ", paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""),sep = "")
  system(command)
  base_in_epi <- fread(paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  base_in_epi <- data.table(base_ID = base_in_epi[,4], epi = 1)
  colnames(base_in_epi) <- c("base_ID", "epi")
  coverage_noncoding_mutrate_adjusted <- base_in_epi[coverage_noncoding_mutrate_adjusted, on = "base_ID"]
  coverage_noncoding_mutrate_adjusted[is.na(get("epi")), ("epi"):=0]
  system(paste("rm ", prefix, "_temp*", sep = ""))
  coverage_noncoding_mutrate_adjusted
}


#This version works in RCC, or any cluster with a reasonablly big memory. 
# One thing different is that, when generating base-level mutation rate and features, the 50-bp windows data.frame will be partitioned to 20 parts
# and each part will be processed to get base-level information and combined into a txt file as the output of this function
adjust_mutation_rate_window_to_base_RCC <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks = "no", overlap = 0.5, rm_nonsyn = FALSE, annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", annovar_input = "no.txt", 
                                                    gene_assign_file, gene_prior_file, report_proportion,  mutrate_ref_file, node_n =6, output_file){
  # [mut_file] is the mutation location file from one dataset. e.g.,../analysis/160229_data_for_analysis.bed, chr \t start \t end \t index,  1-based start and end
  # [window_file] is the file with genome partitioned into windows. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed
  # [cg_file] is the file that has CG-content calculated, matching [window_file]. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.cg
  # [mutrate_file] is the file that has mutrate data. e.g., "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate"
  # [sample_size] is the number of individuals
  # [epitenomic_marks] if it is "no", then glm won't use if a mutation is in h3k27ac as a categorical variable, if it is a epigenomic annotation name.
  # for example, "../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed". It is fine to directly use original epigenomics bed files
  # instead of the partitioned files. Because the window_file already has annotations for each window in terms of whetehr they are promoter or coding. 
  # [overlap] the fraction of overlap (the paramter f in bedtools intersect), default is 0.5, For the ease of TADA modeling, it is better to assign it to be 1E-9(1bp)
  # [rm_nonsyn] TRUE or FALSE. Indicating whether nonysnonymous mutations should be removed when fitting the model. Use TRUE when fitting model using ASD mutations. 
  # [annovar_folder] is where the annovar program is, e.g., . e.g., /media/yuwen/Elements/ANNOVAR/annovar/
  # [annovar_input] is a input file generated from mutations that will be annotated by annovar. The format is chr \t start \t end \t ref_allele \t alt_allele \t index 1-based for starts and ends
  # [gene_assign_file] is the file that assigns each window to a gene. for example ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename.bed
  # [gene_prior_file], a file that has prior (derived from posterior and prior)for a gene as a risk gene. genename \t probability e.g. ../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.94.txt
  # [report_proportion] Choose the top X% TADA genes to report. 
  # [node_n] is the number of nodes used to run , default is 6
  # [mutrate_ref_file] the file with base level mutation rate reference. e.g., "/media/yuwen/Elements/mutation_rate_raw/MarkDaly/Daly_mutrate.bw"
  # [output_file] the name of the file that stores the output from the function (i.e., adjusted mutrate for each base, genename, mutcount....)
  prefix=system("date +%s", intern = TRUE) # prefix for temporary files that will be deleted at the end of the pipeline
  
  mut = read.delim(mut_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  if(rm_nonsyn){
    # command to transform the start base to 0-based in order to use bedtools to do overlap
    command = paste("awk {'print $1\"\t\"$2-1\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$7'} ", annovar_input, " > ", paste(prefix, "_annovar_for_bedtools.bed",sep = ""),sep = "")
    system(command)
    command = paste("grep coding ", window_file, " > ", paste(prefix, "_coding_window_temp.bed", sep = ""), sep = "")
    system(command)
    command = paste("bedtools intersect -a ",paste(prefix, "_annovar_for_bedtools.bed",sep = ""), " -b ", paste(prefix, "_coding_window_temp.bed", sep = ""), " -wa | sort | uniq | awk {'print $1\"\t\"$2+1\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$7'} > ",paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""), sep = "")
    system(command)
    command = paste("perl ", annovar_folder, "annotate_variation.pl -geneanno -buildver hg19 ", paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""), " ", annovar_folder, "humandb/", sep = "")
    system(command)
    coding_anno = read.delim(paste(paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""),".exonic_variant_function", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    removed_id = coding_anno[coding_anno[,2] != "synonymous SNV",10]
    # now remove non-synonymous mutations (every coding mutations that are not synonymous mutations)
    mut = mut[!is.element(mut[,4], removed_id),]
  }
  mut[,2] = mut[,2] - 1  # change 1-based start to 0-based start. 
  write.table(mut, paste(prefix,"_temp.bed", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  command = paste("bedtools coverage -a ", paste(prefix,"_temp.bed", sep = ""), " -b ", window_file, " > ", paste(prefix,"_temp_coverage.bed", sep = ""), sep = "")
  system(command)
  coverage = read.delim(paste(prefix,"_temp_coverage.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  coverage = coverage[,1:5]
  colnames(coverage) = c("chr","start","end","site_index","mut_count")
  coverage$window_size = coverage$end - coverage$start # some regions don't have full length of, say, 50bp
  window = read.delim(window_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  cg = read.delim(cg_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(cg) = c("site_index", "cg")
  coverage = merge(coverage, cg, by.x = "site_index", by.y ="site_index")
  coverage$coding = as.numeric(mapply(function(x,y) grepl(y,x),coverage$site_index, MoreArgs = list(y="coding")))
  coverage$promoter = as.numeric(mapply(function(x,y) grepl(y,x),coverage$site_index, MoreArgs = list(y="promoter")))
  coverage$nf = as.numeric(mapply(function(x,y) grepl(y,x),coverage$site_index, MoreArgs = list(y="nf")))
  mutrate = read.delim(mutrate_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  mutrate = data.frame(site_index = mutrate[,1], mutrate = mutrate[,4])
  coverage = merge(coverage, mutrate, by.x = "site_index", by.y = "site_index")
  coverage = coverage[coverage$mutrate !=0,]
  #fit mutation rate model using a fixed set of features not including histone modification marks. 
  out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
  # get the adjusted mutation rate per base per individual
  coverage = data.frame(coverage, adjusted_mutrate = out.offset$fitted.values/(2*sample_size))
  # assign gene name to windows
  gene_assign = read.delim(gene_assign_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_assign) = c("site_index", "genename")
  coverage = merge(coverage, gene_assign, by.x = "site_index", by.y = "site_index")
  # get the piror probability of genes.
  gene_prior = read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_prior) = c("genename", "prior")
  # change prior to the 
  gene_prior$prior = 1- gene_prior$prior
  
  if(report_proportion !=1){
    genes_for_report = gene_prior[order(gene_prior[,2],decreasing = TRUE),1]
    genes_for_report = genes_for_report[1:floor(length(genes_for_report)*report_proportion)]
    coverage = coverage[is.element(coverage$genename, genes_for_report),]
  }
  
  # get the adjustion effect for each window by dividing adjusted_mutrate by mutrate
  coverage$adjusting_effect = coverage$adjusted_mutrate/coverage$mutrate
  
  # now need to extropolate window mutation file to base level mutation file, now removes coding region, and don't differntiate between promoter and nf
  coverage_noncoding <-coverage[coverage$coding != 1,c("chr","start","end","adjusting_effect","genename")]
  coverage_noncoding$ID <- paste(coverage_noncoding$genename, coverage_noncoding$start, sep = "_")
  
  #funtion to expand windows to base level
  window_expansion <- function(table_row){
    start <- seq(as.numeric(table_row[2]),as.numeric(table_row[3])-1)
    data.frame(table_row[1], start, start+1, paste(table_row[5],start,sep = "_"), table_row[6])
  }
  
  options(warn=-1)
  #options(scipen=999)
  cl <- makeCluster(node_n)
  #clusterExport(cl, "coverage_noncoding",  environment()) # need to tell clusterExport to look for window_expansion in the current environment inside of the function
  # use parallel computing and rbindlist to increase the speed by thousands of times. 
  environment(window_expansion) <- .GlobalEnv
  total_rows <- nrow(coverage_noncoding)
  interval <- floor(total_rows/20)
  data_bins <- c(rep(seq(1,20), each = interval),rep(20, total_rows -interval*20))
  coverage_noncoding <- split(coverage_noncoding, data_bins)
  
  system(paste("rm -f ", output_file, sep = ""))
  for(i in 1:20){             
    # split coverage_noncoding to 10 parts, and for each part, generate feature table (which will be used for )
    coverage_noncoding_for_base_mutrate <- rbindlist(parApply(cl, coverage_noncoding[[i]], 1, window_expansion))
    
    #options(scipen=0)
    colnames(coverage_noncoding_for_base_mutrate) <- c("chr","start","end","base_ID","ID")
    coverage_noncoding_for_base_mutrate$start <- as.integer(coverage_noncoding_for_base_mutrate$start)
    coverage_noncoding_for_base_mutrate$end <- as.integer(coverage_noncoding_for_base_mutrate$end)
    
    # write out a bed file to get base-level mutation rates
    fwrite(coverage_noncoding_for_base_mutrate[,1:4],paste(prefix, "_temp_for_mutrate.bed", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    command <- paste("../lib/bigWigAverageOverBed ", mutrate_ref_file, " ", paste(prefix, "_temp_for_mutrate.bed", sep = ""), " ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = "" ), sep = "")
    system(command)
    command <- paste("awk {'print $1\"\t\"$4'} ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = ""), " > ", paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), sep = "")
    system(command)
    coverage_noncoding_base_mutrate <-fread(paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE, sep = "\t")
    colnames(coverage_noncoding_base_mutrate) <- c("base_ID","base_mutrate")
    a_temp <- coverage_noncoding_base_mutrate[coverage_noncoding_for_base_mutrate, on = "base_ID"]
    a_temp <- a_temp[as.data.table(coverage_noncoding[[i]]), on = "ID"]
    coverage_noncoding_mutrate_adjusted <- a_temp[,c("base_ID","chr","start","end","genename")]
    coverage_noncoding_mutrate_adjusted$adjusted_base_mutrate <- a_temp[,c("base_mutrate")]*a_temp[,c("adjusting_effect")]
    rm(a_temp)
    
    # get the mutation count for each base
    command = paste("bedtools coverage -a ", paste(prefix,"_temp.bed", sep = ""), " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""), " > ", paste(prefix,"_temp_base_level_coverage.bed", sep = ""), sep = "")
    system(command)
    mutation_overlap_base <-fread(paste(prefix,"_temp_base_level_coverage.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    mutation_overlap_base <- mutation_overlap_base[,c("V4","V5")]
    colnames(mutation_overlap_base) <- c("base_ID","mut_count")
    coverage_noncoding_mutrate_adjusted <- coverage_noncoding_mutrate_adjusted[mutation_overlap_base, on = "base_ID"]
    rm(mutation_overlap_base)
    
    # overlap with epi feature.
    command <- paste("bedtools intersect -f ", overlap, " -a ", paste(prefix, "_temp_for_mutrate.bed", sep = "") , " -b ", epigenomic_marks, " -wa > ", paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""),sep = "")
    system(command)
    base_in_epi <- fread(paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    base_in_epi <- data.table(base_ID = base_in_epi[,4], epi = 1)
    colnames(base_in_epi) <- c("base_ID", "epi")
    coverage_noncoding_mutrate_adjusted <- base_in_epi[coverage_noncoding_mutrate_adjusted, on = "base_ID"]
    coverage_noncoding_mutrate_adjusted[is.na(get("epi")), ("epi"):=0]
    if(i ==1){
      fwrite(coverage_noncoding_mutrate_adjusted, output_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)
    }
    else{
      fwrite(coverage_noncoding_mutrate_adjusted, output_file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)
    }
  }
  stopCluster(cl)
  options(warn = 0)
  system(paste("rm ", prefix, "_temp*", sep = ""))
  
}


#This version works will output base level information in a compact way. e.g., categorization will be done at gene level and feature configuration level.
# When generating base-level mutation rate and features, the 50-bp windows data.frame will be partitioned based on genename
# and each partition will be processed to get base-level information and combined into a txt file as the output of this function
adjust_mutation_rate_window_to_base_compact <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks = "no", overlap = 0.5, rm_nonsyn = FALSE, annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", annovar_input = "no.txt", 
                                                        gene_assign_file, gene_prior_file, report_proportion,  mutrate_ref_file, node_n =6, feature_start = 5, feature_end = 5, feature_number = 1, chunk_partition_num =20){
  # [mut_file] is the mutation location file from one dataset. e.g.,../analysis/160229_data_for_analysis.bed, chr \t start \t end \t index,  1-based start and end
  # [window_file] is the file with genome partitioned into windows. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed
  # [cg_file] is the file that has CG-content calculated, matching [window_file]. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.cg
  # [mutrate_file] is the file that has mutrate data. e.g., "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate"
  # [sample_size] is the number of individuals
  # [epitenomic_marks] if it is "no", then glm won't use if a mutation is in h3k27ac as a categorical variable, if it is a epigenomic annotation name.
  # for example, "../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed". It is fine to directly use original epigenomics bed files
  # instead of the partitioned files. Because the window_file already has annotations for each window in terms of whetehr they are promoter or coding. 
  # [overlap] the fraction of overlap (the paramter f in bedtools intersect), default is 0.5, For the ease of TADA modeling, it is better to assign it to be 1E-9(1bp)
  # [rm_nonsyn] TRUE or FALSE. Indicating whether nonysnonymous mutations should be removed when fitting the model. Use TRUE when fitting model using ASD mutations. 
  # [annovar_folder] is where the annovar program is, e.g., . e.g., /media/yuwen/Elements/ANNOVAR/annovar/
  # [annovar_input] is a input file generated from mutations that will be annotated by annovar. The format is chr \t start \t end \t ref_allele \t alt_allele \t index 1-based for starts and ends
  # [gene_assign_file] is the file that assigns each window to a gene. for example ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename.bed
  # [gene_prior_file], a file that has prior (derived from posterior and prior)for a gene as a risk gene. genename \t probability e.g. ../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.94.txt
  # [report_proportion] Choose the top X% TADA genes to report. 
  # [node_n] is the number of nodes used to run , default is 6
  # [mutrate_ref_file] the file with base level mutation rate reference. e.g., "/media/yuwen/Elements/mutation_rate_raw/MarkDaly/Daly_mutrate.bw"
  # [feature_start = 5, feature_end = 5, feature_number = 1] are the arguments to define the start position and end position of features that might have relative risk. 5 is for the current setting where only 1 binary feature is considered.
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
  coverage$coding <- as.numeric(mapply(function(x,y) grepl(y,x),coverage$site_index, MoreArgs = list(y="coding")))
  coverage$promoter <- as.numeric(mapply(function(x,y) grepl(y,x),coverage$site_index, MoreArgs = list(y="promoter")))
  coverage$nf <- as.numeric(mapply(function(x,y) grepl(y,x),coverage$site_index, MoreArgs = list(y="nf")))
  mutrate <- fread(mutrate_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  mutrate<-mutrate[,c("V1","V4"),with = FALSE]
  colnames(mutrate) <- c("site_index", "mutrate")
  coverage <-coverage[mutrate, on = "site_index"]
  coverage <- coverage[mutrate !=0]
  #fit mutation rate model using a fixed set of features not including histone modification marks. 
  out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
  # get the adjusted mutation rate per base per individual
  coverage$adjusted_mutrate = out.offset$fitted.values/(2*sample_size)
  # assign gene name to windows
  gene_assign <- fread(gene_assign_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_assign) <- c("site_index", "genename")
  coverage <- coverage[gene_assign, on = "site_index"]
  coverage<-coverage[complete.cases(coverage)]
  # get the piror probability of genes.
  gene_prior <- read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_prior) <- c("genename", "prior")
  # change prior to the posterior probability of alternative (risk) hypothesis
  gene_prior$prior <- 1- gene_prior$prior
  # merge gene prior info
  coverage <- coverage[gene_prior, on = "genename"]
  coverage<-coverage[complete.cases(coverage)]
  
  # remove genes based on TADA prior probability
  if(report_proportion !=1){
    genes_for_report <- gene_prior[order(gene_prior[,2],decreasing = TRUE),1]
    genes_for_report <- genes_for_report[1:floor(length(genes_for_report)*report_proportion)]
    coverage <- coverage[is.element(coverage$genename, genes_for_report),]
  }
  
  # get the adjustion effect for each window by dividing adjusted_mutrate by mutrate
  coverage$adjusting_effect = coverage$adjusted_mutrate/coverage$mutrate
  
  # now need to extropolate window mutation file to base level mutation file, now removes coding region, and don't differntiate between promoter and nf
  coverage_noncoding <-coverage[coding != 1,c("chr","start","end","adjusting_effect","genename"),with = FALSE]
  coverage_noncoding$ID <- paste(coverage_noncoding$genename, coverage_noncoding$start, sep = "_")
  
  total_rows <- nrow(coverage_noncoding)
  interval <- floor(total_rows/chunk_partition_num)
  data_bins <- c(rep(seq(1,chunk_partition_num), each = interval),rep(chunk_partition_num, total_rows -interval*chunk_partition_num))
  
  # split into 20 different chunks, then for each chunk split by genes, then by feature configuration. If split by gene at the first level, implementation would take too long
  
  coverage_noncoding <- split(coverage_noncoding, data_bins)
  
  # function to get effective information of each element of partition_by_gene
  # These information are those necessary to compute log-likelihood in the optimization function
  partition_feature <- function(pbg){
    # input is one element of the list of partition_by_gene
    pbg_split <- split(pbg, pbg[,feature_start:feature_end],drop = TRUE)
    feature_combination_number <- length(pbg_split)
    # this function below is different from the function used in dealing with dataset without reading by chunk. Here, prior is not incoporated at this step.
    info_for_each_feature <- function(feature_set){
      list(feature_vector = c(as.numeric(feature_set[1,feature_start:feature_end])), sum_mut_rate_count = sum(feature_set$mut_count*log(feature_set$adjusted_base_mutrate)), sum_mut_rate = sum(feature_set$adjusted_base_mutrate), sum_mut_count = sum(feature_set$mut_count), log_fcount = sum(log(factorial(feature_set$mut_count))))
    }
    sapply(pbg_split, info_for_each_feature,simplify = FALSE)
  }
  
  
  #funtion to expand windows to base level
  window_expansion <- function(table_row){
    start <- seq(as.numeric(table_row[2]),as.numeric(table_row[3])-1)
    data.frame(table_row[1], start, start+1, paste(table_row[5],start,sep = "_"), table_row[6])
  }
  
  options(warn=-1)
  cl <- makeCluster(node_n)
  # use parallel computing and rbindlist to increase the speed by thousands of times. 
  environment(window_expansion) <- .GlobalEnv
  
  # build vectors and list to store results
  partition_time1 <-c()
  partition_time2 <-c()
  data_partition <-list()
  
  
  for(i in 1:chunk_partition_num){
    # split coverage_noncoding to 10 parts, and for each part, generate feature table (which will be used for )
    coverage_noncoding_for_base_mutrate <- rbindlist(parApply(cl, coverage_noncoding[[i]], 1, window_expansion))
    
    colnames(coverage_noncoding_for_base_mutrate) <- c("chr","start","end","base_ID","ID")
    coverage_noncoding_for_base_mutrate$start <- as.integer(coverage_noncoding_for_base_mutrate$start)
    coverage_noncoding_for_base_mutrate$end <- as.integer(coverage_noncoding_for_base_mutrate$end)
    
    # write out a bed file to get base-level mutation rates
    fwrite(coverage_noncoding_for_base_mutrate[,1:4],paste(prefix, "_temp_for_mutrate.bed", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    command <- paste("../lib/bigWigAverageOverBed ", mutrate_ref_file, " ", paste(prefix, "_temp_for_mutrate.bed", sep = ""), " ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = "" ), sep = "")
    system(command)
    command <- paste("awk {'print $1\"\t\"$4'} ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = ""), " > ", paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), sep = "")
    system(command)
    coverage_noncoding_base_mutrate <-fread(paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE, sep = "\t")
    colnames(coverage_noncoding_base_mutrate) <- c("base_ID","base_mutrate")
    a_temp <- coverage_noncoding_for_base_mutrate[coverage_noncoding_base_mutrate, on = "base_ID"]
    a_temp <- a_temp[coverage_noncoding[[i]], on = "ID"]
    coverage_noncoding_mutrate_adjusted <- a_temp[,c("base_ID","chr","start","end","genename")]
    coverage_noncoding_mutrate_adjusted$adjusted_base_mutrate <- a_temp[,c("base_mutrate")]*a_temp[,c("adjusting_effect")]*2*sample_size # remember need to add back sample_size
    rm(a_temp)
    
    # get the mutation count for each base
    command = paste("bedtools coverage -a ", paste(prefix,"_temp.bed", sep = ""), " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""), " > ", paste(prefix,"_temp_base_level_coverage.bed", sep = ""), sep = "")
    system(command)
    mutation_overlap_base <-fread(paste(prefix,"_temp_base_level_coverage.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    mutation_overlap_base <- mutation_overlap_base[,c("V4","V5")]
    colnames(mutation_overlap_base) <- c("base_ID","mut_count")
    coverage_noncoding_mutrate_adjusted <- coverage_noncoding_mutrate_adjusted[mutation_overlap_base, on = "base_ID"]
    rm(mutation_overlap_base)
    
    # overlap with epi feature.
    command <- paste("bedtools coverage -a ", overlap, " -a ", epigenomic_marks, " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""),sep = "")
    system(command)
    base_in_epi <- fread(paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    base_in_epi <- base_in_epi[,c("V4","V5"), with = FALSE]
    colnames(base_in_epi) <- c("base_ID", "epi")
    coverage_noncoding_mutrate_adjusted <- coverage_noncoding_mutrate_adjusted[base_in_epi, on = "base_ID"]
    
    coverage_noncoding_mutrate_adjusted<- coverage_noncoding_mutrate_adjusted[,c("genename","base_ID","mut_count","adjusted_base_mutrate","epi")]
    
    # remove bases without mutation rates
    coverage_noncoding_mutrate_adjusted <- coverage_noncoding_mutrate_adjusted[adjusted_base_mutrate!=0]
    
    # remove rows with NA
    coverage_noncoding_mutrate_adjusted <- coverage_noncoding_mutrate_adjusted[complete.cases(coverage_noncoding_mutrate_adjusted)]
    
    #then partition by gene for the current chunk
    partition_time1[i] <- system.time(coverage_noncoding_mutrate_adjusted<- split(coverage_noncoding_mutrate_adjusted, coverage_noncoding_mutrate_adjusted$genename))
    # then partition by feature configuration for each gene in the current chunk
    partition_time2[i] <- system.time(coverage_noncoding_mutrate_adjusted <- sapply(coverage_noncoding_mutrate_adjusted, partition_feature, simplify = FALSE))
    data_partition <- append(data_partition, coverage_noncoding_mutrate_adjusted)
  }
  
  stopCluster(cl)
  options(warn = 0)
  system(paste("rm ", prefix, "_temp*", sep = ""))
  
  list(partition_time1 = partition_time1, partition_time2 = partition_time2, base_info = data_partition)
  
}


#This version works will output base level information in a compact way. e.g., categorization will be done at gene level and feature configuration level.
# When generating base-level mutation rate and features, the 50-bp windows data.frame will be partitioned based on genename
# and each partition will be processed to get base-level information and combined into a txt file as the output of this function
# more features that might have effect size is added and binaried. 
adjust_mutation_rate_window_to_base_compact_more_feature <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks = "no", sequence_annotation = c("phylop_100way"), 
                                                                     sequence_annotation_cutoff = list(phylop_100way=2), sequence_annotation_ref = list(phylop_100way = "../other_annotation/conservation_annotation/hg19.100way.phyloP100way.bw"),
                                                                     overlap = 0.5, rm_nonsyn = FALSE, annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", annovar_input = "no.txt", 
                                                                     gene_assign_file, gene_prior_file, report_proportion,  mutrate_ref_file = "../other_annotation/mutation_rate/Daly_mutrate.bw", 
                                                                     node_n =6, feature_start = 5, feature_end = 6, feature_number = 2, chunk_partition_num =20){
  # [mut_file] is the mutation location file from one dataset. e.g.,../analysis/160229_data_for_analysis.bed, chr \t start \t end \t index,  1-based start and end
  # [window_file] is the file with genome partitioned into windows. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed
  # [cg_file] is the file that has CG-content calculated, matching [window_file]. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.cg
  # [mutrate_file] is the file that has mutrate data. e.g., "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate"
  # [sample_size] is the number of individuals
  # [epitenomic_marks] if it is "no", then glm won't use if a mutation is in h3k27ac as a categorical variable, if it is a epigenomic annotation name.
  # for example, "../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed". It is fine to directly use original epigenomics bed files
  # instead of the partitioned files. Because the window_file already has annotations for each window in terms of whetehr they are promoter or coding. 
  # [sequence_annotation] the type of base-level annotation that will be added to estimate its relative risk
  # [sequence_annotation_cutoff] the cutoff for one type of sequence_annotation
  # [sequence_annoation_ref] the base level annotation, should be a .bw file e.g. "../other../other_annotation/conservation_annotation/hg19.100way.phyloP100way.bw"
  
  # [overlap] the fraction of overlap (the paramter f in bedtools intersect), default is 0.5, For the ease of TADA modeling, it is better to assign it to be 1E-9(1bp)
  # [rm_nonsyn] TRUE or FALSE. Indicating whether nonysnonymous mutations should be removed when fitting the model. Use TRUE when fitting model using ASD mutations. 
  # [annovar_folder] is where the annovar program is, e.g., . e.g., /media/yuwen/Elements/ANNOVAR/annovar/
  # [annovar_input] is a input file generated from mutations that will be annotated by annovar. The format is chr \t start \t end \t ref_allele \t alt_allele \t index 1-based for starts and ends
  # [gene_assign_file] is the file that assigns each window to a gene. for example ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename.bed
  # [gene_prior_file], a file that has prior (derived from posterior and prior)for a gene as a risk gene. genename \t probability e.g. ../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.94.txt
  # [report_proportion] Choose the top X% TADA genes to report. 
  # [node_n] is the number of nodes used to run , default is 6
  # [mutrate_ref_file] the file with base level mutation rate reference. e.g., "/media/yuwen/Elements/mutation_rate_raw/MarkDaly/Daly_mutrate.bw"
  # [feature_start = 5, feature_end = 5, feature_number = 1] are the arguments to define the start position and end position of features that might have relative risk. 5 is for the current setting where only 1 binary feature is considered.
  # [chunk_partition_num] is the chunk number that will be used to partition the 50-bp window info file, so that each partition will be extended to base-level sequencially to decrease RAM burden
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
  # get the adjusted mutation rate per base per individual
  coverage$adjusted_mutrate = out.offset$fitted.values/(2*sample_size)
  # assign gene name to windows
  gene_assign <- fread(gene_assign_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_assign) <- c("site_index", "genename")
  coverage <- coverage[gene_assign, on = "site_index"]
  coverage<-coverage[complete.cases(coverage)]
  # get the piror probability of genes.
  gene_prior <- read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_prior) <- c("genename", "prior")
  # change prior to the posterior probability of alternative (risk) hypothesis
  gene_prior$prior <- 1- gene_prior$prior
  # merge gene prior info
  coverage <- coverage[gene_prior, on = "genename"]
  coverage<-coverage[complete.cases(coverage)]
  
  # remove genes based on TADA prior probability
  if(report_proportion !=1){
    genes_for_report <- gene_prior[order(gene_prior[,2],decreasing = TRUE),1]
    genes_for_report <- genes_for_report[1:floor(length(genes_for_report)*report_proportion)]
    coverage <- coverage[is.element(coverage$genename, genes_for_report),]
  }
  
  # get the adjustion effect for each window by dividing adjusted_mutrate by mutrate
  coverage$adjusting_effect = coverage$adjusted_mutrate/coverage$mutrate
  
  # now need to extropolate window mutation file to base level mutation file, now removes coding region, and don't differntiate between promoter and nf
  coverage_noncoding <-coverage[coding != 1,c("chr","start","end","adjusting_effect","genename"),with = FALSE]
  coverage_noncoding$ID <- paste(coverage_noncoding$genename, coverage_noncoding$start, sep = "_")
  
  total_rows <- nrow(coverage_noncoding)
  interval <- floor(total_rows/chunk_partition_num)
  data_bins <- c(rep(seq(1,chunk_partition_num), each = interval),rep(chunk_partition_num, total_rows -interval*chunk_partition_num))
  
  # split into 20 different chunks, then for each chunk split by genes, then by feature configuration. If split by gene at the first level, implementation would take too long
  
  coverage_noncoding <- split(coverage_noncoding, data_bins)
  
  # function to get effective information of each element of partition_by_gene
  # These information are those necessary to compute log-likelihood in the optimization function
  partition_feature <- function(pbg){
    # input is one element of the list of partition_by_gene
    pbg_split <- split(pbg, pbg[,feature_start:feature_end],drop = TRUE)
    feature_combination_number <- length(pbg_split)
    # this function below is different from the function used in dealing with dataset without reading by chunk. Here, prior is not incoporated at this step.
    info_for_each_feature <- function(feature_set){
      list(feature_vector = c(as.numeric(feature_set[1,feature_start:feature_end])), sum_mut_rate_count = sum(feature_set$mut_count*log(feature_set$adjusted_base_mutrate)), sum_mut_rate = sum(feature_set$adjusted_base_mutrate), sum_mut_count = sum(feature_set$mut_count), log_fcount = sum(log(factorial(feature_set$mut_count))))
    }
    sapply(pbg_split, info_for_each_feature,simplify = FALSE)
  }
  
  
  #funtion to expand windows to base level
  window_expansion <- function(table_row){
    start <- seq(as.numeric(table_row[2]),as.numeric(table_row[3])-1)
    data.frame(table_row[1], start, start+1, paste(table_row[5],start,sep = "_"), table_row[6])
  }
  
  options(warn=-1)
  cl <- makeCluster(node_n)
  # use parallel computing and rbindlist to increase the speed by thousands of times. 
  environment(window_expansion) <- .GlobalEnv
  
  # build vectors and list to store results
  partition_time1 <-c()
  partition_time2 <-c()
  data_partition <-list()
  
  
  for(i in 1:chunk_partition_num){
    # split coverage_noncoding to 10 parts, and for each part, generate feature table (which will be used for )
    coverage_noncoding_for_base_mutrate <- rbindlist(parApply(cl, coverage_noncoding[[i]], 1, window_expansion))
    
    colnames(coverage_noncoding_for_base_mutrate) <- c("chr","start","end","base_ID","ID")
    coverage_noncoding_for_base_mutrate$start <- as.integer(coverage_noncoding_for_base_mutrate$start)
    coverage_noncoding_for_base_mutrate$end <- as.integer(coverage_noncoding_for_base_mutrate$end)
    
    # write out a bed file to get base-level mutation rates
    fwrite(coverage_noncoding_for_base_mutrate[,1:4],paste(prefix, "_temp_for_mutrate.bed", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    command <- paste("../lib/bigWigAverageOverBed ", mutrate_ref_file, " ", paste(prefix, "_temp_for_mutrate.bed", sep = ""), " ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = "" ), sep = "")
    system(command)
    command <- paste("awk {'print $1\"\t\"$4'} ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = ""), " > ", paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), sep = "")
    system(command)
    coverage_noncoding_base_mutrate <-fread(paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE, sep = "\t")
    colnames(coverage_noncoding_base_mutrate) <- c("base_ID","base_mutrate")
    a_temp <- coverage_noncoding_for_base_mutrate[coverage_noncoding_base_mutrate, on = "base_ID"]
    a_temp <- a_temp[coverage_noncoding[[i]], on = "ID"]
    coverage_noncoding_mutrate_adjusted <- a_temp[,c("base_ID","chr","start","end","genename")]
    coverage_noncoding_mutrate_adjusted$adjusted_base_mutrate <- a_temp[,c("base_mutrate")]*a_temp[,c("adjusting_effect")]*2*sample_size # remember need to add back sample_size
    rm(a_temp)
    
    # get the mutation count for each base
    command = paste("bedtools coverage -a ", paste(prefix,"_temp.bed", sep = ""), " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""), " > ", paste(prefix,"_temp_base_level_coverage.bed", sep = ""), sep = "")
    system(command)
    mutation_overlap_base <-fread(paste(prefix,"_temp_base_level_coverage.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    mutation_overlap_base <- mutation_overlap_base[,c("V4","V5")]
    colnames(mutation_overlap_base) <- c("base_ID","mut_count")
    coverage_noncoding_mutrate_adjusted <- coverage_noncoding_mutrate_adjusted[mutation_overlap_base, on = "base_ID"]
    rm(mutation_overlap_base)
    
    # overlap with epi feature.
    command <- paste("bedtools coverage -a ", overlap, " -a ", epigenomic_marks, " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""),sep = "")
    system(command)
    base_in_epi <- fread(paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    base_in_epi <- base_in_epi[,c("V4","V5"), with = FALSE]
    colnames(base_in_epi) <- c("base_ID", "epi")
    coverage_noncoding_mutrate_adjusted <- coverage_noncoding_mutrate_adjusted[base_in_epi, on = "base_ID"]
    coverage_noncoding_mutrate_adjusted<- coverage_noncoding_mutrate_adjusted[,c("genename","base_ID","mut_count","adjusted_base_mutrate","epi")]
    
    # get the sequence annotation feature.
    if(sequence_annotation != "no"){
      if(is.element("phylop_100way", sequence_annotation)){ # if phylop annotation needs to be added
        cutoff <- sequence_annotation_cutoff[["phylop_100way"]]
        # run commands to get base level phylop annotation
        ref_file <- sequence_annotation_ref[["phylop_100way"]]
        command <- paste("../lib/bigWigAverageOverBed ", ref_file, " ", paste(prefix, "_temp_for_mutrate.bed", sep = ""), " ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = "" ), sep = "")
        system(command)
        command <- paste("awk {'print $1\"\t\"$4'} ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = ""), " > ", paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), sep = "")
        system(command)
        coverage_noncoding_base_anno <-fread(paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE, sep = "\t")
        colnames(coverage_noncoding_base_anno) <- c("base_ID","phylop_100way")
        coverage_noncoding_mutrate_adjusted <- coverage_noncoding_mutrate_adjusted[coverage_noncoding_base_anno, on = "base_ID"]
        # binarize annotation score
        coverage_noncoding_mutrate_adjusted[,phylop_100way:= as.numeric(.SD > cutoff ), .SDcol = c("phylop_100way")]
      }
    }
    
    # remove bases without mutation rates
    coverage_noncoding_mutrate_adjusted <- coverage_noncoding_mutrate_adjusted[adjusted_base_mutrate!=0]
    
    # remove rows with NA
    coverage_noncoding_mutrate_adjusted <- coverage_noncoding_mutrate_adjusted[complete.cases(coverage_noncoding_mutrate_adjusted)]
    
    #then partition by gene for the current chunk
    partition_time1[i] <- system.time(coverage_noncoding_mutrate_adjusted<- split(coverage_noncoding_mutrate_adjusted, coverage_noncoding_mutrate_adjusted$genename))
    # then partition by feature configuration for each gene in the current chunk
    partition_time2[i] <- system.time(coverage_noncoding_mutrate_adjusted <- sapply(coverage_noncoding_mutrate_adjusted, partition_feature, simplify = FALSE))
    data_partition <- append(data_partition, coverage_noncoding_mutrate_adjusted)
  }
  
  stopCluster(cl)
  options(warn = 0)
  system(paste("rm ", prefix, "_temp*", sep = ""))
  
  list(partition_time1 = partition_time1, partition_time2 = partition_time2, base_info = data_partition)
  
}