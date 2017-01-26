
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
    all_rr <- c(log(data$base_mutrate[1]),x)
    
    sum_cate <-function(data_subset){
      exp(as.matrix(data_subset)[1,]%*%all_rr)*nrow(data_subset)
    }
    # summation of part1 of Equation 4
    
    sum_logP_1 <- sum(data_with_response$response*(as.matrix(data_with_response)[,-1]%*%all_rr))
    # summation of part2 of Equation 4
    sum_logP_2 <- sum(sapply(data.split, sum_cate))
    sum_logP_1 - sum_logP_2 
  }
  optimization_time <- system.time(mle <- optim(c(0.1), fr,control=list("fnscale"=-1, maxit = 5000), hessian = TRUE))  
  list(split_time = split_time,optimization_time = optimization_time, mle = mle)
}


# only process generated from sim_table_simple, use glm to fit the parameters
adjust_muatation_rate_optimization_test_simple_glm <- function(data){
  optimization_time <- system.time(out.offset <- glm(count ~ ep1+offset(log(data$base_mutrate)), family = poisson, data = data))

  list(optimization_time = optimization_time, mle = out.offset$coefficients)  
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