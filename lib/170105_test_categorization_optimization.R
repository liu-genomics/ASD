# The formal location of these two function files are in 170102_test_glm_functions.R

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

sample_size <- 3e8
data <- sim_table(sample_size)

write.table("data_finished","../analysis/170105_data_finished_temp.txt")

output <- adjust_muatation_rate_optimization_test(data)

rm(data)
save.image("../analysis/170105_simulate_3e8_data_points_for_5_binary_2_continuous_for_categorization.Rdata")
quit("no")


