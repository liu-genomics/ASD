# library(data.table)
# 
# set.seed(1234)
# n = 1e8
# 
# data_frame_1 = data.frame(id=paste("id_", 1:n, sep=""),
#                           factor1=sample(c("A", "B", "C"), n, replace=TRUE))
# data_frame_2 = data.frame(id=sample(data_frame_1$id),
#                           value1=rnorm(n))
# 
# data_table_1 = data.table(data_frame_1, key="id")
# data_table_2 = data.table(data_frame_2, key="id")
# 
# system.time(df.merged <- merge(data_frame_1, data_frame_2))
# #   user  system elapsed 
# # 17.983   0.189  18.063 
# 
# 
# system.time(dt.merged <- merge(data_table_1, data_table_2))
# #   user  system elapsed 
# #  0.729   0.099   0.821 

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



sample_size <- 3e8
data <- sim_table(sample_size)

running_time <- system.time(data.fit <- glm(count~ ep1 + ep2 + ep3 + ep4 + ep5 + annotation1 + annotation2, data = data, family = c("poisson")))

rm(data.fit)
save.image("../analysis/170102_simulate_3e8_data_points_for_5_binary_2_continuous_for_GLM.Rdata")
