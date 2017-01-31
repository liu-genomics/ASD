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