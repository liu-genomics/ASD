#function to take in a table like the following and other information that has been loaded into the workspace
enhancer_to_mutation <-function(enhancers,effective_mut_ID ,recur_number = 2){ # cutoff for the number of recurrences for each enhancers
  #effective_mut_ID could be 
  recur_enh = enhancers[enhancers$Freq >= recur_number,]
  
}

#generate the distribution of total number of enhancers with recurrent SNVs
recurrent_enhancers_number <- function(x){
  sum(x>=2)
}


# function to do simulations and get p-value
sim_recurrent <- function(mutation_bed, annotation_folder, suffix, simulation_times){
  #mutation_bed is the file name of mutations
  # annotation_folder is where the annotation files for recurrent enhancers are, default is ../other_annotation/epigenomic_annotation/bed_files_for_recurrent_analysis
  # suffix is the identifer to know which epigenomic mark is going to be used. e.g., Noonan_brain_roadmap_union.bed
  # simulation_times is the number of simulations to be done
  enhancer_file <- paste(annotation_folder, "recurrent_enhancer_", suffix, sep = "")
  enhancer_mutrate_file <- paste(annotation_folder, "recurrent_enhancer_mutrate_", suffix, sep = "")
  enhancer_mutrate <- read.delim(enhancer_mutrate_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  command = paste("bedtools intersect -a ",mutation_bed, " -b ", enhancer_file," -wa -wb | awk {'print $4\"\t\"$8'} | sort | uniq ",sep = "")
  temp = system(command, intern = TRUE)
  temp = as.data.frame(t(mapply(function(x) strsplit(x,split="\t")[[1]],temp)))
  #the 1st column of temp is mutation ID, the 2nd is enhancer ID.
  ASD_mutation_in_enhancer <- length(temp[is.element(temp[,1],mutation$ASD_effective_SNV_ID),1])
  control_mutation_in_enhancer <- length(temp[is.element(temp[,1],mutation$control_effective_SNV_ID),1])
  
  sample_mut_number <- min(ASD_mutation_in_enhancer,control_mutation_in_enhancer)
  temp1 <- temp[is.element(temp[,1],mutation$ASD_effective_SNV_ID),]
  ASD_mut_sample <- temp1[sample(seq(1:length(temp1[,1])),sample_mut_number),]
  temp2 <- temp[is.element(temp[,1],mutation$control_effective_SNV_ID),]
  control_mut_sample <- temp2[sample(seq(1:length(temp2[,1])),sample_mut_number),]
  
  
  enhancer_with_ASD <- as.data.frame(table(ASD_mut_sample[,2]))
  enhancer_with_ASD <- enhancer_with_ASD[enhancer_with_ASD[,2] > 0,]
  enhancer_with_control <- as.data.frame(table(control_mut_sample[,2]))
  enhancer_with_control <- enhancer_with_control[enhancer_with_control[,2] > 0,]
  ASD_recur_enhancer_number <- length(enhancer_with_ASD[enhancer_with_ASD[,2] >=2,2])
  control_recur_enhancer_number <- length(enhancer_with_control[enhancer_with_control[,2] >=2,2])
  
  #list(ASD_recur_enhancer_number, control_recur_enhancer_number, sample_mut_number, ASD_mutation_in_enhancer,control_mutation_in_enhancer)
  
 
   SNV_recur_count = c()
   for(i in 1:1) { # 1000 runs of simulations
     SNV_recur_count = c(SNV_recur_count,apply(data.frame(rmultinom(simulation_times,sample_mut_number,enhancer_mutrate[,4]/sum(enhancer_mutrate[,4]))),2,recurrent_enhancers_number))
   }
  ASD_p <- length(SNV_recur_count[SNV_recur_count >= ASD_recur_enhancer_number])/simulation_times
  control_p <- length(SNV_recur_count[SNV_recur_count >= control_recur_enhancer_number])/simulation_times
  list(ASD_pvalue = ASD_p, control_pvalue = control_p, ASD_observed = ASD_recur_enhancer_number,
       control_observed = control_recur_enhancer_number, effective_mutation_sampled = sample_mut_number,
       simulated_mean = mean(SNV_recur_count), simulated_lowerbound = quantile(SNV_recur_count, c(0.025)), 
       simulated_upperbound = quantile(SNV_recur_count, c(0.975)))
  
}