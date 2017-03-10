
#This version works will output base level information in a compact way. e.g., categorization will be done at gene level and feature configuration level.
# When generating base-level mutation rate and features, the 50-bp windows data.frame will be partitioned based on genename
# and each partition will be processed to get base-level information and combined into a txt file as the output of this function
# more features that might have effect size is added and binaried. 
# include even more features for epigenomic annotations
adjust_mutation_rate_window_to_base_compact_more_feature_v2 <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks_list = c("../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed","../other_annotation/epigenomic_annotation/Epigenome_E081_E082_intersection.bed","../other_annotation/epigenomic_annotation/encode_DHS_union.bed"),
                                                                        sequence_annotation = c("phylop_100way", "gerp"), 
                                                                     sequence_annotation_cutoff = list(phylop_100way=2, gerp =2 ), sequence_annotation_ref = list(phylop_100way = "../other_annotation/conservation_annotation/hg19.100way.phyloP100way.bw", gerp = "../other_annotations/conservation_annotation/hg19_gerp_score.bw"),
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
    
    coverage_noncoding_mutrate_adjusted<- coverage_noncoding_mutrate_adjusted[,c("genename","base_ID","mut_count","adjusted_base_mutrate")]
    
    # overlap with multiple epi feature.
    if (epigenomic_marks != "no"){ # then epigenomic_marks must be a vector of epigenomic bed files that need to be compard with the mutation data
      epi_ID = 1
      for(epi in epigenomic_marks){
        command <- paste("bedtools coverage -a ", overlap, " -a ", epi, " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""),sep = "")
        system(command)
        base_in_epi <- fread(paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        base_in_epi <- base_in_epi[,c("V4","V5"), with = FALSE]
        colnames(base_in_epi) <- c("base_ID", paste("epi",epi_ID, sep = "_"))
        coverage_noncoding_mutrate_adjusted <- coverage_noncoding_mutrate_adjusted[base_in_epi, on = "base_ID"]
        epi_ID <- epi_ID + 1
      }
    }
    
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
      if(is.element("gerp", sequence_annotation)){ # if phylop annotation needs to be added
        cutoff <- sequence_annotation_cutoff[["gerp"]]
        # run commands to get base level phylop annotation
        ref_file <- sequence_annotation_ref[["gerp"]]
        command <- paste("../lib/bigWigAverageOverBed ", ref_file, " ", paste(prefix, "_temp_for_mutrate.bed", sep = ""), " ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = "" ), sep = "")
        system(command)
        command <- paste("awk {'print $1\"\t\"$4'} ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = ""), " > ", paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), sep = "")
        system(command)
        coverage_noncoding_base_anno <-fread(paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE, sep = "\t")
        colnames(coverage_noncoding_base_anno) <- c("base_ID","gerp")
        coverage_noncoding_mutrate_adjusted <- coverage_noncoding_mutrate_adjusted[coverage_noncoding_base_anno, on = "base_ID"]
        # binarize annotation score
        coverage_noncoding_mutrate_adjusted[,gerp:= as.numeric(.SD > cutoff ), .SDcol = c("gerp")]
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



new_control = new.env()
load("../data/debug_region_list_073116_using_258_control_data_matrix.Rdata", new_control)

prefix <-system("date +%s", intern = TRUE) # prefix for temporary files that will be deleted at the end of the pipeline

write.table(new_control$data_matrix[new_control$data_matrix$Prediction == "dnv_proband" & new_control$data_matrix$Type == "SNV",c(1,2,3,7)], paste(prefix,"temp.bed",sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(new_control$data_matrix[new_control$data_matrix$Prediction == "dnv_proband" ,c(1:5,8,7)], paste(prefix,"annovar_input_temp.bed",sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)


mut_file = paste(prefix,"temp.bed",sep = "")
window_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed"
cg_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.fasta.cg"
mutrate_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.mutrate"
sample_size = 314
epigenomic_marks = c("../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed","../other_annotation/epigenomic_annotation/Epigenome_E081_E082_intersection.bed","../other_annotation/epigenomic_annotation/encode_DHS_union.bed")
sequence_annotation = c("phylop_100way", "gerp")
sequence_annotation_cutoff = list(phylop_100way=2, gerp =2 )
sequence_annotation_ref = list(phylop_100way = "../other_annotation/conservation_annotation/hg19.100way.phyloP100way.bw", gerp = "../other_annotations/conservation_annotation/hg19_gerp_score.bw")
overlap = 0.5
rm_nonsyn = FALSE
annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/"
annovar_input = "no.txt"
gene_assign_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename.bed"
gene_prior_file = "../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.94.txt"
report_proportion = 0.06
mutrate_ref_file = "/media/yuwen/Elements/mutation_rate_raw/MarkDaly/Daly_mutrate.bw"
node_n = 6
feature_start = 5
feature_end = 9
feature_number = 5
chunk_partition_num =1