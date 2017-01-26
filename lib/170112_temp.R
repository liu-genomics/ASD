# use verified parameter estimate methods (code) to estiamte effect sizes of active promoters and active enhancers. 
# The v4 version is different from that in the first version in that only a fixed set of parameters (cg, promote, coding, will be used to adjust for mutation rate. )
# The v4 version also included splicing mutations in estimating relative risks.
verified_effect_size_estimate_noncoding_mutations_v4 <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks_list = "no", overlap = 0.5, gene_assign_file,rm_nonsyn = FALSE, annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", annovar_input = "no.txt", gene_prior_file, 
                                                                 sequence_annotation_list = "no", optimization_prop = 1, ){
  # [mut_file] is the mutation location file from one dataset. e.g.,../analysis/160229_data_for_analysis.bed, chr \t start \t end \t index,  1-based start and end
  # [window_file] is the file with genome partitioned into windows. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed
  # [cg_file] is the file that has CG-content calculated, matching [window_file]. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.cg
  # [mutrate_file] is the file that has mutrate data. e.g., "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate"
  # [sample_size] is the number of individuals
  # [epitenomic_marks_list] if it is "no", then glm won't use if a mutation is in h3k27ac as a categorical variable, if it is a vector of epigenomic annotation names, epigenomic mark will be added one by one.
  # for example, c("../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed","../other_annotation/epigenomic_annotation/Epigenome_E081_E082_intersection.bed","../other_annotation/epigenomic_annotation/encode_DHS_union.bed"). It is fine to directly use original epigenomics bed files
  # instead of the partitioned files. Because the window_file already has annotations for each window in terms of whetehr they are promoter or coding. 
  # [overlap] the fraction of overlap (the paramter f in bedtools intersect), default is 0.5, For the ease of TADA modeling, it is better to assign it to be 1E-9(1bp)
  # [gene_assign_file] is the file that assigns each window to a gene. for example ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename.bed
  # [rm_nonsyn] TRUE or FALSE. Indicating whether nonysnonymous mutations should be removed when fitting the model. Use TRUE when fitting model using ASD mutations. 
  # [annovar_folder] is where the annovar program is, e.g., . e.g., /media/yuwen/Elements/ANNOVAR/annovar/
  # [annovar_input] is a input file generated from mutations that will be annotated by annovar. The format is chr \t start \t end \t ref_allele \t alt_allele \t index 1-based for starts and ends
  # [gene_prior_file], a file that has prior (derived from posterior and prior)for a gene as a risk gene. genename \t probability e.g. ../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.05.txt
  # [sequence_annotation_list] a list of annotation files for each window in the window_file, e.g., after extracting conservation score for each window. c("../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.phastcons100way")
  # [optimization_prop], the top percentage of genes that are used for optimization. genes are ranked high to low based on their Posterior probability of the alternative hypothesis 
  
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
  if(epigenomic_marks_list == "no"){
    out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
  }
  else{
    epi_index = 1 # start from the first epi marker, will add 1 to itself after incorporate each epigenomics mark. 
    for(epigenomic_marks in epigenomic_marks_list){
      command = paste("bedtools intersect -f ", overlap, " -a ", window_file , " -b ", epigenomic_marks, " -wa > ", paste(prefix,"_temp_overlap_epi.bed", sep = ""),sep = "")
      system(command)
      window_in_epi = read.delim(paste(prefix,"_temp_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      window_in_epi = data.frame(site_index = window_in_epi[,4])
      window_in_epi[,paste("epi", epi_index, sep="")] = 1
      coverage = merge(coverage, window_in_epi, by.x = "site_index", by.y = "site_index", all.x = TRUE)
      coverage[is.na(coverage[,paste("epi", epi_index, sep="")]),paste("epi", epi_index, sep="")] = 0
      epi_index = epi_index + 1
    }
    annotation_index = 1
    for(seq_annotation in sequence_annotation_list){
      annotation_raw = read.delim(seq_annotation, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      annotation = data.frame(site_index = annotation_raw[,1])
      annotation[,paste("annotation",annotation_index,sep = "")] = annotation_raw[,4]
      coverage = merge(coverage, annotation, by.x = "site_index", by.y = "site_index")
      annotation_index = annotation_index + 1
    }
    num_addition_par = length(epigenomic_marks_list) + length(sequence_annotation_list) # number of additional feature parameters 
    
    # only use a set of fixed predictors to run glm in order to adjust for mutation rates. 
    out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
    
    coverage = data.frame(coverage, adjusted_mutrate = out.offset$fitted.values)
    
    # get the number of synonymous mutations that are captured in the windows, which will be used to scale factor for splicing mutations
    syn_num = sum(coverage[coverage$coding == 1 & coverage$mut_count ==1,]$coding)
    
    # get mutations that have dpsi_zscore under some certain cutoff, and these mutations will be defined as splicing mutations
    splicing_threshold = quantile(mut_spidex$dpsi_zscore,seq(0,1,0.1),na.rm = TRUE)[2] # lower 10% 
    splicing_mut = mut_spedix[mut_spedix$dpsi_zscore <= splicing_threshold & !is.na(mut_spedix$dpsi_zscore),c("Chrom","Start","End","ID")]
    # transform data types: integer to numeric
    splicing_mut[,2] = as.numeric(as.character(splicing_mut[,2]))
    splicing_mut[,3] = as.numeric(as.character(splicing_mut[,3]))
    # change start site to 0-based
    splicing_mut[,2] = splicing_mut[,2] - 1
    write.table(splicing_mut, paste(prefix,"_temp_splicing_mut.bed",sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    
    command = paste("bedtools coverage -a ", paste(prefix,"_temp_splicing_mut.bed",sep = ""), " -b ", window_file, " > ", paste(prefix,"_temp_splicing_mut_overlap.bed", sep = ""), sep = "")
    system(command)
    splicing_mut_overlap_window = read.delim(paste(prefix,"_temp_splicing_mut_overlap.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    splicing_mut_overlap_window = splicing_mut_overlap_window[,4:5]
    colnames(splicing_mut_overlap_window) = c("site_index","raw_splicing_mut") # raw_splicing_mut is any mutations that are considered to be splicing mutations, regulatory mutations haven't been removed
    coverage = merge(coverage, splicing_mut_overlap_window, by = "site_index")
    
    #get splicing mutation number
    splicing_num = sum(coverage$raw_splicing_mut) # only get 51 for ASD_manuscript_cases
    
    
    gene_assign = read.delim(gene_assign_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(gene_assign) = c("site_index", "genename")
    
    coverage = merge(coverage, gene_assign, by.x = "site_index", by.y = "site_index")
    # now I need to use the prior probability of each gene to generate risk gene indicator and mutation numbers
    gene_prior = read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    #test_count = mapply(function(x){rpois(1,x)},mutrate_temp$mutrate2)
    
    if(optimization_prop !=1){
      genes_for_optimization = gene_prior[order(gene_prior[,2],decreasing = FALSE),1]
      genes_for_optimization = genes_for_optimization[1:floor(length(genes_for_optimization)*optimization_prop)]
      coverage = coverage[is.element(coverage$genename, genes_for_optimization),]
    }
    
    fr<-function(x){ # the function that will be optimized later for promoter_beta, enhancer_beta
      # the order of relative risk in x is epigenomic_relative_risk followed by anno_relative_risk
      all_rr = x[1:num_addition_par] # relative risks for features other than splicing
      splicing_rr = x[(num_addition_par+1) : (num_addition_par+1)]
      logP_Zg1 = by(coverage, coverage[,"genename"],  
                    function(x) sum(x$mut_count*(log(x$adjusted_mutrate)+(x$promoter+x$nf)*(as.matrix(x[,12:(12+num_addition_par-1)])%*%all_rr))-x$adjusted_mutrate*exp((x$promoter+x$nf)*(as.matrix(x[,12:(12+num_addition_par-1)])%*%all_rr))-log(factorial(x$mut_count))))
      logP_Zg0 = by(coverage, coverage[,"genename"],
                    function(x) sum(x$mut_count*log(x$adjusted_mutrate)-x$adjusted_mutrate-log(factorial(x$mut_count))))
      
      zg1_splicing_logp <- function(x){
        splicing_mutate_per_gene = sum(x[x$coding == 1,]$adjusted_mutrate)*splicing_num/syn_num # use coding mutation (syn) rate for each as the reference to get splicing mutation rate
        splicing_count_per_gene = sum(x$raw_splicing_mut)
        splicing_count_per_gene*(log(splicing_mutate_per_gene)+splicing_rr)-splicing_mutate_per_gene*exp(splicing_rr)-log(factorial(splicing_count_per_gene))
      }
      
      zg0_splicing_logp <- function(x){
        splicing_mutate_per_gene = sum(x[x$coding == 1,]$adjusted_mutrate)*splicing_num/syn_num # use coding mutation (syn) rate for each as the reference to get splicing mutation rate
        splicing_count_per_gene = sum(x$raw_splicing_mut)
        splicing_count_per_gene*log(splicing_mutate_per_gene) - splicing_mutate_per_gene - log(factorial(splicing_count_per_gene))
      }
      
      splicing_logP_Zg1 = by(coverage, coverage[,"genename"], zg1_splicing_logp)
      splicing_logP_Zg0 = by(coverage, coverage[,"genename"], zg0_splicing_logp)
      
      gene_prob_table = data.frame(genename = names(logP_Zg1), Zg1 = as.vector(logP_Zg1), Zg0 = as.vector(logP_Zg0), splicing_Zg1 = as.vector(splicing_logP_Zg1), splicing_Zg0 = as.vector(splicing_logP_Zg0))
      gene_prior = read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      gene_prior = data.frame(genename = gene_prior[,1], prob = 1-gene_prior[,2])
      gene_mle_table = merge(gene_prob_table, gene_prior, by.x = "genename", by.y = "genename")
      
      gene_mle_table[is.na(gene_mle_table$splicing_Zg1),]$splicing_Zg1 = 0 # there are about 35 genes that don't have coding windows which will give NA values, change logP to 0 for these cases
      gene_mle_table[is.na(gene_mle_table$splicing_Zg0),]$splicing_Zg0 = 0 # there are about 35 genes that don't have coding windows which will give NA values, change logP to 0 for these cases
      
      sum(log((gene_mle_table$prob*exp(gene_mle_table$Zg1+gene_mle_table$splicing_Zg1)+(1-gene_mle_table$prob)*exp(gene_mle_table$Zg0+gene_mle_table$splicing_Zg0)))) # minimization
    }
    mle = optim(rep(0.1, num_addition_par+1), fr,control=list("fnscale"=-1), hessian = TRUE)
  }
  system(paste("rm ", prefix, "_temp*", sep = ""))
  list(mle = mle, regression_coef = summary(out.offset)$coefficients)
}
