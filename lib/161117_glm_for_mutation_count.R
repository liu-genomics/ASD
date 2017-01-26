# function to run glm regression for window based mutation count
run_glm_for_mutation <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks = "no", overlap = 0.5, rm_nonsyn = FALSE, annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", annovar_input = "no.txt"){
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
  if(epigenomic_marks == "no"){
    out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
  }
  else{
    command = paste("bedtools intersect -f ", overlap, " -a ", window_file , " -b ", epigenomic_marks, " -wa > ", paste(prefix,"_temp_overlap_epi.bed", sep = ""),sep = "")
    system(command)
    window_in_epi = read.delim(paste(prefix,"_temp_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    window_in_epi = data.frame(site_index = window_in_epi[,4], epi = 1)
    coverage = merge(coverage, window_in_epi, by.x = "site_index", by.y = "site_index", all.x = TRUE)
    coverage[is.na(coverage$epi),]$epi = 0
    out.offset <- glm(mut_count ~ coding+promoter+cg+epi+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
  }
  system(paste("rm ", prefix, "_temp*", sep = ""))
  out.offset
}


# the v2 run glm_for_mutation functions allows for adding multiple number of epigenomics dataset and window-level scores. 
run_glm_for_mutation_v2 <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks_list = "no", overlap = 0.5, gene_assign_file,rm_nonsyn = FALSE, annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", annovar_input = "no.txt", gene_prior_file, 
                                    sequence_annotation_list = "no"){
  # [mut_file] is the mutation location file from one dataset. e.g.,../analysis/160229_data_for_analysis.bed, chr \t start \t end \t index,  1-based start and end
  # [window_file] is the file with genome partitioned into windows. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed
  # [cg_file] is the file that has CG-content calculated, matching [window_file]. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.cg
  # [mutrate_file] is the file that has mutrate data. e.g., "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate"
  # [sample_size] is the number of individuals
  # [epitenomic_marks_list] if it is "no", then glm won't use if a mutation is in h3k27ac as a categorical variable, if it is a vector of epigenomic annotation names, epigenomic mark will be added one by one.
  # for example, c("../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed","../other_annotation/epigenomic_annotation/Epigenome_E081_E082_intersection.bed","../other_annotation/epigenomic_annotation/encode_DHS_union.bed"). It is fine to directly use original epigenomics bed files
  # instead of the partitioned files. Because the window_file already has annotations for each window in terms of whetehr they are promoter or coding. 
  # [overlap] the fraction of overlap (the paramter f in bedtools intersect), default is 0.5, For the ease of TADA modeling, it is better to assign it to be 1E-9(1bp)
  # [rm_nonsyn] TRUE or FALSE. Indicating whether nonysnonymous mutations should be removed when fitting the model. Use TRUE when fitting model using ASD mutations. 
  # [annovar_folder] is where the annovar program is, e.g., . e.g., /media/yuwen/Elements/ANNOVAR/annovar/
  # [annovar_input] is a input file generated from mutations that will be annotated by annovar. The format is chr \t start \t end \t ref_allele \t alt_allele \t index 1-based for starts and ends
  # [sequence_annotation_list] a list of annotation files for each window in the window_file, e.g., after extracting conservation score for each window. c("../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.phastcons100way")
  
  
  
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
  if(epigenomic_marks_list == "no" & sequence_annotation_list == "no"){
    out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
  }
  else if(epigenomic_marks_list != "no" & sequence_annotation_list == "no"){
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
    num_addition_par = length(epigenomic_marks_list)# number of additional feature parameters 
    f = paste("mut_count ~ coding + promoter + cg + ",paste(tail(colnames(coverage), num_addition_par), collapse = " + "), "+ offset(log(2*mutrate*sample_size))", sep = "")
    out.offset <- do.call("glm", list(as.formula(f), data = coverage, family = poisson))
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
    f = paste("mut_count ~ coding + promoter + cg + ",paste(tail(colnames(coverage), num_addition_par), collapse = " + "), "+ offset(log(2*mutrate*sample_size))", sep = "")
    out.offset <- do.call("glm", list(as.formula(f), data = coverage, family = poisson))
  }
  
  system(paste("rm ", prefix, "_temp*", sep = ""))
  out.offset
  summary(out.offset)$coefficients
}


# simulate mutation count from real gene and histone modificationd data
simulate_noncoding_mutations <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks = "no", overlap = 0.5, gene_assign_file,rm_nonsyn = FALSE, annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", annovar_input = "no.txt", gene_prior_file,
                                         promoter_effect=2, enhancer_effect=1.5, optimization_prop = 1){
  # [mut_file] is the mutation location file from one dataset. e.g.,../analysis/160229_data_for_analysis.bed, chr \t start \t end \t index,  1-based start and end
  # [window_file] is the file with genome partitioned into windows. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed
  # [cg_file] is the file that has CG-content calculated, matching [window_file]. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.cg
  # [mutrate_file] is the file that has mutrate data. e.g., "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate"
  # [sample_size] is the number of individuals
  # [epitenomic_marks] if it is "no", then glm won't use if a mutation is in h3k27ac as a categorical variable, if it is a epigenomic annotation name.
  # for example, "../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed". It is fine to directly use original epigenomics bed files
  # instead of the partitioned files. Because the window_file already has annotations for each window in terms of whetehr they are promoter or coding. 
  # [overlap] the fraction of overlap (the paramter f in bedtools intersect), default is 0.5, For the ease of TADA modeling, it is better to assign it to be 1E-9(1bp)
  # [gene_assign_file] is the file that assigns each window to a gene. for example ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename.bed
  # [rm_nonsyn] TRUE or FALSE. Indicating whether nonysnonymous mutations should be removed when fitting the model. Use TRUE when fitting model using ASD mutations. 
  # [annovar_folder] is where the annovar program is, e.g., . e.g., /media/yuwen/Elements/ANNOVAR/annovar/
  # [annovar_input] is a input file generated from mutations that will be annotated by annovar. The format is chr \t start \t end \t ref_allele \t alt_allele \t index 1-based for starts and ends
  # [gene_prior_file], a file that has prior (derived from posterior and prior)for a gene as a risk gene. genename \t probability e.g. ../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.05.txt
  # [promoter_effect] the relative risk of mutations in active promoters
  # [enhancer_effect] the relative risk of mutations in active enhancers
  #[optimization_prop], the top percentage of genes that are used for optimization. genes are ranked high to low based on their Posterior probability of the alternative hypothesis 
  
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
  if(epigenomic_marks == "no"){
    out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
  }
  else{
    command = paste("bedtools intersect -f ", overlap, " -a ", window_file , " -b ", epigenomic_marks, " -wa > ", paste(prefix,"_temp_overlap_epi.bed", sep = ""),sep = "")
    system(command)
    window_in_epi = read.delim(paste(prefix,"_temp_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    window_in_epi = data.frame(site_index = window_in_epi[,4], epi = 1)
    coverage = merge(coverage, window_in_epi, by.x = "site_index", by.y = "site_index", all.x = TRUE)
    coverage[is.na(coverage$epi),]$epi = 0
    gene_assign = read.delim(gene_assign_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(gene_assign) = c("site_index", "genename")
    out.offset <- glm(mut_count ~ coding+promoter+cg+epi+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
    coverage = data.frame(coverage, adjusted_mutrate = out.offset$fitted.values)
    mutrate_temp = merge(coverage[,c("site_index","mut_count","epi","coding", "promoter","nf","adjusted_mutrate")], gene_assign, by.x = "site_index", by.y = "site_index")
    # now I need to use the prior probability of each gene to generate risk gene indicator and mutation numbers
    gene_prior = read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    gene_risk_assign = mapply(function(x){rbinom(1,1,1-x)}, gene_prior[,2])
    gene_prior = data.frame(gene_prior, risk_indicator = gene_risk_assign)
    mutrate_temp_risk = mutrate_temp[is.element(mutrate_temp$genename, gene_prior[gene_prior$risk_indicator == 1,1]),]
    mutrate_temp_nonrisk = mutrate_temp[!is.element(mutrate_temp$genename, gene_prior[gene_prior$risk_indicator == 1,1]),]
    mutrate_temp_risk = data.frame(mutrate_temp_risk, mutrate2 = mutrate_temp_risk$adjusted_mutrate*exp(mutrate_temp_risk$epi*mutrate_temp_risk$promoter*log(promoter_effect)+mutrate_temp_risk$epi*mutrate_temp_risk$nf*log(enhancer_effect)))
    mutrate_temp_nonrisk = data.frame(mutrate_temp_nonrisk, mutrate2 = mutrate_temp_nonrisk$adjusted_mutrate)
    mutrate_temp = rbind(mutrate_temp_risk,mutrate_temp_nonrisk)
    mutrate_temp = data.frame(mutrate_temp, sim_mut_count = mapply(function(x){rpois(1,x)}, mutrate_temp$mutrate2))
    if(optimization_prop !=1){
      genes_for_optimization = gene_prior[order(gene_prior[,2],decreasing = FALSE),1]
      genes_for_optimization = genes_for_optimization[1:floor(length(genes_for_optimization)*optimization_prop)]
      mutrate_temp = mutrate_temp[is.element(mutrate_temp$genename, genes_for_optimization),]
    }
    #test_count = mapply(function(x){rpois(1,x)},mutrate_temp$mutrate2)
    fr<-function(x){ # the function that will be optimized later for promoter_beta, enhancer_beta
      promoter_beta = x[1]
      enhancer_beta = x[2]
      logP_Zg1 = by(mutrate_temp[,c("sim_mut_count","coding","promoter","nf","epi","adjusted_mutrate")], mutrate_temp[,"genename"],
                    function(x) sum(x$sim_mut_count*(log(x$adjusted_mutrate)+x$promoter*x$epi*promoter_beta+x$nf*x$epi*enhancer_beta)-x$adjusted_mutrate*exp(x$promoter*x$epi*promoter_beta+x$nf*x$epi*enhancer_beta)-log(factorial(x$sim_mut_count))))
      logP_Zg0 = by(mutrate_temp[,c("sim_mut_count","coding","promoter","nf","epi","adjusted_mutrate")], mutrate_temp[,"genename"],
                    function(x) sum(x$sim_mut_count*log(x$adjusted_mutrate)-x$adjusted_mutrate-log(factorial(x$sim_mut_count))))
      gene_prob_table = data.frame(genename = names(logP_Zg1), Zg1 = as.vector(logP_Zg1), Zg0 = as.vector(logP_Zg0))
      gene_prior = read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      gene_prior = data.frame(genename = gene_prior[,1], prob = 1-gene_prior[,2])
      gene_mle_table = merge(gene_prob_table, gene_prior, by.x = "genename", by.y = "genename")
      sum(log((gene_mle_table$prob*exp(gene_mle_table$Zg1)+(1-gene_mle_table$prob)*exp(gene_mle_table$Zg0)))) # minimization
    }
    mle = optim(c(0.1,0.1), fr, control=list("fnscale"=-1), hessian = TRUE)
  }
  system(paste("rm ", prefix, "_temp*", sep = ""))
  list(mutrate_rate_table = mutrate_temp, mle = mle, true_effect = c(promoter_effect, enhancer_effect))
  #list(promoter = active_promoter, enhancer = active_enhancer, coding = coding_part)
}




# use verified parameter estimate methods (code) to estiamte effect sizes of active promoters and active enhancers. 
verified_effect_size_estimate_noncoding_mutations <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks = "no", overlap = 0.5, gene_assign_file,rm_nonsyn = FALSE, annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", annovar_input = "no.txt", gene_prior_file){
  # [mut_file] is the mutation location file from one dataset. e.g.,../analysis/160229_data_for_analysis.bed, chr \t start \t end \t index,  1-based start and end
  # [window_file] is the file with genome partitioned into windows. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed
  # [cg_file] is the file that has CG-content calculated, matching [window_file]. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.cg
  # [mutrate_file] is the file that has mutrate data. e.g., "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate"
  # [sample_size] is the number of individuals
  # [epitenomic_marks] if it is "no", then glm won't use if a mutation is in h3k27ac as a categorical variable, if it is a epigenomic annotation name.
  # for example, "../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed". It is fine to directly use original epigenomics bed files
  # instead of the partitioned files. Because the window_file already has annotations for each window in terms of whetehr they are promoter or coding. 
  # [overlap] the fraction of overlap (the paramter f in bedtools intersect), default is 0.5, For the ease of TADA modeling, it is better to assign it to be 1E-9(1bp)
  # [gene_assign_file] is the file that assigns each window to a gene. for example ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename.bed
  # [rm_nonsyn] TRUE or FALSE. Indicating whether nonysnonymous mutations should be removed when fitting the model. Use TRUE when fitting model using ASD mutations. 
  # [annovar_folder] is where the annovar program is, e.g., . e.g., /media/yuwen/Elements/ANNOVAR/annovar/
  # [annovar_input] is a input file generated from mutations that will be annotated by annovar. The format is chr \t start \t end \t ref_allele \t alt_allele \t index 1-based for starts and ends
  # [gene_prior_file], a file that has prior (derived from posterior and prior)for a gene as a risk gene. genename \t probability e.g. ../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.05.txt
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
  if(epigenomic_marks == "no"){
    out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
  }
  else{
    command = paste("bedtools intersect -f ", overlap, " -a ", window_file , " -b ", epigenomic_marks, " -wa > ", paste(prefix,"_temp_overlap_epi.bed", sep = ""),sep = "")
    system(command)
    window_in_epi = read.delim(paste(prefix,"_temp_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    window_in_epi = data.frame(site_index = window_in_epi[,4], epi = 1)
    coverage = merge(coverage, window_in_epi, by.x = "site_index", by.y = "site_index", all.x = TRUE)
    coverage[is.na(coverage$epi),]$epi = 0
    gene_assign = read.delim(gene_assign_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(gene_assign) = c("site_index", "genename")
    out.offset <- glm(mut_count ~ coding+promoter+cg+epi+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
    coverage = data.frame(coverage, adjusted_mutrate = out.offset$fitted.values)
    coverage = merge(coverage[,c("site_index","mut_count","epi","coding", "promoter","nf","adjusted_mutrate")], gene_assign, by.x = "site_index", by.y = "site_index")
    # now I need to use the prior probability of each gene to generate risk gene indicator and mutation numbers
    gene_prior = read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    #test_count = mapply(function(x){rpois(1,x)},mutrate_temp$mutrate2)
    fr<-function(x){ # the function that will be optimized later for promoter_beta, enhancer_beta
      promoter_beta = x[1]
      enhancer_beta = x[2]
      logP_Zg1 = by(coverage[,c("mut_count","coding","promoter","nf","epi","adjusted_mutrate")], coverage[,"genename"],
                    function(x) sum(x$mut_count*(log(x$adjusted_mutrate)+x$promoter*x$epi*promoter_beta+x$nf*x$epi*enhancer_beta)-x$adjusted_mutrate*exp(x$promoter*x$epi*promoter_beta+x$nf*x$epi*enhancer_beta)-log(factorial(x$mut_count))))
      logP_Zg0 = by(coverage[,c("mut_count","coding","promoter","nf","epi","adjusted_mutrate")], coverage[,"genename"],
                    function(x) sum(x$mut_count*log(x$adjusted_mutrate)-x$adjusted_mutrate-log(factorial(x$mut_count))))
      gene_prob_table = data.frame(genename = names(logP_Zg1), Zg1 = as.vector(logP_Zg1), Zg0 = as.vector(logP_Zg0))
      gene_prior = read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      gene_prior = data.frame(genename = gene_prior[,1], prob = 1-gene_prior[,2])
      gene_mle_table = merge(gene_prob_table, gene_prior, by.x = "genename", by.y = "genename")
      sum(log((gene_mle_table$prob*exp(gene_mle_table$Zg1)+(1-gene_mle_table$prob)*exp(gene_mle_table$Zg0)))) # minimization
    }
    mle = optim(c(0.1,0.1), fr,control=list("fnscale"=-1), hessian = TRUE)
  }
  system(paste("rm ", prefix, "_temp*", sep = ""))
  mle
  #list(promoter = active_promoter, enhancer = active_enhancer, coding = coding_part)
}




# function to run glm regression for window based mutation count, and then get adjusted mutation rate and mutation count for active promoters and enhancers per gene
predict_sum_mutation_rate_per_gene <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks = "no", overlap = 0.5, gene_assign_file,rm_nonsyn = FALSE, annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", annovar_input = "no.txt"){
  # [mut_file] is the mutation location file from one dataset. e.g.,../analysis/160229_data_for_analysis.bed, chr \t start \t end \t index,  1-based start and end
  # [window_file] is the file with genome partitioned into windows. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed
  # [cg_file] is the file that has CG-content calculated, matching [window_file]. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.cg
  # [mutrate_file] is the file that has mutrate data. e.g., "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate"
  # [sample_size] is the number of individuals
  # [epitenomic_marks] if it is "no", then glm won't use if a mutation is in h3k27ac as a categorical variable, if it is a epigenomic annotation name.
  # for example, "../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed". It is fine to directly use original epigenomics bed files
  # instead of the partitioned files. Because the window_file already has annotations for each window in terms of whetehr they are promoter or coding. 
  # [overlap] the fraction of overlap (the paramter f in bedtools intersect), default is 0.5, For the ease of TADA modeling, it is better to assign it to be 1E-9(1bp)
  # [gene_assign_file] is the file that assigns each window to a gene. for example ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename.bed
  # [rm_nonsyn] TRUE or FALSE. Indicating whether nonysnonymous mutations should be removed when fitting the model. Use TRUE when fitting model using ASD mutations. 
  # [annovar_folder] is where the annovar program is, e.g., . e.g., /media/yuwen/Elements/ANNOVAR/annovar/
  # [annovar_input] is a input file generated from mutations that will be annotated by annovar. The format is chr \t start \t end \t ref_allele \t alt_allele \t index 1-based for starts and ends
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
  if(epigenomic_marks == "no"){
    out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
  }
  else{
    command = paste("bedtools intersect -f ", overlap, " -a ", window_file , " -b ", epigenomic_marks, " -wa > ", paste(prefix,"_temp_overlap_epi.bed", sep = ""),sep = "")
    system(command)
    window_in_epi = read.delim(paste(prefix,"_temp_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    window_in_epi = data.frame(site_index = window_in_epi[,4], epi = 1)
    coverage = merge(coverage, window_in_epi, by.x = "site_index", by.y = "site_index", all.x = TRUE)
    coverage[is.na(coverage$epi),]$epi = 0
    gene_assign = read.delim(gene_assign_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(gene_assign) = c("site_index", "genename")
    out.offset <- glm(mut_count ~ coding+promoter+cg+epi+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
    coverage = data.frame(coverage, adjusted_mutrate = out.offset$fitted.values)
    mutrate_temp = merge(coverage[,c("site_index","mut_count","epi","coding", "promoter","nf","adjusted_mutrate")], gene_assign, by.x = "site_index", by.y = "site_index")
    temp = by(mutrate_temp[mutrate_temp$promoter ==1 & mutrate_temp$coding == 0  & mutrate_temp$nf ==0 & mutrate_temp$epi ==1 ,]$adjusted_mutrate, 
              mutrate_temp[mutrate_temp$promoter ==1 & mutrate_temp$coding == 0  & mutrate_temp$nf ==0 & mutrate_temp$epi ==1 ,]$genename, sum)
    active_promoter = data.frame(genename = names(temp),adjusted_mutrate_sum = as.vector(temp))
    temp = by(mutrate_temp[mutrate_temp$promoter ==1 & mutrate_temp$coding == 0  & mutrate_temp$nf ==0 & mutrate_temp$epi ==1 ,]$mut_count, 
              mutrate_temp[mutrate_temp$promoter ==1 & mutrate_temp$coding == 0  & mutrate_temp$nf ==0 & mutrate_temp$epi ==1 ,]$genename, sum)
    active_promoter = data.frame(active_promoter, mut_count_sum = as.vector(temp))
    temp = by(mutrate_temp[mutrate_temp$promoter ==0 & mutrate_temp$coding == 0  & mutrate_temp$nf ==1 & mutrate_temp$epi ==1 ,]$adjusted_mutrate, 
              mutrate_temp[mutrate_temp$promoter ==0 & mutrate_temp$coding == 0  & mutrate_temp$nf ==1 & mutrate_temp$epi ==1 ,]$genename, sum)
    active_enhancer = data.frame(genename = names(temp),adjusted_mutrate_sum = as.vector(temp))
    temp = by(mutrate_temp[mutrate_temp$promoter ==0 & mutrate_temp$coding == 0  & mutrate_temp$nf ==1 & mutrate_temp$epi ==1 ,]$mut_count, 
              mutrate_temp[mutrate_temp$promoter ==0 & mutrate_temp$coding == 0  & mutrate_temp$nf ==1 & mutrate_temp$epi ==1 ,]$genename, sum)
    active_enhancer = data.frame(active_enhancer, mut_count_sum = as.vector(temp))
    temp = by(mutrate_temp[mutrate_temp$promoter ==0 & mutrate_temp$coding == 1  & mutrate_temp$nf ==0,]$adjusted_mutrate, 
              mutrate_temp[mutrate_temp$promoter ==0 & mutrate_temp$coding == 1  & mutrate_temp$nf ==0,]$genename, sum)
    coding_part = data.frame(genename = names(temp),adjusted_mutrate_sum = as.vector(temp))
    temp = by(mutrate_temp[mutrate_temp$promoter ==0 & mutrate_temp$coding == 1  & mutrate_temp$nf ==0,]$mut_count, 
              mutrate_temp[mutrate_temp$promoter ==0 & mutrate_temp$coding == 1  & mutrate_temp$nf ==0,]$genename, sum)
    coding_part = data.frame(coding_part, mut_count_sum = as.vector(temp))
    rm(temp)
  }
  system(paste("rm ", prefix, "_temp*", sep = ""))
  list(promoter = active_promoter, enhancer = active_enhancer, coding = coding_part)
}


# function to calculate log-likelihood from estimated effect sizes. 
get_10fold_log_likelihood <- function(prediction_model,fold_n = 10, sample_size){
  # [prediction_model] is the return value from function{run_glm_for_mutation}
  # [fold_n] is the number of partition of mutation data. 
  # [sample_size] is the number of samples that are used to train the prediction_model
  all_windows = data.frame(prediction_model$data, adjusted_mutrate = prediction_model$fitted.values, plain_mutrate = prediction_model$data$mutrate*sum(prediction_model$data$mut_count)/sum(prediction_model$data$mutrate))
  # in the above data.frame, mutrate is the tri-nucleotide based mutation rate without adjusting for total sample size and total number of mutations in considered regions (within 10kb of TSS and coding)
  all_windows = all_windows[sample(nrow(all_windows), nrow(all_windows)),] # permute the ordering of rows. 
  sub_size = ceiling(nrow(all_windows)/fold_n)
  output = data.frame(model1 = numeric(),model_null = numeric())
  for(i in 1:fold_n){
    sub_data = all_windows[((i-1)*sub_size+1):min(i*sub_size,nrow(all_windows)),]
    a = sum(mapply(function(x,y){x*log(y)-y-factorial(x)}, sub_data$mut_count, sub_data$adjusted_mutrate))
    b = sum(mapply(function(x,y){x*log(y)-y-factorial(x)}, sub_data$mut_count, sub_data$plain_mutrate))
    output[i,] = c(a,b)
  }
  output
}

# function to run glm regression for window based mutation count, get the effect size of active enhancer and active promoter
# a mixture model for each gene, probability for each gene based on multiplying probabilities across windows (using corrected mutation rates), rather than get the corrected mutation rate per gene
# needs to set overlap to be 1e-9 in most of the cases. 
estimate_effect_size <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks = "no", overlap = 0.5, gene_assign_file, gene_prior_file){
  # [mut_file] is the mutation location file from one dataset. e.g.,../analysis/160229_data_for_analysis.bed, chr \t start \t end \t index,  1-based start and end
  # [window_file] is the file with genome partitioned into windows. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed
  # [cg_file] is the file that has CG-content calculated, matching [window_file]. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.cg
  # [mutrate_file] is the file that has mutrate data. e.g., "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate"
  # [sample_size] is the number of individuals
  # [epitenomic_marks] if it is "no", then glm won't use if a mutation is in h3k27ac as a categorical variable, if it is a epigenomic annotation name.
  # for example, "../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed". It is fine to directly use original epigenomics bed files
  # instead of the partitioned files. Because the window_file already has annotations for each window in terms of whetehr they are promoter or coding. 
  # [overlap] the fraction of overlap (the paramter f in bedtools intersect), default is 0.5, For the ease of TADA modeling, it is better to assign it to be 1E-9(1bp)
  # [gene_assign_file] is the file that assigns each window to a gene. for example ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename.bed
  # [gene_prior_file], a file that has prior probability(FDR) for a gene as a risk gene. genename \t probability e.g. ../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_combined_qvalue.txt
  prefix=system("date +%s", intern = TRUE) # prefix for temporary files that will be deleted at the end of the pipeline
  mut = read.delim(mut_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
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
  if(epigenomic_marks == "no"){
    out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
  }
  else{
    command = paste("bedtools intersect -f ", overlap, " -a ", window_file , " -b ", epigenomic_marks, " -wa > ", paste(prefix,"_temp_overlap_epi.bed", sep = ""),sep = "")
    system(command)
    window_in_epi = read.delim(paste(prefix,"_temp_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    window_in_epi = data.frame(site_index = window_in_epi[,4], epi = 1)
    coverage = merge(coverage, window_in_epi, by.x = "site_index", by.y = "site_index", all.x = TRUE)
    coverage[is.na(coverage$epi),]$epi = 0
    gene_assign = read.delim(gene_assign_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(gene_assign) = c("site_index", "genename")
    out.offset <- glm(mut_count ~ coding+promoter+cg+epi+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
    coverage = data.frame(coverage, adjusted_mutrate = out.offset$fitted.values)
    coverage = merge(coverage, gene_assign, by.x = "site_index", by.y = "site_index")
    fr<-function(x){ # the function that will be optimized later for promoter_beta, enhancer_beta
      promoter_beta = x[1]
      enhancer_beta = x[2]
      logP_Zg1 = by(coverage[,c("mut_count","coding","promoter","nf","epi","adjusted_mutrate")], coverage[,"genename"],
         function(x) sum(x$mut_count*(log(x$adjusted_mutrate)+x$promoter*x$epi*promoter_beta+x$nf*x$epi*enhancer_beta)-x$adjusted_mutrate*exp(x$promoter*x$epi*promoter_beta+x$nf*x$epi*enhancer_beta)-log(factorial(x$mut_count))))
      logP_Zg0 = by(coverage[,c("mut_count","coding","promoter","nf","epi","adjusted_mutrate")], coverage[,"genename"],
         function(x) sum(x$mut_count*log(x$adjusted_mutrate)-x$adjusted_mutrate-log(factorial(x$mut_count))))
      gene_prob_table = data.frame(genename = names(logP_Zg1), Zg1 = as.vector(logP_Zg1), Zg0 = as.vector(logP_Zg0))
      gene_prior = read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      gene_prior = data.frame(genename = gene_prior[,1], prob = 1-gene_prior[,2])
      gene_mle_table = merge(gene_prob_table, gene_prior, by.x = "genename", by.y = "genename")
      sum(log((gene_mle_table$prob*exp(gene_mle_table$Zg1)+(1-gene_mle_table$prob)*exp(gene_mle_table$Zg0))))
    }
    mle = optim(c(0.1,0.1), fr)
  }
  system(paste("rm ", prefix, "_temp*", sep = ""))
  mle
}


# function to get the number of coding mutations from window files that I have used to train mutation rate models. 
get_coding_mut_number<-function(mut_info, ref_alt_allele, window_file, annovar_folder, phenotype, geneset = "no"){
  # [mut_info] is the mutation information in the data.frame {mutation} in Rdata {*transformed_for_old_code.Rdata}
  # [ref_alt_allele] is the information of ref/mutation alleles stored as a data.frame in {ref_alt_allele} in Rdata {*transformed_for_old_code.Rdata}
  # [window_file] is the file with genome partitioned into windows. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed
  # [annovar_folder] is the folder that has annovar.pl. e.g., /media/yuwen/Elements/ANNOVAR/annovar/
  # [geneset] if "no", then mutation counts of different types for all the genes will be returned. Otherwise, could specify strigent_ASD, relaxed_ASD or other genesets. 
  # use Annovar to annotate coding mutations. 
  # [phenotype] is the phenotype that I want to count, for example "ASD" or "control"
  
  prefix=system("date +%s", intern = TRUE) # prefix for temporary files that will be deleted at the end of the pipeline
  temp_for_annovar = data.frame(mut_info[,1:3], ref_alt_allele, mut_info$mut_type, mut_info$index, mut_info$phenotype)
  temp_for_annovar = temp_for_annovar[temp_for_annovar$mut_info.phenotype == phenotype,1:7]
  write.table(temp_for_annovar, paste(prefix, "_for_annovar_temp.txt",sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  
  # only keep mutations that overlap with coding portions from [window_file]
  command = paste("grep coding ", window_file, " > ", paste(prefix, "_coding_window_temp.bed", sep = ""), sep = "")
  system(command)
  command = paste("bedtools intersect -a ",paste(prefix, "_for_annovar_temp.txt",sep = ""), " -b ", paste(prefix, "_coding_window_temp.bed", sep = ""), " -wa > ",paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""), sep = "")
  system(command)
  command = paste("perl ", annovar_folder, "annotate_variation.pl -geneanno -buildver hg19 ", paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""), " ", annovar_folder, "humandb/", sep = "")
  system(command)
  coding_anno = read.delim(paste(paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""),".exonic_variant_function", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  coding_gene = unlist(lapply(strsplit(coding_anno[,3],split = ":"),`[[`,1))
  coding_mut_for_gene = data.frame(mut_type_raw = coding_anno[,2], genename = coding_gene)
  coding_mut_for_gene = data.frame(coding_mut_for_gene, mut_type = "not_defined")
  coding_mut_for_gene$mut_type = as.character( coding_mut_for_gene$mut_type)
  coding_mut_for_gene[coding_mut_for_gene$mut_type_raw == "synonymous SNV",]$mut_type = "synonymous"
  coding_mut_for_gene[coding_mut_for_gene$mut_type_raw == "nonsynonymous SNV",]$mut_type = "nonsynonymous"
  coding_mut_for_gene[coding_mut_for_gene$mut_type_raw == "stopgain" | coding_mut_for_gene$mut_type_raw == "frameshift deletion"| coding_mut_for_gene$mut_type_raw == "frameshift insertion" ,]$mut_type = "LoF"
  
  if(geneset != "no"){
    coding_mut_for_gene = coding_mut_for_gene[is.element(coding_mut_for_gene$gene, geneset), ]
  }
  #exon_mut_type = rep("unknown", length(mutation[,1]))
  #exon_mut_type[is.element(mutation$index, coding_anno[coding_anno[,2] == "frameshift deletion",10])] = "frameshift deletion"
  #exon_mut_type[is.element(mutation$index, coding_anno[coding_anno[,2] == "frameshift insertion",10])] = "frameshift insertion"
  #exon_mut_type[is.element(mutation$index, coding_anno[coding_anno[,2] == "nonframeshift deletion",10])] = "nonframeshift deletion"
  #exon_mut_type[is.element(mutation$index, coding_anno[coding_anno[,2] == "nonsynonymous SNV",10])] = "nonsynonymous"
  #exon_mut_type[is.element(mutation$index, coding_anno[coding_anno[,2] == "stopgain",10])] = "stopgain"
  #exon_mut_type[is.element(mutation$index, coding_anno[coding_anno[,2] == "synonymous SNV",10])] = "synonymous"

  system(paste("rm ", prefix, "_for_annovar_temp.txt*",sep = ""))
  
  mutation_number = list(syn_number = nrow(coding_mut_for_gene[coding_mut_for_gene[,3] == "synonymous",]),
       nonsyn_number = nrow(coding_mut_for_gene[coding_mut_for_gene[,3] == "nonsynonymous",]),
       lof_number = nrow(coding_mut_for_gene[coding_mut_for_gene[,3] == "LoF",]))
  list(mutation_number = mutation_number, coding_mut_for_gene = coding_mut_for_gene)
}


# function to compare predicted counts and observed counts for each type of mutations in different gene sets for a mutational model
compare_pre_with_obs <- function(mutation_model_coding, coding_count, geneset = c("stringent_ASD", "relaxed_ASD", "nonASD_genes","constraint_union","all_reseq_genes"), mutation_type = c("synonymous", "nonsynonymous", "LoF"), 
                                 mutation_scailing = list(synonymous = 0.2, nonsynonymous = 0.5, LoF = 0.01)){
#[mutation_model_coding] is the object that has predicted coding rate from {161120_mutation_rate_calibration_sum_per_gene_with_Noonan_brain_roadmap_union_1bp_overlap.Rdata}$ASD_model_manuscript_mutrate$coding
#[coding_count] is the {coding_mut_for_gene} element of the list returned from function{get_coding_mut_number}
#[geneset] is the name of the genesets that were defined in "../lib/160930_screen_burden_functions_test_more_genelist_focused.R"  
#[mutation_type] is the type of mutations that need to be looked at. They sould be those that are defined in function{get_coding_mut_number}
#[mutation_scailing] is a list of scaling factor that should be multiplied to total coding mutation rate to get the mutation rate of a specific type. 
  report = data.frame(geneset = character(), mut_type = character(), predict = integer(), observed = integer(), pvalue = numeric(),stringsAsFactors=FALSE)
  m = 1
  for(gs in geneset){
    for(mt in mutation_type){
      observed = nrow(coding_count[coding_count$mut_type == mt & is.element(coding_count$genename, eval(parse(text = gs))),])
      predict = sum(mutation_model_coding[is.element(mutation_model_coding$genename, eval(parse(text = gs))),]$adjusted_mutrate_sum)*eval(parse(text = paste("mutation_scailing$",mt)))
      pvalue = ppois(observed-1, predict, lower.tail = FALSE)
      report[m,] = c(as.character(gs), as.character(mt), predict, observed, pvalue)
      m = m + 1
    }
  }
  report
}



# v2 is different from the original version in that it could take multiple epigenomic marks, and multiple window-based annotation scores, such phastcons
verified_effect_size_estimate_noncoding_mutations_v2 <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks_list = "no", overlap = 0.5, gene_assign_file,rm_nonsyn = FALSE, annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", annovar_input = "no.txt", gene_prior_file, 
                                                                 sequence_annotation_list = "no", optimization_prop = 1){
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
    f = paste("mut_count ~ coding + promoter + cg + ",paste(tail(colnames(coverage), num_addition_par), collapse = " + "), "+ offset(log(2*mutrate*sample_size))", sep = "")
    out.offset <- do.call("glm", list(as.formula(f), data = coverage, family = poisson))
    
    coverage = data.frame(coverage, adjusted_mutrate = out.offset$fitted.values)
    
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
      all_rr = x
      logP_Zg1 = by(coverage, coverage[,"genename"],  
                    function(x) sum(x$mut_count*(log(x$adjusted_mutrate)+(x$promoter+x$nf)*(as.matrix(x[,12:(12+num_addition_par-1)])%*%all_rr))-x$adjusted_mutrate*exp((x$promoter+x$nf)*(as.matrix(x[,12:(12+num_addition_par-1)])%*%all_rr))-log(factorial(x$mut_count))))
      logP_Zg0 = by(coverage, coverage[,"genename"],
                    function(x) sum(x$mut_count*log(x$adjusted_mutrate)-x$adjusted_mutrate-log(factorial(x$mut_count))))
      gene_prob_table = data.frame(genename = names(logP_Zg1), Zg1 = as.vector(logP_Zg1), Zg0 = as.vector(logP_Zg0))
      gene_prior = read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      gene_prior = data.frame(genename = gene_prior[,1], prob = 1-gene_prior[,2])
      gene_mle_table = merge(gene_prob_table, gene_prior, by.x = "genename", by.y = "genename")
      sum(log((gene_mle_table$prob*exp(gene_mle_table$Zg1)+(1-gene_mle_table$prob)*exp(gene_mle_table$Zg0)))) # minimization
    }
    mle = optim(rep(0.1, num_addition_par), fr,control=list("fnscale"=-1), hessian = TRUE)
  }
  system(paste("rm ", prefix, "_temp*", sep = ""))
  list(mle = mle, regression_coef = summary(out.offset)$coefficients)
}




# use verified parameter estimate methods (code) to estiamte effect sizes of active promoters and active enhancers. 
# The V3 version is different from that in the first version in that only a fixed set of parameters (cg, promote, coding, will be used to adjust for mutation rate. )
verified_effect_size_estimate_noncoding_mutations_v3 <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks_list = "no", overlap = 0.5, gene_assign_file,rm_nonsyn = FALSE, annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", annovar_input = "no.txt", gene_prior_file, 
                                                                 sequence_annotation_list = "no", optimization_prop = 1){
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
      all_rr = x
      logP_Zg1 = by(coverage, coverage[,"genename"],  
                    function(x) sum(x$mut_count*(log(x$adjusted_mutrate)+(x$promoter+x$nf)*(as.matrix(x[,12:(12+num_addition_par-1)])%*%all_rr))-x$adjusted_mutrate*exp((x$promoter+x$nf)*(as.matrix(x[,12:(12+num_addition_par-1)])%*%all_rr))-log(factorial(x$mut_count))))
      logP_Zg0 = by(coverage, coverage[,"genename"],
                    function(x) sum(x$mut_count*log(x$adjusted_mutrate)-x$adjusted_mutrate-log(factorial(x$mut_count))))
      gene_prob_table = data.frame(genename = names(logP_Zg1), Zg1 = as.vector(logP_Zg1), Zg0 = as.vector(logP_Zg0))
      gene_prior = read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      gene_prior = data.frame(genename = gene_prior[,1], prob = 1-gene_prior[,2])
      gene_mle_table = merge(gene_prob_table, gene_prior, by.x = "genename", by.y = "genename")
      sum(log((gene_mle_table$prob*exp(gene_mle_table$Zg1)+(1-gene_mle_table$prob)*exp(gene_mle_table$Zg0)))) # minimization
    }
    mle = optim(rep(0.1, num_addition_par), fr,control=list("fnscale"=-1), hessian = TRUE)
  }
  system(paste("rm ", prefix, "_temp*", sep = ""))
  list(mle = mle, regression_coef = summary(out.offset)$coefficients)
}


# use verified parameter estimate methods (code) to estiamte effect sizes of active promoters and active enhancers. 
# The v4 version is different from that in the first version in that only a fixed set of parameters (cg, promote, coding, will be used to adjust for mutation rate. )
# The v4 version also included splicing mutations in estimating relative risks.
verified_effect_size_estimate_noncoding_mutations_v4 <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks_list = "no", overlap = 0.5, gene_assign_file,rm_nonsyn = FALSE, annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", annovar_input = "no.txt", gene_prior_file, 
                                                                 sequence_annotation_list = "no", optimization_prop = 1, mut_spidex){
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
  # [mut_spidex] is the merge of data_matrix table and mut_spidex table sotored the Rdata that holds the mutation information. e.g., "../data/debug_region_list_073116_data_matrixtransformed_for_old_code.Rdata"
  
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
    splicing_mut = mut_spidex[mut_spidex$dpsi_zscore <= splicing_threshold & !is.na(mut_spidex$dpsi_zscore),c("Chrom","Start","End","ID")]
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

# simulate mutation count from real gene and histone modificationd data, would incoporate multiple epigenomics data set and sequence-level annotations if possible, also combine promoter and enhancers within 10-kb
simulate_noncoding_mutations_v2 <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks_list = "no", overlap = 0.5, gene_assign_file,rm_nonsyn = FALSE, annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", annovar_input = "no.txt", gene_prior_file, 
                                            epi_effect_size = c(log(3),log(2),log(1)), sequence_annotation_list = "no", anno_effect_size = c(0.03), optimization_prop = 1){
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
  # [epi_effect_size] is a vector of the relative risk of mutations with different epigenomic marks for risk genes. For now, promoter and enhancers 10-kb within TSS are combined. The length should be the same as epigenomic_marks
  # [sequence_annotation_list] a list of annotation files for each window in the window_file, e.g., after extracting conservation score for each window. c("../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.phastcons100way")
  # [anno_effect_size] is a vector of the relative risk of mutations in regons with certain annotations, e.g., conservation score.length should be the same as the length of [sequence_annotation_list]. The value here is on the log-scale 
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
  if(epigenomic_marks == "no"){
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
    f = paste("mut_count ~ coding + promoter + cg + ",paste(tail(colnames(coverage), num_addition_par), collapse = " + "), "+ offset(log(2*mutrate*sample_size))", sep = "")
    out.offset <- do.call("glm", list(as.formula(f), data = coverage, family = poisson))
    
    gene_assign = read.delim(gene_assign_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(gene_assign) = c("site_index", "genename")
    coverage = data.frame(coverage, adjusted_mutrate = out.offset$fitted.values)
    mutrate_temp = merge(coverage, gene_assign, by.x = "site_index", by.y = "site_index")
    # now I need to use the prior probability of each gene to generate risk gene indicator and mutation numbers
    gene_prior = read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    gene_risk_assign = mapply(function(x){rbinom(1,1,1-x)}, gene_prior[,2])
    gene_prior = data.frame(gene_prior, risk_indicator = gene_risk_assign)
    mutrate_temp_risk = mutrate_temp[is.element(mutrate_temp$genename, gene_prior[gene_prior$risk_indicator == 1,1]),]
    mutrate_temp_nonrisk = mutrate_temp[!is.element(mutrate_temp$genename, gene_prior[gene_prior$risk_indicator == 1,1]),]
    
    # for now assume the relative risk of different epigenomic marks add up independently. 
    # for simplicity, assume the effects of conservation scores also add up independently. 
    # assume the relative risk is the same for promoters and enhancers. 
    addition_par_rr = c(epi_effect_size, anno_effect_size)
    active_effect = as.matrix(mutrate_temp_risk[,12:(12+num_addition_par-1)])%*%addition_par_rr
    mutrate_temp_risk = data.frame(mutrate_temp_risk, mutrate2 = mutrate_temp_risk$adjusted_mutrate*exp((mutrate_temp_risk$promoter + mutrate_temp_risk$nf)*active_effect))
    mutrate_temp_nonrisk = data.frame(mutrate_temp_nonrisk, mutrate2 = mutrate_temp_nonrisk$adjusted_mutrate)
    mutrate_temp = rbind(mutrate_temp_risk,mutrate_temp_nonrisk)
    mutrate_temp = data.frame(mutrate_temp, sim_mut_count = mapply(function(x){rpois(1,x)}, mutrate_temp$mutrate2))
    if(optimization_prop !=1){
      genes_for_optimization = gene_prior[order(gene_prior[,2],decreasing = FALSE),1]
      genes_for_optimization = genes_for_optimization[1:floor(length(genes_for_optimization)*optimization_prop)]
      mutrate_temp = mutrate_temp[is.element(mutrate_temp$genename, genes_for_optimization),]
    }
    #test_count = mapply(function(x){rpois(1,x)},mutrate_temp$mutrate2)
    fr<-function(x){ # the function that will be optimized later for promoter_beta, enhancer_beta
      # the order of relative risk in x is epigenomic_relative_risk followed by anno_relative_risk
      all_rr = x
      logP_Zg1 = by(mutrate_temp, mutrate_temp[,"genename"],  
                    function(x) sum(x$sim_mut_count*(log(x$adjusted_mutrate)+(x$promoter+x$nf)*(as.matrix(x[,12:(12+num_addition_par-1)])%*%all_rr))-x$adjusted_mutrate*exp((x$promoter+x$nf)*(as.matrix(x[,12:(12+num_addition_par-1)])%*%all_rr))-log(factorial(x$sim_mut_count))))
      logP_Zg0 = by(mutrate_temp, mutrate_temp[,"genename"],
                    function(x) sum(x$sim_mut_count*log(x$adjusted_mutrate)-x$adjusted_mutrate-log(factorial(x$sim_mut_count))))
      gene_prob_table = data.frame(genename = names(logP_Zg1), Zg1 = as.vector(logP_Zg1), Zg0 = as.vector(logP_Zg0))
      gene_prior = read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      gene_prior = data.frame(genename = gene_prior[,1], prob = 1-gene_prior[,2])
      gene_mle_table = merge(gene_prob_table, gene_prior, by.x = "genename", by.y = "genename")
      sum(log((gene_mle_table$prob*exp(gene_mle_table$Zg1)+(1-gene_mle_table$prob)*exp(gene_mle_table$Zg0)))) # minimization
    }
    mle = optim(rep(0.1,num_addition_par), fr, control=list("fnscale"=-1), hessian = TRUE)
  }
  system(paste("rm ", prefix, "_temp*", sep = ""))
  list(mutrate_rate_table = mutrate_temp, mle = mle, true_effect = c(promoter_effect, enhancer_effect))
  #list(promoter = active_promoter, enhancer = active_enhancer, coding = coding_part)
}


# use estimated relative risks to get the bayes factors for each gene
bayes_factor_for_each_gene <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks_list = "no", overlap = 0.5, gene_assign_file,rm_nonsyn = FALSE, annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", annovar_input = "no.txt", gene_prior_file, 
                                       sequence_annotation_list = "no", mut_spidex, rr, coding_bayes_table_file, TADA_p0 = 0.94){
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
  # [mut_spidex] is the merge of data_matrix table and mut_spidex table sotored the Rdata that holds the mutation information. e.g., "../data/debug_region_list_073116_data_matrixtransformed_for_old_code.Rdata"
  # [rr] is the estimated risk. it is vector with relative risk estimates in order from epigenomic_annotations, sequence_annotations to splicing mutations.
  # [coding_bayes_table_file], a file that has Bayes factors calculated from previous coding studies. e.g., "../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_BF_for_de_novo_mutations.txt"
  # [TADA_prior] the prior proportion of genes being non-risk genes. Will used in TADA, default = 0.94
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
    num_addition_par = 0
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
    num_addition_par = length(epigenomic_marks_list) # number of additional feature parameters
    if(sequence_annotation_list != "no"){
      annotation_index = 1
      for(seq_annotation in sequence_annotation_list){
        annotation_raw = read.delim(seq_annotation, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        annotation = data.frame(site_index = annotation_raw[,1])
        annotation[,paste("annotation",annotation_index,sep = "")] = annotation_raw[,4]
        coverage = merge(coverage, annotation, by.x = "site_index", by.y = "site_index")
        annotation_index = annotation_index + 1
      }
      num_addition_par = length(epigenomic_marks_list) + length(sequence_annotation_list) # number of additional feature parameters
    }
    # only use a set of fixed predictors to run glm in order to adjust for mutation rates. 
    out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
    coverage = data.frame(coverage, adjusted_mutrate = out.offset$fitted.values)
    
    # get the number of synonymous mutations that are captured in the windows, which will be used to scale factor for splicing mutations, could be all mutations if rm_nonsy is set to be FALSE
    syn_num = sum(coverage[coverage$coding == 1 & coverage$mut_count ==1,]$coding)
    
    # get mutations that have dpsi_zscore under some certain cutoff, and these mutations will be defined as splicing mutations
    splicing_threshold = quantile(mut_spidex$dpsi_zscore,seq(0,1,0.1),na.rm = TRUE)[2] # lower 10% 
    splicing_mut = mut_spidex[mut_spidex$dpsi_zscore <= splicing_threshold & !is.na(mut_spidex$dpsi_zscore),c("Chrom","Start","End","ID")]
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
    
    coding_bayes_table = read.delim(coding_bayes_table_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    coding_bayes_table[,2] = log(coding_bayes_table[,2])
    colnames(coding_bayes_table) = c("genename", "logBF_coding")
    
    logBF_noncoding = by(coverage, coverage[,"genename"],
                         function(x) sum(x$mut_count*(x$promoter+x$nf)*(as.matrix(x[,12:(12+num_addition_par-1)])%*%rr[1:num_addition_par]) - x$adjusted_mutrate*(exp((x$promoter+x$nf)*(as.matrix(x[,12:(12+num_addition_par-1)])%*%rr[1:num_addition_par]))-1)))
    
    splicing_rr =rr[num_addition_par+1]
    logBF_splicing = by(coverage, coverage[,"genename"],
                        function(x) sum(x$raw_splicing_mut)*splicing_rr - sum(x[x$coding == 1,]$adjusted_mutrate)*splicing_num/syn_num*(exp(splicing_rr)-1))  
    
    
    gene_BF_table = data.frame(genename = names(logBF_noncoding), logBF_noncoding = as.vector(logBF_noncoding), logBF_splicing = as.vector(logBF_splicing))
    gene_BF_table = merge(gene_BF_table, coding_bayes_table, by = "genename")
    gene_BF_table$logBF_all = gene_BF_table$logBF_noncoding + gene_BF_table$logBF_splicing + gene_BF_table$logBF_coding
    gene_BF_table[,c("BF_noncoding","BF_splicing","BF_coding","BF_all")] = exp(gene_BF_table[,c("logBF_noncoding","logBF_splicing","logBF_coding","logBF_all")])
    gene_BF_table = gene_BF_table[order(gene_BF_table$BF_coding, decreasing = TRUE),]
    gene_BF_table$FDR_coding = Bayesian.FDR(gene_BF_table$BF_coding, TADA_p0)$FDR
    gene_BF_table = gene_BF_table[order(gene_BF_table$BF_all, decreasing = TRUE),]
    gene_BF_table$FDR_all = Bayesian.FDR(gene_BF_table$BF_all, TADA_p0)$FDR
  }
  system(paste("rm ", prefix, "_temp*", sep = ""))
  gene_BF_table
}


# The Version 2 of the function used all splicing mutations, including those mutations that are not overlapped by the window_file.
bayes_factor_for_each_gene_v2 <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks_list = "no", overlap = 0.5, gene_assign_file,rm_nonsyn = FALSE, annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/", annovar_input = "no.txt", gene_prior_file, 
                                          sequence_annotation_list = "no", mut_spidex, rr, coding_bayes_table_file, TADA_p0 = 0.94){
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
  # [mut_spidex] is the proband_SNV_mut_spedix in the R object. e.g., "../data/debug_region_list_073116_data_matrixtransformed_for_old_code.Rdata"
  # [rr] is the estimated risk. it is vector with relative risk estimates in order from epigenomic_annotations, sequence_annotations to splicing mutations. Remeber to use only the mutations that are either in cases or controls depending on the needs. 
  # The papramter is in the log scale. For example, if the real rr is 3, the parameter that needs to be set here is log(3).
  # [coding_bayes_table_file], a file that has Bayes factors calculated from previous coding studies. e.g., "../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_BF_for_de_novo_mutations.txt"
  # [TADA_prior] the prior proportion of genes being non-risk genes. Will used in TADA, default = 0.94
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
    num_addition_par = 0
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
    num_addition_par = length(epigenomic_marks_list) # number of additional feature parameters
    if(sequence_annotation_list != "no"){
      annotation_index = 1
      for(seq_annotation in sequence_annotation_list){
        annotation_raw = read.delim(seq_annotation, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        annotation = data.frame(site_index = annotation_raw[,1])
        annotation[,paste("annotation",annotation_index,sep = "")] = annotation_raw[,4]
        coverage = merge(coverage, annotation, by.x = "site_index", by.y = "site_index")
        annotation_index = annotation_index + 1
      }
      num_addition_par = length(epigenomic_marks_list) + length(sequence_annotation_list) # number of additional feature parameters
    }
    # only use a set of fixed predictors to run glm in order to adjust for mutation rates. 
    out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(2*mutrate*sample_size)), family = poisson, data = coverage)
    coverage = data.frame(coverage, adjusted_mutrate = out.offset$fitted.values)
    
    # get the number of synonymous mutations that are captured in the windows, which will be used to scale factor for splicing mutations, could be all mutations if rm_nonsy is set to be FALSE
    syn_num = sum(coverage[coverage$coding == 1 & coverage$mut_count ==1,]$coding)
    
    # get mutations that have dpsi_zscore under some certain cutoff, and these mutations will be defined as splicing mutations
    splicing_threshold = quantile(mut_spidex$dpsi_zscore,seq(0,1,0.1),na.rm = TRUE)[2] # lower 10% 
    splicing_mut = mut_spidex[mut_spidex$dpsi_zscore <= splicing_threshold & !is.na(mut_spidex$dpsi_zscore),]
    splicing_mut = as.data.frame(table(splicing_mut$gene))
    colnames(splicing_mut) = c("genename","raw_splicing_mut")
    
    #get splicing mutation number
    # splicing_num = sum(splicing_mut$raw_splicing_mut) # only get 51 for ASD_manuscript_cases
    
    
    gene_assign = read.delim(gene_assign_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(gene_assign) = c("site_index", "genename")
    
    coverage = merge(coverage, gene_assign, by.x = "site_index", by.y = "site_index")
    # now I need to use the prior probability of each gene to generate risk gene indicator and mutation numbers
    gene_prior = read.delim(gene_prior_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    #test_count = mapply(function(x){rpois(1,x)},mutrate_temp$mutrate2)
    
    coding_bayes_table = read.delim(coding_bayes_table_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    coding_bayes_table[,2] = log(coding_bayes_table[,2])
    colnames(coding_bayes_table) = c("genename", "logBF_coding")
    
    logBF_noncoding = by(coverage, coverage[,"genename"],
                         function(x) sum(x$mut_count*(x$promoter+x$nf)*(as.matrix(x[,12:(12+num_addition_par-1)])%*%rr[1:num_addition_par]) - x$adjusted_mutrate*(exp((x$promoter+x$nf)*(as.matrix(x[,12:(12+num_addition_par-1)])%*%rr[1:num_addition_par]))-1)))
    gene_with_epi = by(coverage, coverage[,"genename"],
                       function(x) sum(apply(x[,c(9,10)]*x[,12:(12+length(epigenomic_marks_list)-1)],1,sum)))
    
    genes_with_no_epi = names(gene_with_epi[gene_with_epi == 0])
    
    splicing_rr =rr[num_addition_par+1]
    
    splicing_mutrate = by(coverage, coverage[,"genename"],
                          function(x) sum(x[x$coding == 1,]$adjusted_mutrate))
    splicing_mutrate = data.frame(genename = names(splicing_mutrate), adjusted_mutrate = as.vector(splicing_mutrate))
    splicing_data = merge(splicing_mutrate, splicing_mut, by = "genename", all.x = TRUE)
    splicing_data[is.na(splicing_data$raw_splicing_mut),]$raw_splicing_mut = 0
    
    #scale to get the splicing mutation rate for each gene
    splicing_data$adjusted_mutrate = splicing_data$adjusted_mutrate*sum(splicing_data$raw_splicing_mut)/syn_num
    
    logBF_splicing = by(splicing_data, splicing_data[,"genename"],
                        function(x) (sum(x$raw_splicing_mut)*splicing_rr - sum(x$adjusted_mutrate)*(exp(splicing_rr)-1)))  
    logBF_splicing = data.frame(genename = names(logBF_splicing), logBF_splicing = as.vector(logBF_splicing))
    
    gene_BF_table = data.frame(genename = names(logBF_noncoding), logBF_noncoding = as.vector(logBF_noncoding))
    gene_BF_table = merge(gene_BF_table, logBF_splicing, by = "genename")                           
    
    gene_BF_table = merge(gene_BF_table, coding_bayes_table, by = "genename")
    gene_BF_table$logBF_all = gene_BF_table$logBF_noncoding + gene_BF_table$logBF_splicing + gene_BF_table$logBF_coding
    gene_BF_table[,c("BF_noncoding","BF_splicing","BF_coding","BF_all")] = exp(gene_BF_table[,c("logBF_noncoding","logBF_splicing","logBF_coding","logBF_all")])
    gene_BF_table = gene_BF_table[order(gene_BF_table$BF_coding, decreasing = TRUE),]
    gene_BF_table$FDR_coding = Bayesian.FDR(gene_BF_table$BF_coding, TADA_p0)$FDR
    gene_BF_table = gene_BF_table[order(gene_BF_table$BF_all, decreasing = TRUE),]
    gene_BF_table$FDR_all = Bayesian.FDR(gene_BF_table$BF_all, TADA_p0)$FDR
  }
  system(paste("rm ", prefix, "_temp*", sep = ""))
  list(gene_BF_table = gene_BF_table, genes_with_no_epi = genes_with_no_epi)
}

# mut_file = "../analysis/160229_data_for_analysis.bed"
# window_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed"
# cg_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.fasta.cg"
# mutrate_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.mutrate"
# # test = coverage[sample(seq(1,nrow(coverage)),10000),]
# # 
# # test_model = glm(mut_count~cg+coding+promoter+nf, family = poisson(),data = test)
# # 
# # summary(out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(mutrate)), family = poisson, data = coverage))
#  gene_assign_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename.bed"
#  overlap = 1e-9
#  gene_prior_file = "../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_combined_qvalue.txt"
#  epigenomic_marks = "../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed"
# #mut_file = "./temp2.bed"
# #window_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed"
# #cg_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.fasta.cg"
# #mutrate_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate"
# #sample_size = 693
# # test = coverage[sample(seq(1,nrow(coverage)),10000),]
# 
# # test_model = glm(mut_count~cg+coding+promoter+nf, family = poisson(),data = test)
# 
# # summary(out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(mutrate)), family = poisson, data = coverage))
# 
# 
# #  