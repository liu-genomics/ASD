# function to run glm regression for window based mutation count
run_glm_for_mutation <- function(mut_file, window_file, cg_file, mutrate_file, sample_size, epigenomic_marks = "no"){
  # [mut_file] is the mutation location file from one dataset. e.g.,../analysis/160229_data_for_analysis.bed, chr \t start \t end \t index,  1-based start and end
  # [window_file] is the file with genome partitioned into windows. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed
  # [cg_file] is the file that has CG-content calculated, matching [window_file]. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.cg
  # [mutrate_file] is the file that has mutrate data. e.g., "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate"
  # [sample_size] is the number of individuals
  # [epitenomic_marks] if it is "no", then glm won't use if a mutation is in h3k27ac as a categorical variable, if it is a epigenomic annotation name.
  # for example, "../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed". It is fine to directly use original epigenomics bed files
  # instead of the partitioned files. Because the window_file already has annotations for each window in terms of whetehr they are promoter or coding. 
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
    command = paste("bedtools intersect -f 0.5 -a ", window_file , " -b ", epigenomic_marks, " -wa > ", paste(prefix,"_temp_overlap_epi.bed", sep = ""),sep = "")
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




# mut_file = "../analysis/160229_data_for_analysis.bed"
# window_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed"
# cg_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.fasta.cg"
# mutrate_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate"
# test = coverage[sample(seq(1,nrow(coverage)),10000),]
# 
# test_model = glm(mut_count~cg+coding+promoter+nf, family = poisson(),data = test)
# 
# summary(out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(mutrate)), family = poisson, data = coverage))


#mut_file = "../analysis/160229_data_for_analysis.bed"
#window_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed"
#cg_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.fasta.cg"
#mutrate_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate"
#sample_size = 693
#epigenomic_marks = "../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed"
# test = coverage[sample(seq(1,nrow(coverage)),10000),]

# test_model = glm(mut_count~cg+coding+promoter+nf, family = poisson(),data = test)

# summary(out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(mutrate)), family = poisson, data = coverage))


#  