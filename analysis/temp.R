```{r,cache=TRUE}
# load relavant mutation data
setwd("/media/yuwen/F/ASD/160229_new_control_integrating_all_data/systematic_study_burden_bug_fixed/")
source("160229_screen_burden.R")

```

```{r,echo=FALSE}
#use functions that use 2014_jaspar_dummy that have qvalue smaller than 0.1
setwd("/media/yuwen/F/ASD/160229_new_control_integrating_all_data/systematic_study_burden_bug_fixed/")
source("160310_screen_burden_functions_test_more_genelist_focused.R") #doesn't have the full functions yet, because the annotation hasn't been updated fully
```

```{r,echo=FALSE}
setwd("/media/yuwen/F/ASD/160229_new_control_integrating_all_data/systematic_study_burden_bug_fixed/")
table_colname = c("ASD","control","frquency_ratio","burden","odds_ratio","pvalue")
```



### noonan brain +roadmap covered enhancers (10,000bp)+ promoters or 5' utrs
```{r}
system("cat ./union_for_promter_utr5_covered_by_k27ac_for_each_gene_enhancer_pair_overlap_mutation.txt ./mutation_overlap_enhancers_distance_to_TSS/enhancers.bed.10000.bp_within_TSS_overlap_mutation.txt > temp.txt")
system("mv temp.txt noonan_plus_roadmap_brain_covered_promoter_utr5_10000_enhancers_overlap_mutation.txt")
system("rm temp.txt")
file = "./noonan_plus_roadmap_brain_covered_promoter_utr5_10000_enhancers_overlap_mutation.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
#write mutation files that I need to send to Iuliana
getwd()
show_mut_for_Table1(file)
```


### Roadmap nonbrain EXPANDED  (10,000bp)+ promoters or 5' utrs
```{r}
system("cat ./promoter_utr5_covered_by_other_tissues_not_in_Noonan_plus_roadmap_brain_with_previous_overlap_mutation.bed ./mutation_overlap_roadmap_nonbrain_k27ac_expanded_distance_to_TSS/enhancers.bed.10000.bp_within_TSS_overlap_mutation.txt > temp.txt")
file = "./temp.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
system("rm temp.txt")
```

### E081_E082_intersection (10,000bp)+ promoters or 5' utrs
```{r}
system("cat ./promoter_utr5_covered_by_E081_E082_intersection_overlap_mutation.txt ./mutation_overlap_E081_E082_intersection_distance_to_TSS/enhancers.bed.10000.bp_within_TSS_overlap_mutation.txt > temp.txt")
file = "./temp.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
system("rm temp.txt")
```

### non fetal brain Roadmap DNase(10,000bp)+ promoters or 5' utrs
```{r}
system("cat ./promoter_utr5_covered_by_nonbrain_DNAse_overlap_mutation.txt ./mutation_overlap_all_roadmap_dnase_non_intersect_E081_E082_distance_to_TSS/enhancers.bed.10000.bp_within_TSS_overlap_mutation.txt > temp.txt")
file = "./temp.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
system("rm temp.txt")
```

### noonan brain +roadmap covered enhancers (25,000bp)+ promoters or 5' utrs
```{r}
system("cat ./union_for_promter_utr5_covered_by_k27ac_for_each_gene_enhancer_pair_overlap_mutation.txt ./mutation_overlap_enhancers_distance_to_TSS/enhancers.bed.25000.bp_within_TSS_overlap_mutation.txt > temp.txt")
system("mv temp.txt noonan_plus_roadmap_brain_covered_promoter_utr5_25000_enhancers_overlap_mutation.txt")
system("rm temp.txt")
file = "./noonan_plus_roadmap_brain_covered_promoter_utr5_25000_enhancers_overlap_mutation.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
```

### Roadmap nonbrain EXPANDED (25,000bp)+ promoters or 5' utrs
```{r}
system("cat ./promoter_utr5_covered_by_other_tissues_not_in_Noonan_plus_roadmap_brain_with_previous_overlap_mutation.bed ./mutation_overlap_roadmap_nonbrain_k27ac_expanded_distance_to_TSS/enhancers.bed.25000.bp_within_TSS_overlap_mutation.txt > temp.txt")
file = "./temp.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
system("rm temp.txt")
```

### E081_E082_intersection (25,000bp)+ promoters or 5' utrs
```{r}
system("cat ./promoter_utr5_covered_by_E081_E082_intersection_overlap_mutation.txt ./mutation_overlap_E081_E082_intersection_distance_to_TSS/enhancers.bed.25000.bp_within_TSS_overlap_mutation.txt > temp.txt")
file = "./temp.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
system("rm temp.txt")
```

### non fetal brain Roadmap DNase(25,000bp)+ promoters or 5' utrs
```{r}
system("cat ./promoter_utr5_covered_by_nonbrain_DNAse_overlap_mutation.txt ./mutation_overlap_all_roadmap_dnase_non_intersect_E081_E082_distance_to_TSS/enhancers.bed.25000.bp_within_TSS_overlap_mutation.txt > temp.txt")
file = "./temp.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
system("rm temp.txt")
```

### Noonan brain specific 
```{r}
file = "./mutation_overlap_Noonan_brain_specific_distance_to_TSS/enhancers.bed.1000000.bp_within_TSS_overlap_mutation.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
```

### Noonan three other tissues
```{r}
file = "./mutation_overlap_Noonan_nonbrain_distance_to_TSS/enhancers.bed.1000000.bp_within_TSS_overlap_mutation.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
```

### promoters and utr5s covered by H3K27ac
```{r}
file = "./union_for_promter_utr5_covered_by_k27ac_for_each_gene_enhancer_pair_overlap_mutation.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
```

### within 10kb of TSS covered by H3K27ac, not including promoters and utr5s
```{r}
file = "./mutation_overlap_enhancers_distance_to_TSS/enhancers.bed.10000.bp_within_TSS_overlap_mutation.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
```

### within 10kb-25kb of TSS covered by H3K27ac, not including promoters and utr5s
```{r}
system("grep -F -x -v -f ./mutation_overlap_enhancers_distance_to_TSS/enhancers.bed.10000.bp_within_TSS_overlap_mutation.txt ./mutation_overlap_enhancers_distance_to_TSS/enhancers.bed.25000.bp_within_TSS_overlap_mutation.txt > temp.txt")
file = "./temp.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
```

### all promoters in cluding promoters and utr5s not covered by H3K27ac
```{r}
file = "./union_for_promter_utr5_for_each_gene_enhancer_pair_overlap_mutation.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
```


### Noonan 10kb enhancers flanking regions (removing promoter, utr5s, and exons, removing k27ac regions using 1Mb distance cutoff to assign genes)
```{r}
file = "./mutation_overlap_enhancers_10kb_flanking_distance_to_TSS/enhancers.bed.1000000.bp_within_TSS_overlap_mutation.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
```

### Noonan 10kb enhancers (combined with k27ac promoters and utr5)flanking regions  (removing promoter, utr5s, and exons, removing k27ac regions using 1Mb distance cutoff to assign genes)
```{r}
file = "./mutation_overlap_enhancers_10kb_plus_k27ac_promoter_utr5_flanking_distance_to_TSS/enhancers.bed.1000000.bp_within_TSS_overlap_mutation.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
```

### Roadmap nonbrain k27ac batch1 10kb, not including promoters and 5'utr

```{r}
file = "./mutation_overlap_nonBrain_k27ac_batch1_distance_to_TSS/enhancers.bed.10000.bp_within_TSS_overlap_mutation.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
```

### Roadmap nonbrain k27ac batch2 10kb, not including promoters and 5'utr

```{r}
file = "./mutation_overlap_nonBrain_k27ac_batch2_distance_to_TSS/enhancers.bed.10000.bp_within_TSS_overlap_mutation.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
```

### Roadmap nonbrain k27ac batch1, only including promoters and 5'utr

```{r}
file = "./union_for_promoter_utr5_covered_by_nonBrain_k27ac_batch1_for_each_gene_enhancer_pair_overlap_mutation.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
```

### Roadmap nonbrain k27ac batch2, only including promoters and 5'utr

```{r}
file = "./union_for_promoter_utr5_covered_by_nonBrain_k27ac_batch2_for_each_gene_enhancer_pair_overlap_mutation.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
```

### Roadmap nonbrain k27ac batch1 10kb, including enhancers and promoters and 5'utr

```{r}
system("cat ./union_for_promoter_utr5_covered_by_nonBrain_k27ac_batch1_for_each_gene_enhancer_pair_overlap_mutation.txt ./mutation_overlap_nonBrain_k27ac_batch1_distance_to_TSS/enhancers.bed.10000.bp_within_TSS_overlap_mutation.txt > temp.txt")
file = "./temp.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
system("rm temp.txt")
```

### Roadmap nonbrain k27ac batch2 10kb, including enhancers and promoters and 5'utr

```{r}
system("cat ./union_for_promoter_utr5_covered_by_nonBrain_k27ac_batch2_for_each_gene_enhancer_pair_overlap_mutation.txt ./mutation_overlap_nonBrain_k27ac_batch2_distance_to_TSS/enhancers.bed.10000.bp_within_TSS_overlap_mutation.txt > temp.txt")
file = "./temp.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
system("rm temp.txt")
```

### combining batch1 and batch2 10kb, including enhancers and promoters and 5'utr
```{r}
system("cat ./union_for_promoter_utr5_covered_by_nonBrain_k27ac_batch1_plus_batch2_for_each_gene_enhancer_pair_overlap_mutation.txt ./mutation_overlap_nonBrain_k27ac_batch1_plus_batch2_distance_to_TSS/enhancers.bed.10000.bp_within_TSS_overlap_mutation.txt > temp.txt")
file = "./temp.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)
system("rm temp.txt")
```

### combining batch1 and batch2 10kb, including enhancers and promoter and 5'utr, but remove any mutations that are also in brain enhancers
```{r}
system("cat ./union_for_promoter_utr5_covered_by_nonBrain_k27ac_batch1_plus_batch2_for_each_gene_enhancer_pair_overlap_mutation.txt ./mutation_overlap_nonBrain_k27ac_batch1_plus_batch2_distance_to_TSS/enhancers.bed.10000.bp_within_TSS_overlap_mutation.txt > temp.txt")
mut1 = read.delim("./temp.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
mut2 = read.delim("./noonan_plus_roadmap_brain_covered_promoter_utr5_10000_enhancers_overlap_mutation.txt", head = FALSE, sep = "\t", stringsAsFactors = FALSE)

fun.12 <- function(x.1,x.2,...){
  x.1p <- do.call("paste", x.1)
  x.2p <- do.call("paste", x.2)
  x.1[! x.1p %in% x.2p, ]
}

mut = fun.12(mut1,mut2)
knitr::kable(show_burden_for_selected_features(mut, input_type = "table"),digit = c(rep(2,5),5), col.names = table_colname)
```

### brain enhancers 10kb, including enhancers and promoter and 5'utr, but remove any mutations that are also in combining batch1 and batch2 10kb, 
```{r}
system("cat ./union_for_promoter_utr5_covered_by_nonBrain_k27ac_batch1_plus_batch2_for_each_gene_enhancer_pair_overlap_mutation.txt ./mutation_overlap_nonBrain_k27ac_batch1_plus_batch2_distance_to_TSS/enhancers.bed.10000.bp_within_TSS_overlap_mutation.txt > temp.txt")
mut2 = read.delim("./temp.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
mut1 = read.delim("./noonan_plus_roadmap_brain_covered_promoter_utr5_10000_enhancers_overlap_mutation.txt", head = FALSE, sep = "\t", stringsAsFactors = FALSE)

fun.12 <- function(x.1,x.2,...){
  x.1p <- do.call("paste", x.1)
  x.2p <- do.call("paste", x.2)
  x.1[! x.1p %in% x.2p, ]
}

mut = fun.12(mut1,mut2)
knitr::kable(show_burden_for_selected_features(mut, input_type = "table"),digit = c(rep(2,5),5), col.names = table_colname)
```


### batch 2 nonbrain k27ac after removing the three immune cell types, not including promoters and 5'utr
```{r}
file = "./mutation_overlap_nonBrain_k27ac_batch2_no_immune_distance_to_TSS/enhancers.bed.10000.bp_within_TSS_overlap_mutation.txt"
knitr::kable(show_burden_for_selected_features(file),digit = c(rep(2,5),5), col.names = table_colname)

```