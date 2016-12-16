setwd("../analysis/")
source("../lib/161117_glm_for_mutation_count.R")
# remember only use SNV data
# what is different from 161120_mutation_rate_calibration_1.R is that for coding regions, I removed all utrs and exons without a gene name (or
# the genename is n/a)

# also used Noonan brain roadmap H3K27ac as a dummy variable. 
# when predict the effect of coding mutation rate for ASD cases, only use synonymous mutations. 

# only use top 6% genes ranked by posterior probability of being risk genes from high to low. 
# simulation and do estimation of parameters for multiple times

new_control = new.env()
load("../data/debug_region_list_073116_using_258_control_data_matrix.Rdata", new_control)


write.table(new_control$data_matrix[new_control$data_matrix$Prediction == "dnv_proband" & new_control$data_matrix$Type == "SNV",c(1,2,3,7)], "temp.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(new_control$data_matrix[new_control$data_matrix$Prediction == "dnv_proband" ,c(1:5,8,7)], "annovar_input_temp.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

set.seed(100)
estimate = list()
hessian = list()
for(i in 1:10){
ASD_model_manuscript = simulate_noncoding_mutations("./temp.bed", "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed",
                                         "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.fasta.cg",
                                         "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.mutrate",
                                         314,
                                         "../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed",
                                         1e-9,
                                         "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename.bed",
                                         TRUE,
                                         "/media/yuwen/Elements/ANNOVAR/annovar/",
                                         "annovar_input_temp.bed",
                                         "../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.94.txt",
                                         5,
                                         3,
                                         0.06)

estimate[[i]] = ASD_model_manuscript$mle$par
hessian[[i]] = ASD_model_manuscript$mle$hessian
}
rm(new_control)
system("rm annovar_input_temp.bed")
system("rm temp.bed")


save.image("../analysis/161214_simulate_mutation_and_estimate_effect_top6_pct_with_Noonan_brain_roadmap_union_1bp_overlap_using_synony_only_for_ASD_coding.Rdata")
quit("no",0)
