setwd("../analysis/")
source("../lib/161117_glm_for_mutation_count.R")
# remember only use SNV data
# what is different from 161120_mutation_rate_calibration_1.R is that for coding regions, I removed all utrs and exons without a gene name (or
# the genename is n/a)

# also used Noonan brain roadmap H3K27ac as a dummy variable. 
# when predict the effect of coding mutation rate for ASD cases, only use synonymous mutations. 

# Three epigenomic annotation (Noonan and Roadmap Brain H3k27ac, ENCODE DHS, and Roadmap DHS) and one conservation score(phastcons100way) were used to correct mutation rate and to estimate relative risks. 


# only use top 10% TADA genes for estimating relative risks
new_control = new.env()
load("../data/debug_region_list_073116_using_258_control_data_matrix.Rdata", new_control)


write.table(new_control$data_matrix[new_control$data_matrix$Prediction == "dnv_proband" & new_control$data_matrix$Type == "SNV",c(1,2,3,7)], "temp.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(new_control$data_matrix[new_control$data_matrix$Prediction == "dnv_proband" ,c(1:5,8,7)], "annovar_input_temp.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

ASD_model_manuscript = verified_effect_size_estimate_noncoding_mutations_v2(mut_file = "./temp.bed", 
                                                                            window_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed",
                                                                            cg_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.fasta.cg",
                                                                            mutrate_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.mutrate",
                                                                            sample_size = 314,
                                                                            epigenomic_marks_list = c("../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed","../other_annotation/epigenomic_annotation/Epigenome_E081_E082_intersection.bed","../other_annotation/epigenomic_annotation/encode_DHS_union.bed"),
                                                                            overlap = 1e-9,
                                                                            gene_assign_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename.bed",
                                                                            rm_nonsyn = TRUE,
                                                                            annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/",
                                                                            annovar_input = "annovar_input_temp.bed",
                                                                            gene_prior_file = "../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.94.txt",
                                                                            sequence_annotation_list = c("../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.phastcons100way"),
                                                                            optimization_prop = 1)



rm(new_control)
system("rm annovar_input_temp.bed")
system("rm temp.bed")


save.image("../analysis/161229_estimate_effect_with_three_epi_and_one_conservation_1bp_overlap_all_tada_using_synony_only_for_ASD_coding.Rdata")
quit("no",0)
