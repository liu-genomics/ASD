setwd("/media/yuwen/F/ASD/analysis/")
source("../lib/161117_glm_for_mutation_count.R")
# remember only use SNV data
# what is different from 161120_mutation_rate_calibration_1.R is that for coding regions, I removed all utrs and exons without a gene name (or
# the genename is n/a)

# also used Noonan brain roadmap H3K27ac as a dummy variable. 

new_control = new.env()
load("../data/debug_region_list_073116_using_258_control_data_matrix.Rdata", new_control)
write.table(new_control$data_matrix[new_control$data_matrix$Prediction == "dnv_sibling" & new_control$data_matrix$Type == "SNV",c(1,2,3,7)], "temp.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
control_model_258 = run_glm_for_mutation("./temp.bed", "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed",
                                         "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.fasta.cg",
                                         "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.mutrate",
                                         258,
                                         "../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed")
control_model_258_summary = summary(control_model_258)$coefficients
rm(control_model_258)

write.table(new_control$data_matrix[new_control$data_matrix$Prediction == "dnv_proband" & new_control$data_matrix$Type == "SNV",c(1,2,3,7)], "temp.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
ASD_model_manuscript = run_glm_for_mutation("./temp.bed", "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed",
                                         "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.fasta.cg",
                                         "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.mutrate",
                                         314,
                                         "../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed")
ASD_model_manuscript_summary = summary(ASD_model_manuscript)$coefficients
rm(ASD_model_manuscript)
rm(new_control)

real_new_control = new.env()
load("../data/debug_region_list_073116_data_matrix.Rdata", real_new_control)
write.table(real_new_control$data_matrix[real_new_control$data_matrix$Prediction == "dnv_sibling" & real_new_control$data_matrix$Type == "SNV",c(1,2,3,7)], "temp.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
control_model_693 = run_glm_for_mutation("./temp.bed", "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed",
                                         "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.fasta.cg",
                                         "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.mutrate",
                                         693,
                                         "../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed")

control_model_693_summary = summary(control_model_693)$coefficients
rm(control_model_693)
rm(real_new_control)

scherer_case = new.env()
load("../data/161026_Stephen_Scherer_Gemomic_Medicine_ASD_de_novo_SNVs_germline_with_NC_control_data_matrix.Rdata", scherer_case)
write.table(scherer_case$data_matrix[scherer_case$data_matrix$Prediction == "dnv_proband" & scherer_case$data_matrix$Type == "SNV" ,c(1,2,3,7)], "temp.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
ASD_model_Scherer = run_glm_for_mutation("./temp.bed", "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed",
                                         "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.fasta.cg",
                                         "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.mutrate",
                                         200,
                                         "../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed")
ASD_model_Scherer_summary = summary(ASD_model_Scherer)$coefficients
rm(ASD_model_Scherer)
rm(scherer_case)


save.image("../analysis/161119_mutation_rate_calibration_with_Noonan_brain_roadmap_union.Rdata")
quit("no",0)
