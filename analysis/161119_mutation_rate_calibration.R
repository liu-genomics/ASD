setwd("/media/yuwen/F/ASD/analysis/")
source("../lib/161117_glm_for_mutation_count.R")
# remeber only use SNV data

new_control = new.env()
load("../data/debug_region_list_073116_using_258_control_data_matrix.Rdata", new_control)
write.table(new_control$data_matrix[new_control$data_matrix$Prediction == "dnv_sibling" & new_control$data_matrix$Type == "SNV",c(1,2,3,7)], "temp.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
control_model_258 = run_glm_for_mutation("./temp.bed", "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed",
                                         "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.fasta.cg",
                                         "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate",
                                         258)

write.table(new_control$data_matrix[new_control$data_matrix$Prediction == "dnv_proband" & new_control$data_matrix$Type == "SNV",c(1,2,3,7)], "temp.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
ASD_model_manuscript = run_glm_for_mutation("./temp.bed", "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed",
                                         "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.fasta.cg",
                                         "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate",
                                         314)

real_new_control = new.env()
load("../data/debug_region_list_073116_data_matrix.Rdata", real_new_control)
write.table(real_new_control$data_matrix[real_new_control$data_matrix$Prediction == "dnv_sibling" & real_new_control$data_matrix$Type == "SNV",c(1,2,3,7)], "temp.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
control_model_693 = run_glm_for_mutation("./temp.bed", "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed",
                                         "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.fasta.cg",
                                         "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate",
                                         693)



scherer_case = new.env()
load("../data/161026_Stephen_Scherer_Gemomic_Medicine_ASD_de_novo_SNVs_germline_with_NC_control_data_matrix.Rdata", scherer_case)
write.table(scherer_case$data_matrix[scherer_case$data_matrix$Prediction == "dnv_proband" & scherer_case$data_matrix$Type == "SNV" ,c(1,2,3,7)], "temp.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
ASD_model_Scherer = run_glm_for_mutation("./temp.bed", "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed",
                                         "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.fasta.cg",
                                         "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed.mutrate",
                                         200)


rm(new_control)
rm(real_new_control)
rm(scherer_case)
save.image("161119_mutation_rate_calibration_new.Rdata")
quit(0)
