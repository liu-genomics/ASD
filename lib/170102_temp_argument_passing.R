new_control = new.env()
load("../data/debug_region_list_073116_using_258_control_data_matrix.Rdata", new_control)
write.table(new_control$data_matrix[new_control$data_matrix$Prediction == "dnv_proband" & new_control$data_matrix$Type == "SNV",c(1,2,3,7)], "temp.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(new_control$data_matrix[new_control$data_matrix$Prediction == "dnv_proband" ,c(1:5,8,7)], "annovar_input_temp.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

mut_file = "./temp.bed"
window_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed"
cg_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.fasta.cg"
mutrate_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.mutrate"
# test = coverage[sample(seq(1,nrow(coverage)),10000),]
# 
# test_model = glm(mut_count~cg+coding+promoter+nf, family = poisson(),data = test)
# 
# summary(out.offset <- glm(mut_count ~ coding+promoter+cg+offset(log(mutrate)), family = poisson, data = coverage))
gene_assign_file = "../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed_genename.bed"
overlap = 1e-9
gene_prior_file = "../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.94.txt"
epigenomic_marks_list = c("../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.bed","../other_annotation/epigenomic_annotation/Epigenome_E081_E082_intersection.bed","../other_annotation/epigenomic_annotation/encode_DHS_union.bed")
rm_nonsyn = TRUE
annovar_input = "annovar_input_temp.bed"
annovar_folder = "/media/yuwen/Elements/ANNOVAR/annovar/"
sample_size = 314
effect_size = c(5,3,2)
sequence_annotation_list = c("../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons_no_utr_no_na_genes.50bp_window.bed.phastcons100way")
gene_prior_file = "../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.94.txt"
set.seed(100)
