old = new.env()
load("/media/yuwen/F/ASD/integrating_all_data/150812_all_data_burden_analysis_v2.Rdata", old)
load("/media/yuwen/F/ASD/data/debug_region_list_073116_data_matrix.Rdata")
ASD_data_matrix = data_matrix[data_matrix$Prediction == "dnv_proband",]
> dim(ASD_data_matrix)
[1] 8069   44
control_data_matrix = data.frame(Chrom = old$mutation$chr,
                                 Start = factor(old$mutation$start),
                                 End = factor(old$mutation$end),
                                 Ref = old$all[,4],
                                 Alt = old$all[,5],
                                 Prediction = old$mutation$phenotype,
                                 ID = old$mutation$index,
                                 Type = old$mutation$mut_type,
                                 Sample = "NA",
                                 GERP = old$gerp[,8],
                                 Eigen = 0,
                                 CADD13_PHRED = old$CADD_score[,1],
                                 phyloP46wayAllElements = old$phylop[,6])
control_data_matrix = control_data_matrix[control_data_matrix$Prediction == "control",]
control_data_matrix$Prediction = "dnv_sibling"
control_data_matrix = data.frame(control_data_matrix, matrix(NA, nrow(control_data_matrix), 31))
colnames(control_data_matrix)[14:44] = colnames(ASD_data_matrix)[14:44]
data_matrix = rbind(ASD_data_matrix, control_data_matrix)
data_matrix$ID = seq(1,nrow(data_matrix))
data_matrix$Noonan_brain_motif_top50pct_expressed_motif.txt.ref_alt.3.m.pvalue = 0
rm(control_data_matrix, ASD_data_matrix)
rm(old)
