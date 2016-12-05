old = new.env()
scherer = new.env()
load("../data/161026_Stephen_Scherer_Gemomic_Medicine_ASD_de_novo_SNVs_germline_with_NC_control_data_matrix.Rdata", scherer)
load("../data/debug_region_list_073116_data_matrix.Rdata", old)
old_ASD = old$data_matrix[old$data_matrix$Prediction == "dnv_proband",]
old_data_identifier = paste(old_ASD$Chrom, old_ASD$Start, old_ASD$End, old_ASD$Ref, old_ASD$Alt,sep = "_")

scherer_ASD = scherer$data_matrix[scherer$data_matrix$Prediction == "dnv_proband",]
scherer_data_identifier = paste(scherer_ASD$Chrom, scherer_ASD$Start, scherer_ASD$End, scherer_ASD$Ref, scherer_ASD$Alt,sep = "_")

scherer_ASD_filtered = scherer_ASD[!is.element(scherer_data_identifier, old_data_identifier),]

scherer_control = scherer$data_matrix[scherer$data_matrix$Prediction == "dnv_sibling",]

data_matrix = rbind(scherer_ASD_filtered,scherer_control)
data_matrix$ID = c(1:length(data_matrix[,1]))

rm(old)
rm(scherer)
rm(old_ASD)
rm(old_data_identifier)
rm(scherer_control)
rm(scherer_ASD)
rm(scherer_ASD_filtered)
rm(scherer_data_identifier)

family_519 = new.env()
load("../data/0703_region_list_080216_data_matrix.Rdata", family_519)
motif_names = family_519$motif_names
region_names = family_519$region_names
score_names = family_519$score_names

rm(family_519)

save.image("../data/161026_Filtered_Stephen_Scherer_Gemomic_Medicine_ASD_de_novo_SNVs_germline_with_NC_control_data_matrix.Rdata")

q(save = "no")