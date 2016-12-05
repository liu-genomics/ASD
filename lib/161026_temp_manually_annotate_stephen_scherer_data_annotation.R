ASD_mutation = read.delim("/media/yuwen/F/ASD/data/161026_Stephen_Scherer_Gemomic_Medicine_ASD_de_novo_SNVs_germline.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
ASD_mutation <- data.frame(ASD_mutation, seq(1,length(ASD_mutation[,1])))
ASD_output_bed <- ASD_mutation
ASD_output_bed[,2] <- ASD_output_bed[,2]-1
write.table(ASD_output_bed, "/media/yuwen/Elements/genomic_annotation/161026_Stephen_Scherer_Genomic_Medicine_ASD_de_novo_SNVs_germline.bed", col.names = FALSE, row.names = FALSE,
            sep = "\t",quote = FALSE)
gerp <- read.delim("/media/yuwen/Elements/genomic_annotation/161026_Stephen_Scherer_Genomic_Medicine_ASD_de_novo_SNVs_germline_gerp.tab", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
phylop <- read.delim("/media/yuwen/Elements/genomic_annotation/161026_Stephen_Scherer_Genomic_Medicine_ASD_de_novo_SNVs_germline_phylop_100way.tab", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

gerp <- gerp[order(gerp[,1]),]
phylop <- phylop[order(phylop[,1]),]

cadd_raw <- read.delim("/media/yuwen/Elements/genomic_annotation/161026_Stephen_Scherer_Genomic_Medicine_ASD_de_novo_SNVs_germline_raw_CADD.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

cadd <- merge(cadd_raw, ASD_mutation, by.x = "V2", by.y = "V2")

cadd <- cadd[order(cadd[,9]),]


mutation <- new.env()
load("/media/yuwen/F/ASD/data/debug_region_list_073116_data_matrix.Rdata", envir = mutation)

data_matrix <- mutation$data_matrix
data_matrix <- data_matrix[1:9774,]
data_matrix[,1:3] = ASD_mutation[,1:3]
data_matrix[,4:5] = cadd[,3:4]
data_matrix[,6] = "dnv_proband"
data_matrix$GERP = gerp[,6]
data_matrix$Eigen = 0
data_matrix$CADD13_PHRED = cadd[,6]
data_matrix$phyloP46wayAllElements = phylop[,6]

control <- mutation$data_matrix[mutation$data_matrix$Prediction == "dnv_sibling",]
data_matrix <- rbind(data_matrix, control)

data_matrix$ID <- seq(1,length(data_matrix$ID))

motif_names <- mutation$motif_names
regions_names <- mutation$region_names
score_names <- mutation$score_names
rm(mutation)
data_matrix$Eigen = 0
data_matrix[1:9774,]$Type = "SNV"
save.image("/media/yuwen/F/ASD/data/161026_Stephen_Scherer_Gemomic_Medicine_ASD_de_novo_SNVs_germline_with_NC_control_data_matrix.Rdata")
