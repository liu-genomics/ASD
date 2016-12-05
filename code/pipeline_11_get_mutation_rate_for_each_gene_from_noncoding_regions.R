setwd("/media/yuwen/F/ASD/code/")
promoter = read.delim("../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.promoter_1kb_genename.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
promoter = data.frame(promoter, c(1:length(promoter[,1])))
write.table(promoter[,c(1:3,5)], "mutrate_temp.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
system("../lib/bigWigAverageOverBed ../other_annotation/mutation_rate/Daly_mutrate.bw mutrate_temp.bed mutrate_temp.bed.output")
promoter_mutrate = read.delim("mutrate_temp.bed.output", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
promoter_mutrate = promoter_mutrate[order(promoter_mutrate[,1], decreasing = FALSE),]
promoter = data.frame(promoter, promoter_mutrate)

gene_mutrate_promoter = promoter[,c(4,9)]
colnames(gene_mutrate_promoter) = c("genename","mutrate")

enhancer = read.delim("../other_annotation/epigenomic_annotation/Noonan_brain_roadmap_union.yanyu_pipeline_enhancers.10000.bp_within_TSS.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
enhancer = data.frame(enhancer, c(1:length(enhancer[,1])))
write.table(enhancer[,c(1:3,6)], "mutrate_temp.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
system("../lib/bigWigAverageOverBed ../other_annotation/mutation_rate/Daly_mutrate.bw mutrate_temp.bed mutrate_temp.bed.output")
enhancer_mutrate = read.delim("mutrate_temp.bed.output", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
enhancer_mutrate = enhancer_mutrate[order(enhancer_mutrate[,1], decreasing = FALSE),]

enhancer = data.frame(enhancer, enhancer_mutrate)
gene_mutrate_enhancer = enhancer[,c(4,10)]
colnames(gene_mutrate_enhancer) = c("genename","mutrate")


gene_mutrate_promoter_plus_10kb = rbind(gene_mutrate_promoter, gene_mutrate_enhancer)
temp = by(gene_mutrate_promoter_plus_10kb[,'mutrate'], gene_mutrate_promoter_plus_10kb[,'genename'], sum)

gene_mutrate_promoter_plus_10kb_per_gene = data.frame(genename = names(temp),mutrate = as.vector(temp))

write.table(gene_mutrate_promoter_plus_10kb_per_gene,"../other_annotation/mutation_rate/Noonan_brain_roadmap_union.promoter_plus_10kb_yanyu_pipeline.mutrate", col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)

