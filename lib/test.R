## Use jaspar 2014 motif from both ref alleles and mutated allels, with qvalue smaller than 0.1 #####################
## has a lot of comparison with certain genelists that have been associated with autism or brain in some way.
## pvalue is computed for each group using the baseline.
#add rownames to the last column to help visualization. 

coding_tada = read.delim("../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
coding_tada_denovo = data.frame(coding_tada$RefSeqName, coding_tada$TadaName, coding_tada$mut.rate, coding_tada$dn.LoF, coding_tada$dn.mis3, coding_tada$BF.SNV.dn, coding_tada$BF.combined, coding_tada$qvalue.combined)
colnames(coding_tada_denovo) = c("RefSeqName", "TadaName", "mut.rate", "dn.LoF", "dn.mis3", "BF.SNV.dn", "BF.combined", "qvalue.combined")

brainspan = new.env()
load("../other_annotation/brainspan/160206_brainspan_expression.Rdata", envir = brainspan)


##########################Function to calculate burden ############################################
burden <- function (data){ # calculate SNV burden
  #data is a two columns data matrix, where 1st column is mutation ID, 2nd column is assigned gene
  a = intersect(data[,1], mutation$ASD_effective_SNV_ID)
  b = intersect(data[,1], mutation$control_effective_SNV_ID)
  c(length(a), length(b), (length(a)/length(mutation$ASD_effective_SNV_ID))/(length(b)/length(mutation$control_effective_SNV_ID))) 
}

################ Function to calculate burden given a mutation overlap feature file and a exp cutoff ##################
burden_by_brainspan <- function(data, lower_bound, upper_bound){
  genes_in_range = brainspan$gene_with_mutrate[brainspan$gene_with_mutrate$exp_mean >= lower_bound 
                                               & brainspan$gene_with_mutrate$exp_mean <= upper_bound
                                               & !is.na(gene_with_mutrate$exp_mean),]$gene
  data = data[is.element(data[,2],genes_in_range),]
  list(burden = burden(data), mut = data) #also return mutations after filtering. 
}

##############function to calculate burden given a mutatio overlap motif #######################################
burden_with_motif <- function(data, motif_dummy){#motif dummy variable 1/0 for the original mutation file (mutation$mutation)
  motif_index = mutation$mutation[motif_dummy == 1,]$index
  data = data[is.element(data[,1], motif_index),]
  list(burden = burden(data), mut = data) #also return mutations after filtering.
}
 ##############function to calculate burden given a motif_gain or motif break #######################################
burden_with_motif_alt <- function(data, motif_alt_index){
  data = data[is.element(data[,1], motif_alt_index),]
  list(burden = burden(data), mut = data) #also return mutations after filtering.
}

##############function to calculate burden given a Fitcons score cutoff #######################################
burden_with_fitcons <- function(data, fitcons, lowerbound, upperbound){
  fitcons_index = fitcons[fitcons[,6] <= upperbound & fitcons[,6] >= lowerbound, 1]
  data = data[is.element(data[,1], fitcons_index),]
  list(burden = burden(data), mut = data) #also return mutations after filtering.
}

#############function calculate burden given a mutation eigen score cutoff ##############################
burden_with_eigen <-function(data, mut_eigen, lowerbound){
  eigen_index = mutation$mutation[mut_eigen >= lowerbound & !is.na(mut_eigen),]$index
  data = data[is.element(data[,1], eigen_index),]
  list(burden = burden(data), mut = data) #also return mutations after filtering.
}



##############function to calculate burden given a mutation overlap motif or fitcons > 0.125#######################################
burden_with_fitcons_or_motif <- function(data, fitcons,lowerbound, upperbound, motif_dummy){
  fitcons_index = fitcons[fitcons[,6] <= upperbound & fitcons[,6] >= lowerbound, 1]
  motif_index = mutation$mutation[motif_dummy == 1,]$index
  data = data[is.element(data[,1], c(fitcons_index, motif_index)),]
  list(burden = burden(data), mut = data) #also return mutations after filtering.
}

##############function to calculate burden given a mutation overlap motif or CADD_score > certain number  #######################################
burden_with_CADD_or_motif <- function(data, CADD_score,lowerbound, upperbound, motif_dummy){
  CADD_index = CADD_score[CADD_score[,2] >=  lowerbound & CADD_score[,2] <= upperbound & !is.na(CADD_score[,2]),1]
  motif_index = mutation$mutation[motif_dummy == 1,]$index
  data = data[is.element(data[,1], c(CADD_index, motif_index)),]
  list(burden = burden(data), mut = data) #also return mutations after filtering.
}

##############function to calculate burden given a PhyloP score cutoff #######################################
burden_with_phylop <- function(data, phylop, lowerbound, upperbound){
  phylop_index = phylop[phylop[,6] <= upperbound & phylop[,6] >= lowerbound, 1]
  data = data[is.element(data[,1], phylop_index),]
  list(burden = burden(data), mut = data) #also return mutations after filtering.
}

##############function to calculate burden given a CADD score cutoff #######################################
CADD_score_with_index = data.frame(seq(1:length(mutation$mutation[,1])),mutation$CADD_score[,1])
burden_with_CADD <- function(data, CADD_score, lowerbound, upperbound){#CADD score is a data.frame with only one column
  CADD_index = CADD_score[CADD_score[,2] >=  lowerbound & CADD_score[,2] <= upperbound & !is.na(CADD_score[,2]),1]
  data = data[is.element(data[,1], CADD_index),]
  list(burden = burden(data), mut = data) #also return mutations after filtering.
}

##############function to calculate burden given a GERP score cutoff #######################################
burden_with_gerp <- function(data, gerp_score, lowerbound, upperbound){#gerp score is a data.frame with 8 columns
  gerp_index = gerp_score[gerp_score[,8] >=  lowerbound & gerp_score[,8] <= upperbound & !is.na(gerp_score[,8]),4]
  data = data[is.element(data[,1], gerp_index),]
  list(burden = burden(data), mut = data) #also return mutations after filtering.
}

##############function to calculate burden for mutations in enhancers that are near a specific set of genes #######################################
burden_with_geneset <- function(data, geneset){
  data = data[is.element(data[,2], geneset),]
  list(burden = burden(data), mut = data) #also return mutations after filtering.
}
############ function to calculate burden for mutations that might affect splicing ###########################################
# this is the code to interagate not all possible spice mutations, but mutations that are covered by k27ac within some distance cutoff or in coding regions
burden_with_spidex <-function(data, spidex, lowerbournd,upperbound){
  index = spidex[spidex$dpsi_zscore >= lowerbound & spidex$dpsi_zscore<= upperbound & !is.na(spidex$dpsi_zscore),]$index
  data = data[is.element(data[,1],index),]
  list(burden = burden(data), mut = data) #also return mutations after filtering
}
## the function below is to just calculate the burden of all mutations that have a splicing score greater than some value
burden_with_spidex2 <- function(spidex, lowerbound, upperbound){
  data = spidex[spidex$dpsi_zscore >= lowerbound & spidex$dpsi_zscore<= upperbound & !is.na(spidex$dpsi_zscore),]
  a = intersect(data$index, mutation$ASD_effective_SNV_ID)
  b = intersect(data$index, mutation$control_effective_SNV_ID)
  c(length(a), length(b), (length(a)/length(mutation$ASD_effective_SNV_ID))/(length(b)/length(mutation$control_effective_SNV_ID))) 
}

## the function then calculate burden of splicing mutations with regard to gene sets
burden_with_spidex2_geneset <- function(spidex, lowerbound, upperbound, geneset){
  data = spidex[is.element(spidex$gene, geneset) & !is.na(spidex$gene),]
  burden_with_spidex2(data, lowerbound, upperbound)
}



## the get genes with splicing mutations of splicing mutations with regard to gene sets
get_mutations_with_spidex2_geneset <- function(spidex, lowerbound, upperbound, geneset){
  data = spidex[is.element(spidex$gene, geneset) & !is.na(spidex$gene) & spidex$dpsi_zscore >= lowerbound & spidex$dpsi_zscore <= upperbound,]
  data = data[is.element(data[,1], mutation$ASD_effective_SNV_ID),]
  data
}

##############function to calculate burden for mutations in genes based on haploinsufficiency score #######################################
burden_with_hpi <- function(data, geneset,lowerbound){
  geneset = geneset[geneset[,2] >= lowerbound,1]
  data = data[is.element(data[,2], geneset),]
  list(burden = burden(data), mut = data) #also return mutations after filtering.
}

############## function to calculate burden for mutations in genes Nick trained ################################
burden_with_nick_gene <- function(data, geneset, minp_cutoff = 0, meanp_cutoff = 0){ 
  geneset = geneset[geneset$minp > minp_cutoff & geneset$meanp > meanp_cutoff,]$gene
  data = data[is.element(data[,2], geneset),]
  list(burden = burden(data), mut = data) #also return mutations after filtering.
}


############# function to calculate burden for coding regions ###########################################################
show_burden_for_coding <- function(data, exon_mut_type){ # a vector indicator the type of noncoding mutations
  #has the same length and mutation index ordering as the main mutation table
  #data is a mutation/gene overlapping table, where the mutations are in coding regions of genes. 
  syn = data[is.element(data[,1], mutation$mutation[exon_mut_type == "synonymous",]$index),]
  nonsyn = data[is.element(data[,1], mutation$mutation[exon_mut_type == "nonsynonymous",]$index),]
  nonsyn_plus_lof = data[is.element(data[,1], mutation$mutation[exon_mut_type != "synonymous" & 
                      exon_mut_type != "unknown" ,]$index),] #noncoding is annotated as unknown
  a = burden(syn)
  b = burden(nonsyn)
  c = burden(nonsyn_plus_lof)
  d = burden_with_geneset(syn, Jinyu_psy_gene[,1])$burden
  e = burden_with_geneset(nonsyn, Jinyu_psy_gene[,1])$burden
  f = burden_with_geneset(nonsyn_plus_lof, Jinyu_psy_gene[,1])$burden
  g = burden_with_geneset(syn,Brain_GO_gene[,1])$burden
  h = burden_with_geneset(nonsyn,Brain_GO_gene[,1])$burden
  i = burden_with_geneset(nonsyn_plus_lof,Brain_GO_gene[,1])$burden
  r = burden_with_geneset(nonsyn, coding_tada[coding_tada$qvalue.combined<0.5,2])$burden
  s = burden_with_geneset(nonsyn, constraint_gene[,1])$burden
  t = burden_with_geneset(nonsyn, Yuen_gene[,1])$burden
  u = burden_with_geneset(burden_with_motif(nonsyn, mutation$jaspar_2014_motif_q0.1_dummy)$mut,coding_tada[coding_tada$qvalue.combined<0.5,2])$burden
  v = burden_with_geneset(burden_with_CADD(nonsyn, CADD_score_with_index, 13.71,53)$mut,coding_tada[coding_tada$qvalue.combined<0.5,2])$burden
  w = burden_with_geneset(burden_with_CADD_or_motif(nonsyn, CADD_score_with_index, 13.71, 53,mutation$jaspar_2014_motif_q0.1_dummy)$mut, coding_tada[coding_tada$qvalue.combined<0.5,2])$burden
  x = burden_with_geneset(burden_with_CADD_or_motif(nonsyn, CADD_score_with_index, 13.71, 53,mutation$jaspar_2014_motif_q0.1_dummy)$mut, constraint_gene[,1])$burden
  y = burden_with_geneset(syn, coding_tada[coding_tada$qvalue.combined<0.5,2])$burden
  z = burden_with_geneset(syn, constraint_gene[,1])$burden
  aa = burden_with_geneset(syn, Yuen_gene[,1])$burden
  ab = burden_with_geneset(burden_with_motif(syn, mutation$jaspar_2014_motif_q0.1_dummy)$mut,coding_tada[coding_tada$qvalue.combined<0.5,2])$burden
  ac = burden_with_geneset(burden_with_CADD(syn, CADD_score_with_index, 13.71,53)$mut,coding_tada[coding_tada$qvalue.combined<0.5,2])$burden
  ad = burden_with_geneset(burden_with_CADD_or_motif(syn, CADD_score_with_index, 13.71, 53,mutation$jaspar_2014_motif_q0.1_dummy)$mut, coding_tada[coding_tada$qvalue.combined<0.5,2])$burden
  ae = burden_with_geneset(burden_with_CADD_or_motif(syn, CADD_score_with_index, 13.71, 53,mutation$jaspar_2014_motif_q0.1_dummy)$mut, constraint_gene[,1])$burden
  af = burden_with_geneset(nonsyn, Purcell_NG_2014_composite_gene[,1])$burden
  ag = burden_with_geneset(syn, Purcell_NG_2014_composite_gene[,1])$burden
  ah = burden_with_geneset(nonsyn, Olfactory_gene[,1])$burden
  ai = burden_with_geneset(nonsyn, DDG2P[,1])$burden
  aj = burden_with_geneset(nonsyn, DDG2P_confirmed_DD[,1])$burden
  output = data.frame(rbind(a,b,c,d,e,f,g,h,i,r,s,t,u,v,w,x,y,z,aa,ab,ac,ad,ae,af,ag,ah,ai,aj))
  rownames(output) = c("syn", "nonsyn", "nonsyn+lof","syn in Jinyu", "nonsyn in Jinyu", "nonsyn+lof in Jinyu",
                       "syn in brain GO", "nonsyn in brain GO", "nonsyn+lof in brain GO", "nonsyn+coding_TADA_0.5", "nonsyn+constraint_top5%", "nonsyn+yuen_NG_2014_gene", 
                       "nonsyn+motif_TADA_q0.5", "nonsyn+CADD_TADA_q0.5","+nonysn+motif_or_CADD_TADA_q0.5", "nonsyn+motif_or_CADD_constraint_top5%",
                       "syn+coding_TADA_0.5", "syn+constraint_top5%", "syn+yuen_NG_2014_gene", 
                       "syn+motif_TADA_q0.5", "syn+CADD_TADA_q0.5","+ysn+motif_or_CADD_TADA_q0.5", "syn+motif_or_CADD_constraint_top5%",
                       "nonsyn+Purcell", "syn+Purcell","nonsyn+olfactory_genes",
                       "nonsyn+DDG2P", "nonsyn+DDG2P_confirmed_DD")
  colnames(output) = c("ASD", "control", "burden")
  output
}


show_burden_for_coding_selected_features <- function(data, exon_mut_type){ # a vector indicator the type of noncoding mutations
  #has the same length and mutation index ordering as the main mutation table
  #data is a mutation/gene overlapping table, where the mutations are in coding regions of genes. 
  syn = data[is.element(data[,1], mutation$mutation[exon_mut_type == "synonymous",]$index),]
  nonsyn = data[is.element(data[,1], mutation$mutation[exon_mut_type == "nonsynonymous",]$index),]
  nonsyn_plus_lof = data[is.element(data[,1], mutation$mutation[exon_mut_type != "synonymous" & 
                                                                  exon_mut_type != "unknown" ,]$index),] #noncoding is annotated as unknown
  synonymous = burden(syn)
  nonsynonymous = burden(nonsyn)
  Tada_q0.1 = burden_with_geneset(nonsyn, coding_tada[coding_tada$qvalue.combined<0.1,2])$burden
  Tada_q0.2 = burden_with_geneset(nonsyn, coding_tada[coding_tada$qvalue.combined<0.2,2])$burden
  Tada_q0.3 = burden_with_geneset(nonsyn, coding_tada[coding_tada$qvalue.combined<0.3,2])$burden
  Tada_q0.4 = burden_with_geneset(nonsyn, coding_tada[coding_tada$qvalue.combined<0.4,2])$burden
  Tada_q0.5 = burden_with_geneset(nonsyn, coding_tada[coding_tada$qvalue.combined<0.5,2])$burden
  ASD_stringent = burden_with_geneset(nonsyn, stringent_ASD)$burden
  ASD_relaxed = burden_with_geneset(nonsyn, relaxed_ASD)$burden
  constraint_union_data1 = burden_with_geneset(nonsyn, constraint_union)$burden
  lowly_expressed = burden_with_geneset(nonsyn, c(brainspan$gene_qt1$gene, brainspan$gene_qt2$gene))$burden
  highly_expressed = burden_with_geneset(nonsyn, c(brainspan$gene_qt4$gene, brainspan$gene_qt3$gene))$burden
  olfactory_burden = burden_with_geneset(nonsyn, Olfactory_gene[,1])$burden
  output = data.frame(rbind(synonymous,
                            nonsynonymous,
                            Tada_q0.1,
                            Tada_q0.2,
                            Tada_q0.3,
                            Tada_q0.4,
                            Tada_q0.5,
                            ASD_stringent,
                            ASD_relaxed,
                            constraint_union_data1,
                            lowly_expressed,
                            highly_expressed,
                            olfactory_burden))
  rownames(output) = c("syn",
                       "nonsyn",
                       "Tada_q0.1",
                       "Tada_q0.2",
                       "Tada_q0.3",
                       "Tada_q0.4",
                       "Tada_q0.5",
                       "ASD stringent",
                       "ASD relaxed",
                       "constraint/haploinsuffient genes",
                       "lowly_expressed",
                       "highly_expressed",
                       "olfactory_burden")
  colnames(output) = c("ASD", "control", "burden")
  relative_burden = output$burden/output[1,]$burden
  get_fisher<-function(a,b){fisher.test(matrix(c(output[1,1],a,output[1,2],b),2,2),alternative = "less")$p.value}
  fisher.p = mapply(get_fisher,output$ASD, output$control)
  output = data.frame(output, relative_burden = relative_burden, fisher.p = fisher.p)
  output
}


# ############ original function of show_burden, with full functions ####################################################################
# show_burden <-function(mut_file_name, input_type = "file", mut_type = "noncoding"){
#   #input_type could be file or a table matching regions with mutations. 
#   # mut_type could be "noncoding" or "coding
#   if(input_type == "file"){
#     mut = read.delim(mut_file_name, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#   }
#   else if (input_type == "table"){
#     mut = mut_file_name
#   }
#   if(mut_type == "noncoding"){
#     a = burden(mut)
#     b = burden_with_geneset(mut, Jinyu_psy_gene[,1])$burden
#     c = burden_with_motif(mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
#     d = burden_with_fitcons(mut, mutation$fitcons, 0.125, max(mutation$fitcons[,6]))$burden
#     e = burden_with_fitcons(burden_with_motif(mut, mutation$jaspar_2014_motif_q0.1_dummy)$mut, mutation$fitcons, 0.125, max(mutation$fitcons[,6]))$burden
#     f = burden_with_phylop(mut, mutation$phylop, 1, max(mutation$phylop[,6]))$burden
#     g = burden_with_phylop(mut, mutation$phylop, 2, max(mutation$phylop[,6]))$burden
#     h = burden_with_geneset(mut, coding_tada[coding_tada$qvalue.combined<0.3,2])$burden
#     i = burden_with_geneset(mut, coding_tada[coding_tada$qvalue.combined<0.1,2])$burden
#     j = burden_with_geneset(mut,Brain_GO_gene[,1])$burden
#     k = burden_with_geneset(mut, brainspan$gene_qt1[,1])$burden
#     l = burden_with_geneset(mut, brainspan$gene_qt4[,1])$burden
#     m = burden_with_CADD(mut, CADD_score_with_index, 12.93,53)$burden
#     n = burden_with_gerp(mut, mutation$gerp, 2, 6.17)$burden
#     o = burden_with_motif_alt(mut, motif_brk_or_gain_index[,1])$burden
#     p = burden_with_geneset(mut, intersect(c(brainspan$gene_qt4[,1],brainspan$gene_qt3[,1]),Brain_GO_gene[,1]))$burden
#     q = burden_with_fitcons_or_motif(mut, mutation$fitcons, 0.125, max(mutation$fitcons[,6]),mutation$jaspar_2014_motif_q0.1_dummy)$burden
#     output = data.frame(rbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q))
#     rownames(output) = c("epigenomics", "+Jinyu_genes", "+motif", "+fitcons_0.125", "+motif+fitcons_0.125", "+phylop_1", "+phylop_2", "+coding_TADA_0.3", "+coding_TADA_0.1",
#                          "Brain_GO_genes", "lowest qt brain exp", "highest qt brain exp", "CADD>12.93", "gerp>2", "motif_br_gain",
#                          "top50% exp genes + brain GO", "motif or fitcons_0.125")
#     colnames(output) = c("ASD", "control", "burden")
#   }
#   else if (mut_type == "coding"){
#     output = show_burden_for_coding(mut, mutation$exon_mut_type)
#   }
#   output 
# }


############ function of show_burden, with some functions deactivated ####################################################################
coding_tada = read.delim("../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
Jinyu_psy_gene = read.delim("../other_annotation/gene_list/151110_Wu_all_gene.genelist", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
Brain_GO_gene = read.delim("../other_annotation/gene_list/160210_GO_brain_genes.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
constraint_gene = read.delim("../other_annotation/gene_list/samocha_NG_2014_top5_percent_genelist_z3.09.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
constraint_gene_3.72 = read.delim("../other_annotation/gene_list/samocha_NG_2014_z3.72_genelist.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#the list below applied certain filters as shown in tableS2 of the corresponding paper, aiming to reduce the false positive calls. 
constraint_gene_tableS2 = read.delim("../other_annotation/gene_list/samocha_NG_2014_constraint_genelist_as_in_tableS2.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

Yuen_gene = read.delim("../other_annotation/gene_list/yuen_NG_2015_genelist.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
FMRP_gene = read.delim("../other_annotation/gene_list/Darnell_Cell_2011_FMRP_targets_genelist.txt", header = FALSE, sep = "\t",stringsAsFactors = FALSE)
Betancur_gene = read.delim("../other_annotation/gene_list/betancur_ASD_genelist.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
SFRAI_gene = read.delim("../other_annotation/gene_list/sfrai_curated_ASD_associated_genelist.txt",header = FALSE, sep = "\t", stringsAsFactors = FALSE)
Kenny_brain_gene = read.delim("../other_annotation/gene_list/betancur_ASD_genelist.txt",header = FALSE, sep = "\t", stringsAsFactors = FALSE)
Purcell_NG_2014_composite_gene = read.delim("../other_annotation/gene_list/Purcell_NG_2014_composite_set_genelist.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
Purcell_calcium_channel = read.delim("../other_annotation/gene_list/Purcell_calcium_channel_genelist.txt",header = FALSE, sep = "\t", stringsAsFactors = FALSE)
AutismKB_gene = read.delim("../other_annotation/gene_list/AutismKB_genelist.txt",header = FALSE, sep = "\t", stringsAsFactors = FALSE)
Irimia_neuron_AS_gene = read.delim("../other_annotation/gene_list/Irimia_cell_2014_neuron_specific_alternative_splicing_genelist.txt",header = FALSE, sep = "\t", stringsAsFactors = FALSE)
Petrovski_RVIS_top25_pct_gene = read.delim("../other_annotation/gene_list/Petrovski_plosgen_RVIS_score_top25_pct_genelist.txt", header = FALSE, sep = "\t",
                                           stringsAsFactors = FALSE)
Petrovski_RVIS_top10_pct_gene = read.delim("../other_annotation/gene_list/Petrovski_plosgen_RVIS_score_top10_pct_genelist.txt", header = FALSE, sep = "\t",
                                           stringsAsFactors = FALSE)
Petrovski_RVIS_top5_pct_gene = read.delim("../other_annotation/gene_list/Petrovski_plosgen_RVIS_score_top5_pct_genelist.txt", header = FALSE, sep = "\t",
                                           stringsAsFactors = FALSE)

Petrovski_RVIS_full_table = read.delim("../other_annotation/gene_list/Petrovski_plosgen_RVIS_percentile_full_table.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
Petrovski_RVIS_bottom_10pct_gene = Petrovski_RVIS_full_table[Petrovski_RVIS_full_table[,3] >=quantile(Petrovski_RVIS_full_table[,3],seq(0,1,0.05))[19],]
Petrovski_RVIS_bottom_5pct_gene = Petrovski_RVIS_full_table[Petrovski_RVIS_full_table[,3] >=quantile(Petrovski_RVIS_full_table[,3],seq(0,1,0.05))[20],]
Petrovski_RVIS_bottom_20pct_gene = Petrovski_RVIS_full_table[Petrovski_RVIS_full_table[,3] >=quantile(Petrovski_RVIS_full_table[,3],seq(0,1,0.05))[17],]

Petrovski_LoF_control_gene = read.delim("../other_annotation/gene_list/Petrovski_LoF_gene_control_set.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

Pinto_ID_genelist = read.delim("../other_annotation/gene_list/Pinto_AJHG_ID_genelist.txt",header = FALSE, sep = "\t", stringsAsFactors = FALSE)
Huang_Hpi_score = read.delim("../other_annotation/gene_list/Huang_Plosgen_haploinsufficiency_score_by_gene.txt", 
                             header = FALSE, sep = "\t",stringsAsFactors = FALSE)
Bayes_postsynaptic_genelist = read.delim("../other_annotation/gene_list/Bayes_postsynaptic_genelist.txt",
                                         header = FALSE, sep = "\t", stringsAsFactors = FALSE)
Bayes_postsynaptic_high_confidence_genelist = read.delim("../other_annotation/gene_list/Bayes_postsynaptic_high_confidence_genelist.txt",
                                                         header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#Petrovski_haploinsufficiency_genelist = read.delim("../other_annotation/gene_list/Petrovski_plosgen_haploinsufficiency_genetlist.txt",
                                                  # header = FALSE, sep = "\t", stringsAsFactors = FALSE)
Petrovski_haploinsufficiency_genelist = read.delim("../other_annotation/gene_list/Petrovski_plosgen_haploinsufficiency_including_all_without_ncscore_genelist.txt",
                                                    header = FALSE, sep = "\t", stringsAsFactors = FALSE)
ncRVIS = read.delim("../other_annotation/gene_list/160708_ncRVIS_score_per_gene.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
ncGerp = read.delim("../other_annotation/gene_list/Petrovski_plosgen_RVIS_ncGERP_table.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
DDG2P = read.delim("../other_annotation/gene_list/DDG2P_2013_november.genelist", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
DDG2P_confirmed_DD = read.delim("../other_annotation/gene_list/DDG2P_2013_november_confirmed_DD_gene.genelist", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

Olfactory_gene = read.delim("../other_annotation/gene_list/Olfactory_genelist.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
SFRAI_gene_syndromic = read.delim("../other_annotation/gene_list/sfrai_genes_sydromic_genelist.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
SFRAI_gene_high_confidence = read.delim("../other_annotation/gene_list/sfrai_genes_high_confidence_genelist.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
SFRAI_gene_strong_candidate = read.delim("../other_annotation/gene_list/sfrai_genes_strong_candidate_genelist.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
SFRAI_gene_suggestive_evidence = read.delim("../other_annotation/gene_list/sfrai_genes_suggestive_evidence_genelist.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
SFRAI_gene_hypothesized_gene = read.delim("../other_annotation/gene_list/sfari_genes_hypothesized_genelist.txt",header = FALSE, sep = "\t", stringsAsFactors = FALSE )
FMRP_intersect_Brain_GO = intersect(FMRP_gene[,1],Brain_GO_gene[,1])
FMRP_intersect_high_brain_exp = intersect(FMRP_gene[,1],brainspan$gene_qt4[,1])
FMRP_intersect_Brain_GO_or_high_brain_exp = intersect(FMRP_gene[,1],c(Brain_GO_gene[,1],brainspan$gene_qt4[,1]))
FRMP_intersect_brain_GO_intersect_top50_pct_exp = intersect(FMRP_gene[,1],intersect(Brain_GO_gene[,1],c(brainspan$gene_qt3[,1],brainspan$gene_qt4[,1])))
FRMP_intersect_brain_GO_or_top50_pct_exp = intersect(FMRP_gene[,1],c(Brain_GO_gene[,1],brainspan$gene_qt3[,1],brainspan$gene_qt4[,1]))
FRMP_intersect_brain_GO_intersect_top75_pct_exp = intersect(FMRP_gene[,1],intersect(Brain_GO_gene[,1],c(brainspan$gene_qt3[,1],brainspan$gene_qt4[,1], brainspan$gene_qt2[,1])))

DAWN_old = read.delim("../other_annotation/gene_list/DAWN_old_q0.05_genelist.txt", header = FALSE, sep = "\t",stringsAsFactors = FALSE)
DAWN_new = read.delim("../other_annotation/gene_list/DAWN_new_q0.05_genelist.txt", header = FALSE, sep = "\t",stringsAsFactors = FALSE)

ExAC_gene = readRDS("../other_annotation/gene_list/ConstraintMat.RDS")

Nick_ASD_gene = read.delim("../other_annotation/gene_list/Nick_machine_learning_ASD_gene_data.txt",
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)


union_1 = c(Jinyu_psy_gene[,1],Betancur_gene[,1],Purcell_NG_2014_composite_gene[,1],AutismKB_gene[,1],SFRAI_gene_high_confidence[,1],SFRAI_gene_strong_candidate[,1],
            Pinto_ID_genelist[,1],Bayes_postsynaptic_high_confidence_genelist[,1],coding_tada[coding_tada$qvalue.combined<0.3,2])

union_2 = c(Betancur_gene[,1],AutismKB_gene[,1],Purcell_NG_2014_composite_gene[,1],SFRAI_gene_high_confidence[,1],SFRAI_gene_strong_candidate[,1],
            Pinto_ID_genelist[,1],Bayes_postsynaptic_high_confidence_genelist[,1],coding_tada[coding_tada$qvalue.combined<0.3,2])

union_3 = c(Betancur_gene[,1],AutismKB_gene[,1],SFRAI_gene_high_confidence[,1],SFRAI_gene_strong_candidate[,1],coding_tada[coding_tada$qvalue.combined<0.3,2])

#stringent_ASD = c(coding_tada[coding_tada$qvalue.combined<0.1,2],SFRAI_gene_high_confidence[,1],SFRAI_gene_strong_candidate[,1],AutismKB_gene[,1])
#stringent_ASD = c(coding_tada[coding_tada$qvalue.combined<0.2,2],SFRAI_gene_high_confidence[,1],SFRAI_gene_strong_candidate[,1],AutismKB_gene[,1])

#stringent_ASD = c(coding_tada[coding_tada$qvalue.combined<0.1,2],SFRAI_gene_high_confidence[,1],SFRAI_gene_strong_candidate[,1])
#stringent_ASD = c(coding_tada[coding_tada$qvalue.combined<0.3,2],SFRAI_gene_high_confidence[,1],SFRAI_gene_strong_candidate[,1],AutismKB_gene[,1],Betancur_gene[,1] )

stringent_ASD = c(coding_tada[coding_tada$qvalue.combined<0.3,2],SFRAI_gene_high_confidence[,1],SFRAI_gene_strong_candidate[,1],AutismKB_gene[,1],Betancur_gene[,1])

#relaxed_ASD = c(coding_tada[coding_tada$qvalue.combined<0.5,2],SFRAI_gene_high_confidence[,1],SFRAI_gene_strong_candidate[,1],Pinto_ID_genelist[,1],
             #   Bayes_postsynaptic_high_confidence_genelist[,1],Purcell_NG_2014_composite_gene[,1],AutismKB_gene[,1],Betancur_gene[,1])

#relaxed_ASD genes are named neuronpsychiatric genes
relaxed_ASD = c(coding_tada[coding_tada$qvalue.combined<0.5,2],SFRAI_gene_high_confidence[,1],SFRAI_gene_strong_candidate[,1],Pinto_ID_genelist[,1],
                Bayes_postsynaptic_high_confidence_genelist[,1],Purcell_NG_2014_composite_gene[,1],AutismKB_gene[,1],Betancur_gene[,1],
                SFRAI_gene_syndromic[,1], SFRAI_gene_suggestive_evidence[,1])

#new definition of putative ASD 
#relaxed_ASD = c(coding_tada[coding_tada$qvalue.combined<0.3,2],SFRAI_gene_high_confidence[,1],SFRAI_gene_strong_candidate[,1],Pinto_ID_genelist[,1],
#                     AutismKB_gene[,1],Betancur_gene[,1],SFRAI_gene_syndromic[,1],SFRAI_gene_suggestive_evidence[,1], SFRAI_gene_hypothesized_gene[,1])

nonASD_genes = coding_tada[coding_tada$qvalue.combined > quantile(coding_tada$qvalue.combined, seq(0,1,0.1))[10],]$TadaName

constraint_union = c(Petrovski_RVIS_top5_pct_gene[,1],Petrovski_haploinsufficiency_genelist[,1],Huang_Hpi_score[Huang_Hpi_score[,2] >= 0.95,1])
constraint_intersection = intersect(Petrovski_RVIS_top5_pct_gene[,1],intersect(Petrovski_haploinsufficiency_genelist[,1],Huang_Hpi_score[Huang_Hpi_score[,2] >= 0.75,1]))
constraint_union2 = c(Petrovski_RVIS_top5_pct_gene[,1],Petrovski_haploinsufficiency_genelist[,1],Huang_Hpi_score[Huang_Hpi_score[,2] >= 0.95,1],constraint_gene_tableS2[,1])

#nonconstraint_union = setdiff(c(Petrovski_RVIS_bottom_5pct_gene[,1],Huang_Hpi_score[Huang_Hpi_score[,2] <= 0.05,1], Petrovski_LoF_control_gene[,1]),constraint_union)

#nonconstraint_union = setdiff(c(Petrovski_RVIS_bottom_20pct_gene[,1],Huang_Hpi_score[Huang_Hpi_score[,2] <= 0.2,1], Petrovski_LoF_control_gene[,1]),constraint_union)

nonconstraint_union = setdiff(c(Petrovski_RVIS_bottom_10pct_gene[,1],Huang_Hpi_score[Huang_Hpi_score[,2] <= 0.1,1], Petrovski_LoF_control_gene[,1]),constraint_union)

show_burden <-function(mut_file_name, input_type = "file", mut_type = "noncoding"){
  #input_type could be file or a table matching regions with mutations. 
  # mut_type could be "noncoding" or "coding
  if(input_type == "file"){
    mut = read.delim(mut_file_name, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  }
  else if (input_type == "table"){
    mut = mut_file_name
  }
  if(mut_type == "noncoding"){
    a = burden(mut)
    b = burden_with_geneset(mut, Jinyu_psy_gene[,1])$burden
    c = burden_with_motif(mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    d = burden_with_fitcons(mut, mutation$fitcons, 0.125, max(mutation$fitcons[,6]))$burden
    e = burden_with_fitcons(burden_with_motif(mut, mutation$jaspar_2014_motif_q0.1_dummy)$mut, mutation$fitcons, 0.125, max(mutation$fitcons[,6]))$burden
    f = burden_with_phylop(mut, mutation$phylop, 1, max(mutation$phylop[,6]))$burden
    g = burden_with_phylop(mut, mutation$phylop, 2, max(mutation$phylop[,6]))$burden
    h = burden_with_geneset(mut, coding_tada[coding_tada$qvalue.combined<0.3,2])$burden
    i = burden_with_geneset(mut, coding_tada[coding_tada$qvalue.combined<0.1,2])$burden
    j = burden_with_geneset(mut,Brain_GO_gene[,1])$burden
    k = burden_with_geneset(mut, brainspan$gene_qt1[,1])$burden
    l = burden_with_geneset(mut, brainspan$gene_qt4[,1])$burden
    m = burden_with_CADD(mut, CADD_score_with_index, 13.71,53)$burden
    n = burden_with_gerp(mut, mutation$gerp, 2, 6.17)$burden
    o = burden_with_motif_alt(mut, mutation$motif_brk_or_gain_index)$burden
    p = burden_with_geneset(mut, intersect(c(brainspan$gene_qt4[,1],brainspan$gene_qt3[,1]),Brain_GO_gene[,1]))$burden
    q = burden_with_fitcons_or_motif(mut, mutation$fitcons, 0.125, max(mutation$fitcons[,6]),mutation$jaspar_2014_motif_q0.1_dummy)$burden
    r = burden_with_geneset(mut, coding_tada[coding_tada$qvalue.combined<0.5,2])$burden
    s = burden_with_geneset(mut, constraint_gene[,1])$burden
    t = burden_with_geneset(mut, Yuen_gene[,1])$burden
    u = burden_with_geneset(burden_with_motif(mut, mutation$jaspar_2014_motif_q0.1_dummy)$mut,coding_tada[coding_tada$qvalue.combined<0.5,2])$burden
    v = burden_with_geneset(burden_with_CADD(mut, CADD_score_with_index, 13.71,53)$mut,coding_tada[coding_tada$qvalue.combined<0.5,2])$burden
    w = burden_with_geneset(burden_with_CADD_or_motif(mut, CADD_score_with_index, 13.71, 53,mutation$jaspar_2014_motif_q0.1_dummy)$mut, coding_tada[coding_tada$qvalue.combined<0.5,2])$burden
    x = burden_with_geneset(burden_with_CADD_or_motif(mut, CADD_score_with_index, 13.71, 53,mutation$jaspar_2014_motif_q0.1_dummy)$mut, constraint_gene[,1])$burden
    y = burden_with_geneset(mut, FMRP_gene[,1])$burden
    z = burden_with_geneset(mut, Betancur_gene[,1])$burden
    a1 = burden_with_geneset(mut, FMRP_intersect_Brain_GO)$burden
    b1 = burden_with_geneset(mut, FMRP_intersect_Brain_GO_or_high_brain_exp)$burden
    c1 = burden_with_geneset(mut, FMRP_intersect_high_brain_exp)$burden
    d1 = burden_with_motif(burden_with_geneset(mut, FMRP_gene[,1])$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    e1 = burden_with_motif(burden_with_geneset(mut, FMRP_intersect_Brain_GO)$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    f1 = burden_with_motif(burden_with_geneset(mut, FMRP_intersect_Brain_GO_or_high_brain_exp)$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    g1 = burden_with_motif(burden_with_geneset(mut, FMRP_intersect_high_brain_exp)$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    h1 = burden_with_geneset(burden_with_CADD_or_motif(mut, CADD_score_with_index, 13.71, 53,mutation$jaspar_2014_motif_q0.1_dummy)$mut,FMRP_gene[,1])$burden
    i1 = burden_with_geneset(burden_with_CADD_or_motif(mut, CADD_score_with_index, 13.71, 53,mutation$jaspar_2014_motif_q0.1_dummy)$mut,FMRP_intersect_Brain_GO)$burden
    j1 = burden_with_geneset(burden_with_CADD_or_motif(mut, CADD_score_with_index, 13.71, 53,mutation$jaspar_2014_motif_q0.1_dummy)$mut,FMRP_intersect_Brain_GO_or_high_brain_exp)$burden
    k1 = burden_with_geneset(burden_with_CADD_or_motif(mut, CADD_score_with_index, 13.71, 53,mutation$jaspar_2014_motif_q0.1_dummy)$mut,FMRP_intersect_high_brain_exp)$burden
    l1 = burden_with_motif(burden_with_geneset(mut,FRMP_intersect_brain_GO_intersect_top50_pct_exp)$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    m1 = burden_with_motif(burden_with_geneset(mut,FRMP_intersect_brain_GO_or_top50_pct_exp)$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    n1 = burden_with_motif(burden_with_geneset(mut,FRMP_intersect_brain_GO_intersect_top75_pct_exp)$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    o1 = burden_with_motif(burden_with_geneset(mut, constraint_gene[,1])$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    p1 = burden_with_motif(burden_with_geneset(mut, intersect(constraint_gene[,1],Brain_GO_gene[,1]))$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    q1 = burden_with_motif(burden_with_geneset(mut, constraint_gene_3.72[,1])$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    r1 = burden_with_motif(burden_with_geneset(mut, intersect(constraint_gene_3.72[,1],Brain_GO_gene[,1]))$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    s1 = burden_with_motif(burden_with_geneset(mut,Betancur_gene[,1])$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    t1 = burden_with_motif(burden_with_geneset(mut, Jinyu_psy_gene[,1])$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    u1 = burden_with_geneset(burden_with_gerp(mut, mutation$gerp, 2, 6.17)$mut,coding_tada[coding_tada$qvalue.combined<0.5,2])$burden
    v1 = burden_with_geneset(burden_with_gerp(mut, mutation$gerp, 2, 6.17)$mut,Jinyu_psy_gene[,1])$burden
    w1 = burden_with_geneset(burden_with_gerp(mut, mutation$gerp, 2, 6.17)$mut,FMRP_intersect_Brain_GO)$burden
    x1 = burden_with_geneset(burden_with_phylop(mut, mutation$phylop, 2, max(mutation$phylop[,6]))$mut,coding_tada[coding_tada$qvalue.combined<0.5,2])$burden
    y1 = burden_with_geneset(burden_with_phylop(mut, mutation$phylop, 2, max(mutation$phylop[,6]))$mut,Jinyu_psy_gene[,1])$burden
    z1 = burden_with_geneset(burden_with_phylop(mut, mutation$phylop, 2, max(mutation$phylop[,6]))$mut,FMRP_intersect_Brain_GO)$burden
    a2 = burden_with_motif(burden_with_geneset(mut, brainspan$gene_qt4[,1])$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    b2 = burden_with_motif(burden_with_geneset(mut, intersect(Brain_GO_gene[,1],c(brainspan$gene_qt4[,1],brainspan$gene_qt3[,1])))$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    c2 = burden_with_motif(burden_with_geneset(mut, intersect(constraint_gene[,1],FMRP_gene[,1]))$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    d2 = burden_with_geneset(mut, SFRAI_gene[,1])$burden
    e2 = burden_with_geneset(mut, Kenny_brain_gene[,1])$burden
    f2 = burden_with_motif(burden_with_geneset(mut, Kenny_brain_gene[,1])$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    g2 = burden_with_geneset(mut, Purcell_NG_2014_composite_gene[,1])$burden
    h2 = burden_with_motif(burden_with_geneset(mut, Purcell_NG_2014_composite_gene[,1])$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    i2 = burden_with_CADD(burden_with_geneset(mut,Purcell_NG_2014_composite_gene[,1])$mut, CADD_score_with_index, 13.71,53)$burden
    j2 = burden_with_geneset(mut,Purcell_calcium_channel[,1])$burden
    k2 = burden_with_geneset(mut, AutismKB_gene[,1])$burden
    l2 = burden_with_motif(burden_with_geneset(mut, AutismKB_gene[,1])$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    m2 = burden_with_geneset(mut, SFRAI_gene_syndromic[,1])$burden
    n2 = burden_with_geneset(mut, SFRAI_gene_high_confidence[,1])$burden
    o2 = burden_with_geneset(mut, SFRAI_gene_strong_candidate[,1])$burden
    p2 = burden_with_geneset(mut, SFRAI_gene_suggestive_evidence[,1])$burden
    q2 = burden_with_geneset(mut, Irimia_neuron_AS_gene[,1])$burden
    r2 = burden_with_motif(burden_with_geneset(mut, Irimia_neuron_AS_gene[,1])$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    s2 = burden_with_geneset(mut, Olfactory_gene[,1])$burden
    t2 = burden_with_motif(burden_with_geneset(mut, Olfactory_gene[,1])$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    u2 = burden_with_geneset(mut,Petrovski_RVIS_top25_pct_gene[,1])$burden
    v2 = burden_with_motif(burden_with_geneset(mut,Petrovski_RVIS_top25_pct_gene[,1])$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    w2 = burden_with_geneset(mut,Petrovski_RVIS_top10_pct_gene[,1])$burden
    x2 = burden_with_motif(burden_with_geneset(mut,Petrovski_RVIS_top10_pct_gene[,1])$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    y2 = burden_with_geneset(mut, Pinto_ID_genelist[,1])$burden
    z2 = burden_with_motif(burden_with_geneset(mut, Pinto_ID_genelist[,1])$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    a3 = burden_with_hpi(mut, Huang_Hpi_score, 0.15)$burden
    b3 = burden_with_hpi(mut, Huang_Hpi_score, 0.35)$burden
    c3 = burden_with_hpi(mut, Huang_Hpi_score, 0.55)$burden
    d3 = burden_with_hpi(mut, Huang_Hpi_score, 0.65)$burden
    e3 = burden_with_hpi(mut, Huang_Hpi_score, 0.75)$burden
    f3 = burden_with_motif(burden_with_hpi(mut, Huang_Hpi_score, 0.15)$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    g3 = burden_with_motif(burden_with_hpi(mut, Huang_Hpi_score, 0.35)$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    h3 = burden_with_motif(burden_with_hpi(mut, Huang_Hpi_score, 0.55)$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    i3 = burden_with_motif(burden_with_hpi(mut, Huang_Hpi_score, 0.65)$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    j3 = burden_with_motif(burden_with_hpi(mut, Huang_Hpi_score, 0.75)$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    k3 = burden_with_geneset(mut, Bayes_postsynaptic_genelist[,1])$burden
    l3 = burden_with_motif(burden_with_geneset(mut, Bayes_postsynaptic_genelist[,1])$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    m3 = burden_with_geneset(mut, Bayes_postsynaptic_high_confidence_genelist[,1])$burden
    n3 = burden_with_motif(burden_with_geneset(mut, Bayes_postsynaptic_high_confidence_genelist[,1])$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    o3 = burden_with_geneset(mut,union_1)$burden
    p3 = burden_with_motif(burden_with_geneset(mut,union_1)$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    q3 = burden_with_geneset(mut, Petrovski_haploinsufficiency_genelist[,1])$burden
    r3 = burden_with_motif(burden_with_geneset(mut, Petrovski_haploinsufficiency_genelist[,1])$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    s3 = burden_with_geneset(mut, constraint_gene_tableS2[,1])$burden
    t3 = burden_with_motif(burden_with_geneset(mut, constraint_gene_tableS2[,1])$mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    u3 = burden_with_geneset(mut, union_2)$burden
    v3 = burden_with_motif(burden_with_geneset(mut,union_2)$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    w3 = burden_with_geneset(mut, constraint_union)$burden
    x3 = burden_with_motif(burden_with_geneset(mut, constraint_union)$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    y3 = burden_with_geneset(mut, constraint_intersection)$burden
    z3 = burden_with_motif(burden_with_geneset(mut, constraint_intersection)$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    a4 = burden_with_geneset(mut, union_3)$burden
    b4 = burden_with_motif(burden_with_geneset(mut,union_3)$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    c4 = burden_with_geneset(mut, constraint_union2)$burden
    d4 = burden_with_motif(burden_with_geneset(mut, constraint_union2)$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    e4 = burden_with_nick_gene(mut,Nick_ASD_gene,0,0.8)$burden #assume meanp cutoff is 0.8
    f4 = burden_with_motif(burden_with_nick_gene(mut,Nick_ASD_gene,0,0.8)$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    g4 = burden_with_nick_gene(mut,Nick_ASD_gene,0.62,0)$burden #assume minp cutoff is 0.62 for high sensitivity mode
    h4 = burden_with_motif(burden_with_nick_gene(mut,Nick_ASD_gene,0.62,0)$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    i4 = burden_with_eigen(mut, mutation$mut_eigen, 0.25304484)$burden
    output = data.frame(rbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,l1,m1,n1,o1,p1,q1,r1,s1,t1,u1,v1,w1,x1,y1,z1,
                              a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,k2,l2,m2,n2,o2,p2,q2,r2,s2,t2,u2,v2,w2,x2,y2,z2,
                              a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,k3,l3,m3,n3,o3,p3,q3,r3,s3,t3,u3,v3,w3,x3,y3,z3,
                              a4,b4,c4,d4,e4,f4,g4,h4,i4))
    rownames(output) = c("epigenomics", "+Jinyu_genes", "+motif", "+fitcons_0.125", "+motif+fitcons_0.125", "+phylop_1", "+phylop_2", "+coding_TADA_0.3", "+coding_TADA_0.1",
                         "Brain_GO_genes", "lowest qt brain exp", "highest qt brain exp", "CADD>12.93", "gerp>2","motif_breaker_funseq2",
                         "top50% exp genes + brain GO", "motif or fitcons_0.125", "+coding_TADA_0.5", "+constraint_top5%", "+yuen_NG_2014_gene", "+motif_TADA_q0.5", "+CADD_TADA_q0.5",
                         "+motif_or_CADD_TADA_q0.5", "+motif_or_CADD_constraint_top5%","+FMRP_gene", "+Betancur_gene","FMRP_in_brainGO",
                         "FMRP_in_brainGO_or_high_exp", "FMRP_in_high_exp", "motif+FMRP_gene", "motif+FMRP_in_brainGO","motif+FMRP_in_brainGO_or_high_exp","motif+FMRP_in_high_exp",
                         "CADD_or_motif+FMRP_gene", "CADD_or_motif+FMRP_in_brainGO","CADD_or_motif+FMRP_in_brainGO_or_high_exp","CADD_or_motif+FMRP_in_high_exp",
                         "motif_AND_FMRP_AND_top50%exp_AND_brainGO", "motif_AND_FRMP_in_top50%exp_or_brainGO","motif_AND_FMRP_AND_top75%exp_AND_brainGO",
                         "motif+5%constraint", "motif+5%constraint+brainGO","motif+2.5%constraint", "motif+2.5%constraint+brainGO",
                         "motif+Betancur","motif_Jinyu",
                         "gerp2+TADA_q0.5","gerp2+Jinyu","gerp2+FMRP_intersect_brainGO",
                         "phylop2+TADA_q0.5","phylop2+Jinyu","phylop2+FMRP_intersect_brainGO",
                         "motif+top25%exp","motif+50%exp+brainGO","motif+constraint+FMRP",
                         "SFRAI_gene","Kenny_brain_gene","motif+Kenny_brain_gene",
                         "Purcell_gene","motif+Purcell_gene","CADD+Purcell_gene",
                         "Purcell_calcium_channel","AutismKB","AutismKB+motif",
                         "sfrai_syndromic","sfrai_high_confidence","sfrai_strong_candidate","sfrai_suggestive_evidence",
                         "Irimia_neuron_AS","Irimia_neuron_AS_motif",
                         "Olfactory_gene","Olfactory_gene+motif",
                         "RVIS_top25pct","RVIS_top25pct+motif",
                         "RVIS_top10pct","RVIS_top10pct+motif",
                         "Pinto_ID_gene","Pinto_ID_gene+motif",
                         "Hpi_0.15","Hpi_0.35","Hpi_0.55","Hpi_0.65","Hpi_0.75","Hpi_0.15+motif","Hpi_0.35+motif","Hpi_0.55+motif","Hpi_0.65+motif","Hpi_0.75+motif",
                         "postsynaptic_genes","postsynaptic_genes+motif","postsynaptic_hconf","postsynaptic_hconf+motif",
                         "union_1","union_1+motif",
                         "Petrovski_haploinsufficiency","Petrovski_haploinsufficiency+motif",
                         "constraint_tableS2", "constriant_tableS2+motif",
                         "union2", "union2+motif",
                         "constraint_union", "constraint_union+motif", "constraint_intersection", "constraint_intersection+motif",
                         "union3","union3+motif",
                         "constraint_union2","constraint_union2+motif",
                         "Nick_ASD_meanp0.8", "Nick_ASD_meanp0.8+motif",
                         "Nick_ASD_minp0.62", "Nick_ASD_minp0.62+motif",
                         "mut_eigen_top10%")
    colnames(output) = c("ASD", "control", "burden")
  }
  else if (mut_type == "coding"){
    output = show_burden_for_coding(mut, mutation$exon_mut_type)
  }
  output 
}


generate_mut_enhancer_pair_for_burden_analysis <-function(bed_file,epigenomic_file_name){
  #bed files is the mutation file
  #epigenomic_file_name is the name of the epigenomic partitions, which are in ../other_annotation/epigenomic_annotation/
  command = paste("bedtools intersect -a ",bed_file, " -b ", epigenomic_file_name," -wa -wb | awk {'print $4\"\t\"$8'} | sort | uniq ",sep = "")
  temp = system(command, intern = TRUE)
  temp = as.data.frame(t(mapply(function(x) strsplit(x,split="\t")[[1]],temp)))
  temp
}

table_colname = c("ASD","control","frquency_ratio","burden","odds_ratio","pvalue","lower_bound","upper_bound")


show_burden_for_selected_features <-function(mut_file_name, input_type = "file", mut_type = "noncoding",ref = "self", syn_A = length(mutation$ASD_effective_SNV_ID), syn_C = length(mutation$control_effective_SNV_ID)){
  #input_type could be file or a table matching regions with mutations. 
  # mut_type could be "noncoding" or "coding
  # ref = "total_ratio" means using total number of SNVs as background, and use fisher.test to 
  # get pvalues
  # syn_A is the number of SNVs in ASD
  # syn_C is the number of SNVs in control
  if(input_type == "file"){
    mut = read.delim(mut_file_name, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  }
  else if (input_type == "table"){
    mut = mut_file_name
  }
  if(mut_type == "noncoding"){
    baseline = burden(mut)
    Tada_q0.1 = burden_with_geneset(mut, coding_tada[coding_tada$qvalue.combined<0.1,2])$burden
    Tada_q0.2 = burden_with_geneset(mut, coding_tada[coding_tada$qvalue.combined<0.2,2])$burden
    Tada_q0.3 = burden_with_geneset(mut, coding_tada[coding_tada$qvalue.combined<0.3,2])$burden
    Tada_q0.4 = burden_with_geneset(mut, coding_tada[coding_tada$qvalue.combined<0.4,2])$burden
    Tada_q0.5 = burden_with_geneset(mut, coding_tada[coding_tada$qvalue.combined<0.5,2])$burden
    Tada_q0.1_motif = burden_with_motif(burden_with_geneset(mut,coding_tada[coding_tada$qvalue.combined<0.1,2])$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    Tada_q0.2_motif = burden_with_motif(burden_with_geneset(mut,coding_tada[coding_tada$qvalue.combined<0.2,2])$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    Tada_q0.3_motif = burden_with_motif(burden_with_geneset(mut,coding_tada[coding_tada$qvalue.combined<0.3,2])$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    Tada_q0.4_motif = burden_with_motif(burden_with_geneset(mut,coding_tada[coding_tada$qvalue.combined<0.4,2])$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    Tada_q0.5_motif = burden_with_motif(burden_with_geneset(mut,coding_tada[coding_tada$qvalue.combined<0.5,2])$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    motif = burden_with_motif(mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    ASD_stringent = burden_with_geneset(mut, stringent_ASD)$burden
    ASD_stringent_motif = burden_with_motif(burden_with_geneset(mut,stringent_ASD)$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    ASD_relaxed = burden_with_geneset(mut, relaxed_ASD)$burden
    ASD_relaxed_motif = burden_with_motif(burden_with_geneset(mut,relaxed_ASD)$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    nonASD = burden_with_geneset(mut, nonASD_genes)$burden
    nonASD_motif = burden_with_motif(burden_with_geneset(mut, nonASD_genes)$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    constraint_union_data1 = burden_with_geneset(mut, constraint_union)$burden
    constraint_union_data1_motif = burden_with_motif(burden_with_geneset(mut, constraint_union)$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    CADD = burden_with_CADD(mut, CADD_score_with_index, as.numeric(quantile(mutation$CADD_score[,1], c(0.95))),as.numeric(quantile(mutation$CADD_score[,1], c(1))))$burden
    CADD_90pct = burden_with_CADD(mut, CADD_score_with_index, as.numeric(quantile(mutation$CADD_score[,1], c(0.9))),as.numeric(quantile(mutation$CADD_score[,1], c(1))))$burden
    ASD_stringent_CADD = burden_with_CADD(burden_with_geneset(mut,stringent_ASD)$mut, CADD_score_with_index, as.numeric(quantile(mutation$CADD_score[,1], c(0.95))),as.numeric(quantile(mutation$CADD_score[,1], c(1))))$burden
    ASD_stringent_CADD_90pct = burden_with_CADD(burden_with_geneset(mut,stringent_ASD)$mut, CADD_score_with_index, as.numeric(quantile(mutation$CADD_score[,1], c(0.9))),as.numeric(quantile(mutation$CADD_score[,1], c(1))))$burden
    ASD_relaxed_CADD = burden_with_CADD(burden_with_geneset(mut,relaxed_ASD)$mut, CADD_score_with_index, as.numeric(quantile(mutation$CADD_score[,1], c(0.95))),as.numeric(quantile(mutation$CADD_score[,1], c(1))))$burden
    ASD_relaxed_CADD_90pct = burden_with_CADD(burden_with_geneset(mut,relaxed_ASD)$mut, CADD_score_with_index, as.numeric(quantile(mutation$CADD_score[,1], c(0.9))),as.numeric(quantile(mutation$CADD_score[,1], c(1))))$burden
    nonASD_CADD = burden_with_CADD(burden_with_geneset(mut,nonASD_genes)$mut, CADD_score_with_index, as.numeric(quantile(mutation$CADD_score[,1], c(0.95))),as.numeric(quantile(mutation$CADD_score[,1], c(1))))$burden
    nonASD_CADD_90pt = burden_with_CADD(burden_with_geneset(mut,nonASD_genes)$mut, CADD_score_with_index, as.numeric(quantile(mutation$CADD_score[,1], c(0.9))),as.numeric(quantile(mutation$CADD_score[,1], c(1))))$burden
    constraint_CADD = burden_with_CADD(burden_with_geneset(mut,constraint_union)$mut, CADD_score_with_index, as.numeric(quantile(mutation$CADD_score[,1], c(0.95))),as.numeric(quantile(mutation$CADD_score[,1], c(1))))$burden
    constraint_CADD_90pct = burden_with_CADD(burden_with_geneset(mut,constraint_union)$mut, CADD_score_with_index, as.numeric(quantile(mutation$CADD_score[,1], c(0.9))),as.numeric(quantile(mutation$CADD_score[,1], c(1))))$burden
    #constraint_union_data2 = burden_with_geneset(mut, constraint_union2)$burden
    #constraint_union_data2_motif = burden_with_motif(burden_with_geneset(mut, constraint_union2)$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    #ASD_predicted = burden_with_nick_gene(mut, Nick_ASD_gene,0,0.8)$burden
    #ASD_predicted_motif = burden_with_motif(burden_with_nick_gene(mut,Nick_ASD_gene,0,0.8)$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    ##ASD_predicted_sensitive = burden_with_nick_gene(mut, Nick_ASD_gene,0.62,0)$burden
    ASD_predicted_motif_sensitive = burden_with_motif(burden_with_nick_gene(mut,Nick_ASD_gene,0.62,0)$mut,mutation$jaspar_2014_motif_q0.1_dummy)$burden
    eigen_gt_95pct = burden_with_eigen(mut, mutation$mut_eigen, as.numeric(quantile(mutation$mut_eigen, c(0.95))))$burden
    eigen_gt_90pct = burden_with_eigen(mut, mutation$mut_eigen, as.numeric(quantile(mutation$mut_eigen, c(0.9))))$burden
    ASD_stringent_eigen_95pct = burden_with_eigen(burden_with_geneset(mut,stringent_ASD)$mut, mutation$mut_eigen, as.numeric(quantile(mutation$mut_eigen, c(0.95))))$burden
    ASD_stringent_eigen_90pct = burden_with_eigen(burden_with_geneset(mut,stringent_ASD)$mut, mutation$mut_eigen, as.numeric(quantile(mutation$mut_eigen, c(0.9))))$burden
    ASD_relaxed_eigen_95pct = burden_with_eigen(burden_with_geneset(mut,relaxed_ASD)$mut, mutation$mut_eigen, as.numeric(quantile(mutation$mut_eigen, c(0.95))))$burden
    ASD_relaxed_eigen_90pct = burden_with_eigen(burden_with_geneset(mut,relaxed_ASD)$mut, mutation$mut_eigen, as.numeric(quantile(mutation$mut_eigen, c(0.9))))$burden
    nonASD_eigen_95pct = burden_with_eigen(burden_with_geneset(mut,nonASD_genes)$mut, mutation$mut_eigen, as.numeric(quantile(mutation$mut_eigen, c(0.95))))$burden
    nonASD_eigen_90pct = burden_with_eigen(burden_with_geneset(mut,nonASD_genes)$mut, mutation$mut_eigen, as.numeric(quantile(mutation$mut_eigen, c(0.9))))$burden
    constraint_eigen_95pct = burden_with_eigen(burden_with_geneset(mut,constraint_union)$mut, mutation$mut_eigen, as.numeric(quantile(mutation$mut_eigen, c(0.95))))$burden
    constraint_eigen_90pct = burden_with_eigen(burden_with_geneset(mut,constraint_union)$mut, mutation$mut_eigen, as.numeric(quantile(mutation$mut_eigen, c(0.9))))$burden
    phylop_gt_95pct = burden_with_phylop(mut, mutation$phylop, as.numeric(quantile(mutation$phylop[,6], c(0.95))), as.numeric(quantile(mutation$phylop[,6], c(1))))$burden
    phylop_gt_90pct = burden_with_phylop(mut, mutation$phylop, as.numeric(quantile(mutation$phylop[,6], c(0.9))), as.numeric(quantile(mutation$phylop[,6], c(1))))$burden
    ASD_stringent_phylop_95pct = burden_with_phylop(burden_with_geneset(mut,stringent_ASD)$mut, mutation$phylop, as.numeric(quantile(mutation$phylop[,6], c(0.95))), as.numeric(quantile(mutation$phylop[,6], c(1))))$burden
    ASD_stringent_phylop_90pct = burden_with_phylop(burden_with_geneset(mut,stringent_ASD)$mut, mutation$phylop, as.numeric(quantile(mutation$phylop[,6], c(0.9))), as.numeric(quantile(mutation$phylop[,6], c(1))))$burden
    ASD_relaxed_phylop_95pct = burden_with_phylop(burden_with_geneset(mut,relaxed_ASD)$mut, mutation$phylop, as.numeric(quantile(mutation$phylop[,6], c(0.95))), as.numeric(quantile(mutation$phylop[,6], c(1))))$burden
    ASD_relaxed_phylop_90pct = burden_with_phylop(burden_with_geneset(mut,relaxed_ASD)$mut, mutation$phylop, as.numeric(quantile(mutation$phylop[,6], c(0.9))), as.numeric(quantile(mutation$phylop[,6], c(1))))$burden
    nonASD_phylop_95pct = burden_with_phylop(burden_with_geneset(mut,nonASD_genes)$mut, mutation$phylop, as.numeric(quantile(mutation$phylop[,6], c(0.95))), as.numeric(quantile(mutation$phylop[,6], c(1))))$burden
    nonASD_phylop_90pct = burden_with_phylop(burden_with_geneset(mut,nonASD_genes)$mut, mutation$phylop, as.numeric(quantile(mutation$phylop[,6], c(0.9))), as.numeric(quantile(mutation$phylop[,6], c(1))))$burden
    constraint_phylop_95pct = burden_with_phylop(burden_with_geneset(mut,constraint_union)$mut, mutation$phylop, as.numeric(quantile(mutation$phylop[,6], c(0.95))), as.numeric(quantile(mutation$phylop[,6], c(1))))$burden
    constraint_phylop_90pct = burden_with_phylop(burden_with_geneset(mut,constraint_union)$mut, mutation$phylop, as.numeric(quantile(mutation$phylop[,6], c(0.9))), as.numeric(quantile(mutation$phylop[,6], c(1))))$burden
    gerp_gt_95pct = burden_with_gerp(mut, mutation$gerp, as.numeric(quantile(mutation$gerp[,8], c(0.95))), as.numeric(quantile(mutation$gerp[,8], c(1))))$burden
    gerp_gt_90pct = burden_with_gerp(mut, mutation$gerp, as.numeric(quantile(mutation$gerp[,8], c(0.9))), as.numeric(quantile(mutation$gerp[,8], c(1))))$burden
    ASD_stringent_gerp_95pct = burden_with_gerp(burden_with_geneset(mut,stringent_ASD)$mut, mutation$gerp, as.numeric(quantile(mutation$gerp[,8], c(0.95))), as.numeric(quantile(mutation$gerp[,8], c(1))))$burden
    ASD_stringent_gerp_90pct = burden_with_gerp(burden_with_geneset(mut,stringent_ASD)$mut, mutation$gerp, as.numeric(quantile(mutation$gerp[,8], c(0.9))), as.numeric(quantile(mutation$gerp[,8], c(1))))$burden
    ASD_relaxed_gerp_95pct = burden_with_gerp(burden_with_geneset(mut,relaxed_ASD)$mut, mutation$gerp, as.numeric(quantile(mutation$gerp[,8], c(0.95))), as.numeric(quantile(mutation$gerp[,8], c(1))))$burden
    ASD_relaxed_gerp_90pct = burden_with_gerp(burden_with_geneset(mut,relaxed_ASD)$mut, mutation$gerp, as.numeric(quantile(mutation$gerp[,8], c(0.9))), as.numeric(quantile(mutation$gerp[,8], c(1))))$burden
    nonASD_gerp_95pct = burden_with_gerp(burden_with_geneset(mut,nonASD_genes)$mut, mutation$gerp, as.numeric(quantile(mutation$gerp[,8], c(0.95))), as.numeric(quantile(mutation$gerp[,8], c(1))))$burden
    nonASD_gerp_90pct = burden_with_gerp(burden_with_geneset(mut,nonASD_genes)$mut, mutation$gerp, as.numeric(quantile(mutation$gerp[,8], c(0.9))), as.numeric(quantile(mutation$gerp[,8], c(1))))$burden
    constraint_gerp_95pct = burden_with_gerp(burden_with_geneset(mut,constraint_union)$mut, mutation$gerp, as.numeric(quantile(mutation$gerp[,8], c(0.95))), as.numeric(quantile(mutation$gerp[,8], c(1))))$burden
    constraint_gerp_90pct = burden_with_gerp(burden_with_geneset(mut,constraint_union)$mut, mutation$gerp, as.numeric(quantile(mutation$gerp[,8], c(0.9))), as.numeric(quantile(mutation$gerp[,8], c(1))))$burden
    
    output = data.frame(rbind(baseline,
                              Tada_q0.1,
                              Tada_q0.1_motif,
                              Tada_q0.2,
                              Tada_q0.2_motif,
                              Tada_q0.3,
                              Tada_q0.3_motif,
                              Tada_q0.4,
                              Tada_q0.4_motif,
                              Tada_q0.5,
                              Tada_q0.5_motif,
                              motif,
                              ASD_stringent,
                              ASD_stringent_motif,
                              ASD_relaxed,
                              ASD_relaxed_motif,
                              nonASD,
                              nonASD_motif,
                              constraint_union_data1,
                              constraint_union_data1_motif,
                              CADD,
                              ASD_stringent_CADD,
                              ASD_relaxed_CADD,
                              nonASD_CADD,
                              constraint_CADD,
                              CADD_90pct,
                              ASD_stringent_CADD_90pct,
                              ASD_relaxed_CADD_90pct,
                              nonASD_CADD_90pt,
                              constraint_CADD_90pct,
                              eigen_gt_95pct,
                              ASD_stringent_eigen_95pct,
                              ASD_relaxed_eigen_95pct,
                              nonASD_eigen_95pct,
                              constraint_eigen_95pct,
                              eigen_gt_90pct,
                              ASD_stringent_eigen_90pct,
                              ASD_relaxed_eigen_90pct,
                              nonASD_eigen_95pct,
                              constraint_eigen_90pct,
                              phylop_gt_95pct,
                              ASD_stringent_phylop_95pct,
                              ASD_relaxed_phylop_95pct,
                              nonASD_phylop_95pct,
                              constraint_phylop_95pct,
                              phylop_gt_90pct,
                              ASD_stringent_phylop_90pct,
                              ASD_relaxed_phylop_90pct,
                              nonASD_phylop_90pct,
                              constraint_phylop_90pct,
                              gerp_gt_95pct,
                              ASD_stringent_gerp_95pct,
                              ASD_relaxed_gerp_95pct,
                              nonASD_gerp_95pct,
                              constraint_gerp_95pct,
                              gerp_gt_90pct,
                              ASD_stringent_gerp_90pct,
                              ASD_relaxed_gerp_90pct,
                              nonASD_gerp_90pct,
                              constraint_gerp_90pct
                              #constraint_union_data2,
                              #constraint_union_data2_motif,
                              #ASD_predicted,
                              #ASD_predicted_motif,
                              #ASD_predicted_sensitive,
                              #ASD_predicted_motif_sensitive
                              ))
    rownames(output) = c("baseline",
                         "Tada_q0.1",
                         "Tada_q0.1_motif",
                         "Tada_q0.2",
                         "Tada_q0.2_motif",
                         "Tada_q0.3",
                         "Tada_q0.3_motif",
                         "Tada_q0.4",
                         "Tada_q0.4_motif",
                         "Tada_q0.5",
                         "Tada_q0.5_motif",
                         "motif",
                         "stringent ASD genes",
                         "stringent ASD genes + motif",
                         "relaxed ASD genes",
                         "relaxed ASD genes + motif",
                         "nonASD",
                         "nonASD+motif",
                         "constraint/haploinsufficient genes",
                         "constraint/haploinsufficient genes + motif",
                         "CADD",
                         "stringent ASD + CADD",
                         "relaxed ASD + CADD",
                         "nonASD + CADD",
                         "constraint + CADD",
                         "CADD_90pct",
                         "stringent ASD + CADD_90pct",
                         "relaxed ASD + CADD+90pct",
                         "nonASD + CADD_90pct",
                         "constraint + CADD+90pct",
                         "eigen_gt_95pct",
                         "stringent ASD + eigen_95pct",
                         "relaxed ASD + eigen_95pct",
                         "nonASD + eigen_95pct",
                         "constraint + eigen_95pct",
                         "eigen_gt_90pct",
                         "stringent ASD + eigen_90pct",
                         "relaxed ASD + eigen_90pct",
                         "nonASD + eigen_90pct",
                         "constraint + eigen_90pct",
                         "phylop_gt_95pct",
                         "stringent ASD + phylop_95pct",
                         "relaxed ASD + phylop_95pct",
                         "nonASD + phylop_95pct",
                         "constraint + phylop_95pct",
                         "phylop_gt_90pct",
                         "stringent ASD + phylop_90pct",
                         "relaxed ASD + phylop_90pct",
                         "nonASD + phylop_90pct",
                         "constraint + phylop_90pct",
                         "gerp_gt_95pct",
                         "stringent ASD + gerp_95pct",
                         "relaxed ASD + gerp_95pct",
                         "nonASD + gerp_95pct",
                         "constraint + gerp_95pct",
                         "gerp_gt_90pct",
                         "stringent ASD + gerp_90pct",
                         "relaxed ASD + gerp_90pct",
                         "nonASD + gerp_90pct",
                         "constraint + gerp_90pct"
                         #"constraint_union_data2",
                         #"constraint_union_data2_motif",
                         #"ASD_precited",
                         #"ASD_predicted_motif",
                         #"ASD_predicted_sensitive",
                         #"ASD_predicted_motif_sensitive"
                         )
    colnames(output) = c("ASD", "control", "frequency_ratio_to_background")
    if( ref == "self"){
      sim_p = t(mapply(get_pvalue_w.burden4,output[1,1],output[1,2],output[,1],output[,2]))
      odds_ratio = (output[,1]/(output[1,1]-output[,1]))/(output[,2]/(output[1,2]-output[,2]))
      burden = (output[,1]/output[,2])/(output[1,1]/output[1,2])
      output = data.frame(output, burden_to_baseline = burden, odds_ratio =odds_ratio, pvalue = as.numeric(sim_p[,1]), lowerbound = as.numeric(sim_p[,2]), upperbound = as.numeric(sim_p[,3]) )
    }
    else if (ref == "total_ratio"){
      sim_p = t(mapply(get_pvalue_w.burden4,syn_A,syn_C,output[,1],output[,2]))
      odds_ratio = (output[,1]/(syn_A-output[,1]))/(output[,2]/(syn_C-output[,2]))
      burden = (output[,1]/syn_A)/(output[,2]/syn_C)
      output = data.frame(output, burden_to_baseline = burden, odds_ratio =odds_ratio, pvalue = as.numeric(sim_p[,1]), lowerbound = as.numeric(sim_p[,2]), upperbound = as.numeric(sim_p[,3]) )
      
    }
  }
  else if (mut_type == "coding"){
    output = show_burden_for_coding_selected_features(mut, mutation$exon_mut_type)
  }
  output 
}



draw_burden_for_selected_features <-function(mut_file_name, input_type = "file", mut_type = "noncoding",ref = "self", syn_A = 0, syn_C = 0){
  #input_type could be file or a table matching regions with mutations. 
  # mut_type could be "noncoding" or "coding
  # this function generates burden table and plot included categories based on discussion on 160720
  # ref = "self" means will use the first row of the output table as baseline, and use negative binomial to calculate pvalues
  # ref = "synon_ratio" means using relative frequency ratio of synonymous mutations as baseline, and use fisher.test to 
  # get pvalues
  # syn_A is the number of synonymous mutations in ASD
  # syn_C is the number of synonymous mutations in control
  if(input_type == "file"){
    mut = read.delim(mut_file_name, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  }
  else if (input_type == "table"){
    mut = mut_file_name
  }
  if(mut_type == "noncoding"){
    baseline = burden(mut)
    ASD_stringent = burden_with_geneset(mut, stringent_ASD)$burden
    ASD_relaxed = burden_with_geneset(mut, relaxed_ASD)$burden
    nonASD = burden_with_geneset(mut, coding_tada[coding_tada$qvalue.combined > quantile(coding_tada$qvalue.combined, seq(0,1,0.1))[10],]$TadaName)$burden #bottom 10%
    constraint_union_data1 = burden_with_geneset(mut, constraint_union)$burden
    nonconstrained_burden = burden_with_geneset(mut, nonconstraint_union)$burden
    motif_burden = burden_with_motif(mut, mutation$jaspar_2014_motif_q0.1_dummy)$burden
    gerp_burden = burden_with_gerp(mut, mutation$gerp, 2, 6.17)$burden
    Phylop_burden = burden_with_phylop(mut, mutation$phylop, 1.49205, 9.87300)$burden # phylop 5%
    Eigen_burden = eigen_gt_95pct = burden_with_eigen(mut, mutation$mut_eigen, 0.82931625)$burden #eigen top 5%
    CADD_burden = burden_with_CADD(mut, CADD_score_with_index, 13.71,53)$burden #CADD top 5%
    output = data.frame(rbind(baseline,
                              ASD_stringent,
                              ASD_relaxed,
                              nonASD,
                              constraint_union_data1,
                              nonconstrained_burden,
                              motif_burden,
                              gerp_burden,
                              Phylop_burden, 
                              Eigen_burden,
                              CADD_burden
    ))
    rownames(output) = c("baseline",
                         "ASD genes",
                         "Neuropsychiatric genes",
                         "nonASD genes",
                         "Intolerant genes",
                         "Tolerant genes",
                         "Motif",
                         "Gerp >=2",
                         "Phylop top5%",
                         "Eigen top5%",
                         "CADD top5%"
                         
    )
    colnames(output) = c("ASD", "control", "frequency_ratio_to_background")
    if( ref == "self"){
      sim_p = t(mapply(get_pvalue_w.burden4,output[1,1],output[1,2],output[,1],output[,2]))
      odds_ratio = (output[,1]/(output[1,1]-output[,1]))/(output[,2]/(output[1,2]-output[,2]))
      burden = (output[,1]/output[,2])/(output[1,1]/output[1,2])
      output = data.frame(output, burden_to_baseline = burden, odds_ratio =odds_ratio, pvalue = as.numeric(sim_p[,1]), lowerbound = as.numeric(sim_p[,2]), upperbound = as.numeric(sim_p[,3]) )
    }
    else if (ref == "synon_ratio"){
      sim_p = t(mapply(get_pvalue_w.burden3,syn_A,syn_C,output[,1],output[,2]))
      odds_ratio = (output[,1]/syn_A)/(output[,2]/syn_C)
      burden = (output[,1]/syn_A)/(output[,2]/syn_C)
      output = data.frame(output, burden_to_baseline = burden, odds_ratio =odds_ratio, pvalue = as.numeric(sim_p[,1]), lowerbound = as.numeric(sim_p[,2]), upperbound = as.numeric(sim_p[,3]) )
      
    }
  }
  else if (mut_type == "coding"){
    output = show_burden_for_coding_selected_features(mut, mutation$exon_mut_type)
  }
  output 
}


show_mut_for_Table1 <-function(mut_file_name, input_type = "file", mut_type = "noncoding"){
  #input_type could be file or a table matching regions with mutations. 
  # mut_type could be "noncoding" or "coding
  if(input_type == "file"){
    mut = read.delim(mut_file_name, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  }
  else if (input_type == "table"){
    mut = mut_file_name
  }
  if(mut_type == "noncoding"){
    write.table(intersect(mut[,1],union(mutation$ASD_effective_SNV_ID,mutation$control_effective_SNV_ID)),"160616_all_genes.txt", col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)
    Tada_q0.3 = intersect(burden_with_geneset(mut, coding_tada[coding_tada$qvalue.combined<0.3,2])$mut[,1],union(mutation$ASD_effective_SNV_ID,mutation$control_effective_SNV_ID))
    write.table(Tada_q0.3,"160616_Tada_q0.3.txt", col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)
    Tada_q0.3_motif = intersect(burden_with_motif(burden_with_geneset(mut,coding_tada[coding_tada$qvalue.combined<0.3,2])$mut,mutation$jaspar_2014_motif_q0.1_dummy)$mut[,1],union(mutation$ASD_effective_SNV_ID,mutation$control_effective_SNV_ID))
    write.table(Tada_q0.3_motif,"160616_Tada_q0.3_motif.txt", col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)
    motif = intersect(burden_with_motif(mut, mutation$jaspar_2014_motif_q0.1_dummy)$mut[,1],union(mutation$ASD_effective_SNV_ID,mutation$control_effective_SNV_ID))
    write.table(motif,"160616_motif.txt", col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)
    ASD_stringent = intersect(burden_with_geneset(mut, stringent_ASD)$mut[,1],union(mutation$ASD_effective_SNV_ID,mutation$control_effective_SNV_ID))
    write.table(ASD_stringent,"160616_Known_ASD.txt", col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)
    ASD_stringent_motif = intersect(burden_with_motif(burden_with_geneset(mut,stringent_ASD)$mut,mutation$jaspar_2014_motif_q0.1_dummy)$mut[,1],union(mutation$ASD_effective_SNV_ID,mutation$control_effective_SNV_ID))
    write.table(ASD_stringent_motif,"160616_Known_ASD_motif.txt", col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)
    ASD_relaxed = intersect(burden_with_geneset(mut, relaxed_ASD)$mut[,1],union(mutation$ASD_effective_SNV_ID,mutation$control_effective_SNV_ID))
    write.table(ASD_relaxed,"160616_putative_ASD.txt", col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)
    ASD_relaxed_motif = intersect(burden_with_motif(burden_with_geneset(mut,relaxed_ASD)$mut,mutation$jaspar_2014_motif_q0.1_dummy)$mut[,1],union(mutation$ASD_effective_SNV_ID,mutation$control_effective_SNV_ID))
    write.table(ASD_relaxed_motif,"160616_putative_ASD_motif.txt", col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)
    constraint_union_data1 = intersect(burden_with_geneset(mut, constraint_union)$mut[,1],union(mutation$ASD_effective_SNV_ID,mutation$control_effective_SNV_ID))
    write.table(constraint_union_data1,"160616_constraint_haplo.txt", col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)
    constraint_union_data1_motif = intersect(burden_with_motif(burden_with_geneset(mut, constraint_union)$mut,mutation$jaspar_2014_motif_q0.1_dummy)$mut[,1],union(mutation$ASD_effective_SNV_ID,mutation$control_effective_SNV_ID))
    write.table(constraint_union_data1_motif,"160616_constraint_haplo_motif.txt", col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  else if (mut_type == "coding"){
    output = show_burden_for_coding_selected_features(mut, mutation$exon_mut_type)
  }
}



show_burden_for_SPIDEX <-function(lowerbound, upperbound){
  a = burden_with_spidex2(SPIDEX, -Inf, Inf)
  b = burden_with_spidex2(SPIDEX, lowerbound,upperbound)
  c = burden_with_spidex2_geneset(SPIDEX, lowerbound,upperbound,stringent_ASD)
  d = burden_with_spidex2_geneset(SPIDEX, lowerbound,upperbound,relaxed_ASD)
  e = burden_with_spidex2_geneset(SPIDEX, lowerbound,upperbound,constraint_union)
  f = burden_with_spidex2_geneset(SPIDEX, lowerbound,upperbound,c(brainspan$gene_qt1$gene, brainspan$gene_qt2$gene))
  g = burden_with_spidex2_geneset(SPIDEX, lowerbound,upperbound,c(brainspan$gene_qt3$gene, brainspan$gene_qt4$gene))
  h = burden_with_spidex2_geneset(SPIDEX, lowerbound,upperbound, Olfactory_gene[,1])
  i = burden_with_spidex2_geneset(SPIDEX, lowerbound,upperbound, DDG2P[,1])
  j = burden_with_spidex2_geneset(SPIDEX, lowerbound, upperbound, DDG2P_confirmed_DD[,1])
  k = burden_with_spidex2_geneset(SPIDEX, lowerbound, upperbound, coding_tada[coding_tada$qvalue.combined<0.1,2])
  l = burden_with_spidex2_geneset(SPIDEX, lowerbound, upperbound, coding_tada[coding_tada$qvalue.combined<0.2,2])
  m = burden_with_spidex2_geneset(SPIDEX, lowerbound, upperbound, coding_tada[coding_tada$qvalue.combined<0.3,2])
  n = burden_with_spidex2_geneset(SPIDEX, lowerbound, upperbound, coding_tada[coding_tada$qvalue.combined<0.4,2])
  o = burden_with_spidex2_geneset(SPIDEX, lowerbound, upperbound, coding_tada[coding_tada$qvalue.combined<0.5,2])
  p = burden_with_spidex2_geneset(SPIDEX, lowerbound,upperbound, nonconstraint_union)
  q = burden_with_spidex2_geneset(SPIDEX, lowerbound,upperbound, coding_tada[coding_tada$qvalue.combined > quantile(coding_tada$qvalue.combined, seq(0,1,0.1))[10],]$TadaName) # bottom 10%
  output = rbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q)
  output = data.frame(ASD = output[,1],control=output[,2],burden=output[,3]/output[1,3], hyper.p = phyper(output[,1]-1,output[1,1],output[1,2], output[,1]+output[,2],lower.tail = FALSE),
                      burden_compared_syn = output[,3]/0.74, burden_p_compare_syn = NA)
  for(i in 1:length(output[,1])){
    output[i,]$burden_p_compare_syn = fisher.test(matrix(c(72,output[i,1],139,output[i,2]),2,2),alternative = "less")$p.value
  }
  rownames(output) = c("all SNVs with a SPIDEX score", "threshold", "Known ASD genes", "Neuron Psychiatric genes", "Intolerant genes","Bottome 50% brain-expressed genes",
                       "Top 50% brain-expressed genes","olfactory genes","DDG2P_gene","DDG2P_gene_confirmed_DD","Tada_q0.1","Tada_q0.2","Tada_q0.3","Tada_q0.4","Tada_q0.5","Tolerant genes","nonASD genes")
  output 
}


# ################## adjusted for CG content near mutations ###################################################
# cg_with_effective_SNV_ID = mutation$cg_pct[is.element(mutation$cg_pct[,1],c(mutation$ASD_effective_SNV_ID, mutation$control_effective_SNV_ID)),]
# #using 3 bins, get indexes in each bin
# cg_with_effective_SNV_ID_bin3_1_index = cg_with_effective_SNV_ID[cg_with_effective_SNV_ID[,2] >= quantile(cg_with_effective_SNV_ID[,2],seq(0,1,1/3))[[1]] & 
#                                                                    cg_with_effective_SNV_ID[,2] < quantile(cg_with_effective_SNV_ID[,2],seq(0,1,1/3))[[2]],1]
# cg_with_effective_SNV_ID_bin3_2_index = cg_with_effective_SNV_ID[cg_with_effective_SNV_ID[,2] >= quantile(cg_with_effective_SNV_ID[,2],seq(0,1,1/3))[[2]] & 
#                                                                    cg_with_effective_SNV_ID[,2] < quantile(cg_with_effective_SNV_ID[,2],seq(0,1,1/3))[[3]],1]
# cg_with_effective_SNV_ID_bin3_3_index = cg_with_effective_SNV_ID[cg_with_effective_SNV_ID[,2] >= quantile(cg_with_effective_SNV_ID[,2],seq(0,1,1/3))[[3]] & 
#                                                                    cg_with_effective_SNV_ID[,2] <= quantile(cg_with_effective_SNV_ID[,2],seq(0,1,1/3))[[4]],1]
# 
# cg_with_effective_SNV_ID_bin3_1_scaling = (length(intersect(cg_with_effective_SNV_ID_bin3_1_index,mutation$ASD_effective_SNV_ID))/length(mutation$ASD_effective_SNV_ID))/
#   (length(intersect(cg_with_effective_SNV_ID_bin3_1_index,mutation$control_effective_SNV_ID))/length(mutation$control_effective_SNV_ID))
# cg_with_effective_SNV_ID_bin3_2_scaling = (length(intersect(cg_with_effective_SNV_ID_bin3_2_index,mutation$ASD_effective_SNV_ID))/length(mutation$ASD_effective_SNV_ID))/
#   (length(intersect(cg_with_effective_SNV_ID_bin3_2_index,mutation$control_effective_SNV_ID))/length(mutation$control_effective_SNV_ID))
# cg_with_effective_SNV_ID_bin3_3_scaling = (length(intersect(cg_with_effective_SNV_ID_bin3_3_index,mutation$ASD_effective_SNV_ID))/length(mutation$ASD_effective_SNV_ID))/
#   (length(intersect(cg_with_effective_SNV_ID_bin3_3_index,mutation$control_effective_SNV_ID))/length(mutation$control_effective_SNV_ID))


# ################## adjusted for CG content near mutations for 4 bins ###################################################
# 
# show_burden_CG_bin3 <- function(mut_file_name){
#   mut = read.delim(mut_file_name, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#   mut1 = mut[is.element(mut[,1], cg_with_effective_SNV_ID_bin3_1_index),]
#   mut2 = mut[is.element(mut[,1], cg_with_effective_SNV_ID_bin3_2_index),]
#   mut3 = mut[is.element(mut[,1], cg_with_effective_SNV_ID_bin3_3_index),]
#   a = show_burden(mut1, input_type = "table")
#   a$burden = a$burden/cg_with_effective_SNV_ID_bin3_1_scaling
#   b = show_burden(mut2, input_type = "table")
#   b$burden = b$burden/cg_with_effective_SNV_ID_bin3_2_scaling
#   c = show_burden(mut3, input_type = "table")
#   c$burden = c$burden/cg_with_effective_SNV_ID_bin3_3_scaling
#   d = cbind(a,b,c)
#   d
# }
# 
# 
# cg_with_effective_SNV_ID = mutation$cg_pct[is.element(mutation$cg_pct[,1],c(mutation$ASD_effective_SNV_ID, mutation$control_effective_SNV_ID)),]
# #using 3 bins, get indexes in each bin
# cg_with_effective_SNV_ID_bin4_1_index = cg_with_effective_SNV_ID[cg_with_effective_SNV_ID[,2] >= quantile(cg_with_effective_SNV_ID[,2],seq(0,1,1/4))[[1]] & 
#                                                                    cg_with_effective_SNV_ID[,2] < quantile(cg_with_effective_SNV_ID[,2],seq(0,1,1/4))[[2]],1]
# cg_with_effective_SNV_ID_bin4_2_index = cg_with_effective_SNV_ID[cg_with_effective_SNV_ID[,2] >= quantile(cg_with_effective_SNV_ID[,2],seq(0,1,1/4))[[2]] & 
#                                                                    cg_with_effective_SNV_ID[,2] < quantile(cg_with_effective_SNV_ID[,2],seq(0,1,1/4))[[3]],1]
# cg_with_effective_SNV_ID_bin4_3_index = cg_with_effective_SNV_ID[cg_with_effective_SNV_ID[,2] >= quantile(cg_with_effective_SNV_ID[,2],seq(0,1,1/4))[[3]] & 
#                                                                    cg_with_effective_SNV_ID[,2] < quantile(cg_with_effective_SNV_ID[,2],seq(0,1,1/4))[[4]],1]
# cg_with_effective_SNV_ID_bin4_4_index = cg_with_effective_SNV_ID[cg_with_effective_SNV_ID[,2] >= quantile(cg_with_effective_SNV_ID[,2],seq(0,1,1/4))[[4]] & 
#                                                                    cg_with_effective_SNV_ID[,2] <= quantile(cg_with_effective_SNV_ID[,2],seq(0,1,1/4))[[5]],1]
# 
# 
# cg_with_effective_SNV_ID_bin4_1_scaling = (length(intersect(cg_with_effective_SNV_ID_bin4_1_index,mutation$ASD_effective_SNV_ID))/length(mutation$ASD_effective_SNV_ID))/
#   (length(intersect(cg_with_effective_SNV_ID_bin4_1_index,mutation$control_effective_SNV_ID))/length(mutation$control_effective_SNV_ID))
# cg_with_effective_SNV_ID_bin4_2_scaling = (length(intersect(cg_with_effective_SNV_ID_bin4_2_index,mutation$ASD_effective_SNV_ID))/length(mutation$ASD_effective_SNV_ID))/
#   (length(intersect(cg_with_effective_SNV_ID_bin4_2_index,mutation$control_effective_SNV_ID))/length(mutation$control_effective_SNV_ID))
# cg_with_effective_SNV_ID_bin4_3_scaling = (length(intersect(cg_with_effective_SNV_ID_bin4_3_index,mutation$ASD_effective_SNV_ID))/length(mutation$ASD_effective_SNV_ID))/
#   (length(intersect(cg_with_effective_SNV_ID_bin4_3_index,mutation$control_effective_SNV_ID))/length(mutation$control_effective_SNV_ID))
# cg_with_effective_SNV_ID_bin4_4_scaling = (length(intersect(cg_with_effective_SNV_ID_bin4_4_index,mutation$ASD_effective_SNV_ID))/length(mutation$ASD_effective_SNV_ID))/
#   (length(intersect(cg_with_effective_SNV_ID_bin4_4_index,mutation$control_effective_SNV_ID))/length(mutation$control_effective_SNV_ID))
# 

#####################function to simulate pvalue for the normalized burdne ######################################################################
get_pvalue_w.burden <- function(total_obs, mean1, mean2, mean3, mean4, sim_rounds = 10000){
  data = matrix(rpois(4*sim_rounds,c(mean1,mean2,mean3,mean4)),4,sim_rounds)
  data_mean = apply(data,2,sum)
  p = length(data_mean[data_mean >= total_obs])/length(data_mean)
  p
}

get_pvalue_w.burden2 <- function(ASD_count_baseline, mutation_count_baseline,ASD_count,mutation_count){# 
  phyper(ASD_count-1,ASD_count_baseline,mutation_count_baseline,(ASD_count + mutation_count),lower.tail = FALSE)
}

get_pvalue_w.burden3 <- function(syn_A, syn_B,ASD_count,control_count){# 
  a = fisher.test(matrix(c(syn_B,control_count,syn_A, ASD_count),2,2), alternative = "greater")
  b = fisher.test(matrix(c(syn_B,control_count,syn_A, ASD_count),2,2))
  list(pvalue = a$p.value, lowerbound = b$conf.int[1], upperbound = b$conf.int[2])
}

# use Fisher one tailed test to estimate confidence interval for odds ratio
#use two-tailed test to get confidence interval
get_pvalue_w.burden4 <- function(ASD_count_baseline, control_count_baseline,ASD_count,control_count){# 
  a = fisher.test(matrix(c(control_count_baseline-control_count, control_count, ASD_count_baseline-ASD_count, ASD_count),2,2), alternative = "greater")
  b = fisher.test(matrix(c(control_count_baseline-control_count, control_count, ASD_count_baseline-ASD_count, ASD_count),2,2))
  list(pvalue = a$p.value, lowerbound = b$conf.int[1], upperbound = b$conf.int[2])
}


# 1010  write a function to translate the output of show_burden_for_selected_features into a data.frame that is going to be used for plotting using Yanyu's code
# the output for mat would be similar to 
# gene_list score cutoff burden CI_low CI_up pvalue asd_sub control_sub bed region_tag
# 1     total  GERP   0.15  1.967  1.168   Inf  0.014      34          17 DHS      3'UTR
# 2   NP_gene  GERP   0.15  2.212  0.734   Inf   0.14       9           4 DHS      3'UTR
# 3     total  GERP      1  1.352  0.751   Inf  0.224      22          16 DHS      3'UTR
# 4   NP_gene  GERP      1  1.721  0.529   Inf  0.284       7           4 DHS      3'UTR
# 5     total  GERP      2  1.202  0.522   Inf  0.427      11           9 DHS      3'UTR
# 6   NP_gene  GERP      2  3.933  0.512   Inf  0.193       4           1 DHS      3'UTR

translate_to_plotting <- function(data,epi_bed, region_tag){
  # data is the output data.frame of show_burden_for_selected_features
  # epi_bed is the epigenomic partition that is used, including "3'UTR", "5'UTR", "Enhancer 0-10kb", "Enhancer 10-20kb", "Promoter"
  # region_tag is the region partition, including "DHS", "Fantom", "Noonan", "Noonan_Roadmap_Intersect", "Noonan_Roadmap_Union", "Roadmap", "Whole Genome"
  temp_table = data.frame(gene_list = NA, score = NA, burden = NA, CI_low = NA, CI_up = NA, pvalue = NA, asd_sub = NA, control_sub = NA, bed = NA, region_tag = NA)
  temp_table[1,] = c("total", "GERP", as.numeric(quantile(mutation$gerp[,8], c(0.95))), data["gerp_gt_95pct",4], data["gerp_gt_95pct",7],data["gerp_gt_95pct",8], data["gerp_gt_95pct",6], data["gerp_gt_95pct",1], data["gerp_gt_95pct",2], epi_bed, region_tag)
  temp_table[2,] = c("total", "GERP", as.numeric(quantile(mutation$gerp[,8], c(0.9))), data["gerp_gt_90pct",4], data["gerp_gt_90pct",7],data["gerp_gt_90pct",8], data["gerp_gt_90pct",6], data["gerp_gt_90pct",1], data["gerp_gt_90pct",2], epi_bed, region_tag)
  temp_table[3,] = c("total", "CADD13_PHRED", as.numeric(quantile(mutation$CADD_score[,1], c(0.95))), data["CADD",4], data["CADD",7],data["CADD",8], data["CADD",6], data["CADD",1], data["CADD",2], epi_bed, region_tag)
  temp_table[4,] = c("total", "CADD13_PHRED", as.numeric(quantile(mutation$CADD_score[,1], c(0.9))), data["CADD_90pct",4], data["CADD_90pct",7],data["CADD_90pct",8], data["CADD_90pct",6], data["CADD_90pct",1], data["CADD_90pct",2], epi_bed, region_tag)
  temp_table[5,] = c("total", "PhyloP", as.numeric(quantile(mutation$phylop[,6], c(0.95))), data["phylop_gt_95pct",4], data["phylop_gt_95pct",7],data["phylop_gt_95pct",8], data["phylop_gt_95pct",6], data["phylop_gt_95pct",1], data["phylop_gt_95pct",2], epi_bed, region_tag)
  temp_table[6,] = c("total", "PhyloP", as.numeric(quantile(mutation$phylop[,6], c(0.9))), data["phylop_gt_90pct",4], data["phylop_gt_90pct",7],data["phylop_gt_90pct",8], data["phylop_gt_90pct",6], data["phylop_gt_90pct",1], data["phylop_gt_90pct",2], epi_bed, region_tag)
  temp_table[7,] = c("total", "Eigen", as.numeric(quantile(mutation$mut_eigen, c(0.95))), data["eigen_gt_95pct",4], data["eigen_gt_95pct",7],data["eigen_gt_95pct",8], data["eigen_gt_95pct",6], data["eigen_gt_95pct",1], data["eigen_gt_95pct",2], epi_bed, region_tag)
  temp_table[8,] = c("total", "Eigen", as.numeric(quantile(mutation$mut_eigen, c(0.9))), data["eigen_gt_90pct",4], data["eigen_gt_90pct",7],data["eigen_gt_90pct",8], data["eigen_gt_90pct",6], data["eigen_gt_90pct",1], data["eigen_gt_90pct",2], epi_bed, region_tag)
  temp_table[9,] = c("total", "Eigen", 0.1, data["motif",4], data["motif",7],data["motif",8], data["motif",6], data["motif",1], data["motif",2], epi_bed, region_tag)
  
  temp_table[10,] = c("NP_genes", "GERP", as.numeric(quantile(mutation$gerp[,8], c(0.95))), data["gerp_gt_95pct",4], data["gerp_gt_95pct",7],data["gerp_gt_95pct",8], data["gerp_gt_95pct",6], data["gerp_gt_95pct",1], data["gerp_gt_95pct",2], epi_bed, region_tag)
  temp_table[11,] = c("NP_genes", "GERP", as.numeric(quantile(mutation$gerp[,8], c(0.9))), data["gerp_gt_90pct",4], data["gerp_gt_90pct",7],data["gerp_gt_90pct",8], data["gerp_gt_90pct",6], data["gerp_gt_90pct",1], data["gerp_gt_90pct",2], epi_bed, region_tag)
  temp_table[12,] = c("NP_genes", "CADD13_PHRED", as.numeric(quantile(mutation$CADD_score[,1], c(0.95))), data["CADD",4], data["CADD",7],data["CADD",8], data["CADD",6], data["CADD",1], data["CADD",2], epi_bed, region_tag)
  temp_table[13,] = c("NP_genes", "CADD13_PHRED", as.numeric(quantile(mutation$CADD_score[,1], c(0.9))), data["CADD_90pct",4], data["CADD_90pct",7],data["CADD_90pct",8], data["CADD_90pct",6], data["CADD_90pct",1], data["CADD_90pct",2], epi_bed, region_tag)
  temp_table[14,] = c("NP_genes", "PhyloP", as.numeric(quantile(mutation$phylop[,6], c(0.95))), data["phylop_gt_95pct",4], data["phylop_gt_95pct",7],data["phylop_gt_95pct",8], data["phylop_gt_95pct",6], data["phylop_gt_95pct",1], data["phylop_gt_95pct",2], epi_bed, region_tag)
  temp_table[15,] = c("NP_genes", "PhyloP", as.numeric(quantile(mutation$phylop[,6], c(0.9))), data["phylop_gt_90pct",4], data["phylop_gt_90pct",7],data["phylop_gt_90pct",8], data["phylop_gt_90pct",6], data["phylop_gt_90pct",1], data["phylop_gt_90pct",2], epi_bed, region_tag)
  temp_table[16,] = c("NP_genes", "Eigen", as.numeric(quantile(mutation$mut_eigen, c(0.95))), data["eigen_gt_95pct",4], data["eigen_gt_95pct",7],data["eigen_gt_95pct",8], data["eigen_gt_95pct",6], data["eigen_gt_95pct",1], data["eigen_gt_95pct",2], epi_bed, region_tag)
  temp_table[17,] = c("NP_genes", "Eigen", as.numeric(quantile(mutation$mut_eigen, c(0.9))), data["eigen_gt_90pct",4], data["eigen_gt_90pct",7],data["eigen_gt_90pct",8], data["eigen_gt_90pct",6], data["eigen_gt_90pct",1], data["eigen_gt_90pct",2], epi_bed, region_tag)
  temp_table[18,] = c("NP_genes", "Eigen", 0.1, data["motif",4], data["motif",7],data["motif",8], data["motif",6], data["motif",1], data["motif",2], epi_bed, region_tag)
  
  
}


# show_burden_CG_bin4 <- function(mut_file_name, mut_type = "noncoding"){ #could also be chosen to be "coding
#   mut = read.delim(mut_file_name, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#   mut1 = mut[is.element(mut[,1], cg_with_effective_SNV_ID_bin4_1_index),]
#   mut2 = mut[is.element(mut[,1], cg_with_effective_SNV_ID_bin4_2_index),]
#   mut3 = mut[is.element(mut[,1], cg_with_effective_SNV_ID_bin4_3_index),]
#   mut4 = mut[is.element(mut[,1], cg_with_effective_SNV_ID_bin4_4_index),]
#   a = show_burden(mut1, input_type = "table", mut_type)
#   a$burden = a$burden/cg_with_effective_SNV_ID_bin4_1_scaling
#   b = show_burden(mut2, input_type = "table", mut_type)
#   b$burden = b$burden/cg_with_effective_SNV_ID_bin4_2_scaling
#   c = show_burden(mut3, input_type = "table", mut_type)
#   c$burden = c$burden/cg_with_effective_SNV_ID_bin4_3_scaling 
#   d = show_burden(mut4, input_type = "table", mut_type)
#   d$burden = d$burden/cg_with_effective_SNV_ID_bin4_4_scaling 
#   e = show_burden(mut, input_type = "table", mut_type)
#   f = cbind(a,b,c,d,e)
#   mean1 = f[,2]*cg_with_effective_SNV_ID_bin4_1_scaling*length(mutation$ASD_effective_SNV_ID)/length(mutation$control_effective_SNV_ID)
#   mean2 = f[,5]*cg_with_effective_SNV_ID_bin4_2_scaling*length(mutation$ASD_effective_SNV_ID)/length(mutation$control_effective_SNV_ID)
#   mean3 = f[,8]*cg_with_effective_SNV_ID_bin4_3_scaling*length(mutation$ASD_effective_SNV_ID)/length(mutation$control_effective_SNV_ID)
#   mean4 = f[,11]*cg_with_effective_SNV_ID_bin4_4_scaling*length(mutation$ASD_effective_SNV_ID)/length(mutation$control_effective_SNV_ID)
#   w.burden = (f[,1]+f[,4]+f[,7]+f[,10])/(mean1+mean2+mean3+mean4)
#   total_obs = mapply(sum,f[,1],f[,4],f[,7],f[,10])
#   #calculate the new simp for 
#   temp = data.frame(total_obs = total_obs, mean1 = mean1, mean2 = mean2, mean3 = mean3, mean4 = mean4)
#   sim_p = mapply(get_pvalue_w.burden2,e[1,1],e[1,2],e[,1],e[,2])
#   f = data.frame(f,w.burden,sim_p,rownames(f))
#   f
# }


# #####################function to draw bar plot of weighted burden fro 4 bins ############################################
# #under 10kb, 25kb, 50kb, 100kb, 500kb and 1000kb cutoff
# library(ggplot2)
# draw_wtd_burden <-function(file, mode = "genebody"){ #or mode = distance to TSS
#   if(mode == "distance to TSS"){
#     a = data.frame(w.burden =show_burden_CG_bin4(paste(file,"10000.bp_within_TSS_overlap_mutation.txt",sep = "."))$w.burden,
#                  label = rownames(show_burden_CG_bin4(paste(file,"10000.bp_within_TSS_overlap_mutation.txt",sep = "."))),
#                  distance = "10kb")
#     b = data.frame(w.burden =show_burden_CG_bin4(paste(file,"25000.bp_within_TSS_overlap_mutation.txt",sep = "."))$w.burden,
#                  label = rownames(show_burden_CG_bin4(paste(file,"25000.bp_within_TSS_overlap_mutation.txt",sep = "."))),
#                  distance = "25kb")
#     c = data.frame(w.burden =show_burden_CG_bin4(paste(file,"50000.bp_within_TSS_overlap_mutation.txt",sep = "."))$w.burden,
#                  label = rownames(show_burden_CG_bin4(paste(file,"50000.bp_within_TSS_overlap_mutation.txt",sep = "."))),
#                  distance = "50kb")
#     d = data.frame(w.burden =show_burden_CG_bin4(paste(file,"100000.bp_within_TSS_overlap_mutation.txt",sep = "."))$w.burden,
#                  label = rownames(show_burden_CG_bin4(paste(file,"100000.bp_within_TSS_overlap_mutation.txt",sep = "."))),
#                  distance = "100kb") 
#     e = data.frame(w.burden =show_burden_CG_bin4(paste(file,"500000.bp_within_TSS_overlap_mutation.txt",sep = "."))$w.burden,
#                  label = rownames(show_burden_CG_bin4(paste(file,"500000.bp_within_TSS_overlap_mutation.txt",sep = "."))),
#                  distance = "500kb")
#     f = data.frame(w.burden =show_burden_CG_bin4(paste(file,"1000000.bp_within_TSS_overlap_mutation.txt",sep = "."))$w.burden,
#                  label = rownames(show_burden_CG_bin4(paste(file,"1000000.bp_within_TSS_overlap_mutation.txt",sep = "."))),
#                  distance = "1000kb")
#   }
#   else if(mode == "genebody with intron"){
#     a = data.frame(w.burden =show_burden_CG_bin4(paste(file,"10kb_with_intron_overlap_mutation.txt",sep = "_"))$w.burden,
#                    label = rownames(show_burden_CG_bin4(paste(file,"10kb_with_intron_overlap_mutation.txt",sep = "_"))),
#                    distance = "10kb")
#     b = data.frame(w.burden =show_burden_CG_bin4(paste(file,"25kb_with_intron_overlap_mutation.txt",sep = "_"))$w.burden,
#                    label = rownames(show_burden_CG_bin4(paste(file,"25kb_with_intron_overlap_mutation.txt",sep = "_"))),
#                    distance = "25kb")
#     c = data.frame(w.burden =show_burden_CG_bin4(paste(file,"50kb_with_intron_overlap_mutation.txt",sep = "_"))$w.burden,
#                    label = rownames(show_burden_CG_bin4(paste(file,"50kb_with_intron_overlap_mutation.txt",sep = "_"))),
#                    distance = "50kb")
#     d = data.frame(w.burden =show_burden_CG_bin4(paste(file,"100kb_with_intron_overlap_mutation.txt",sep = "_"))$w.burden,
#                    label = rownames(show_burden_CG_bin4(paste(file,"100kb_with_intron_overlap_mutation.txt",sep = "_"))),
#                    distance = "100kb") 
#     e = data.frame(w.burden =show_burden_CG_bin4(paste(file,"500kb_with_intron_overlap_mutation.txt",sep = "_"))$w.burden,
#                    label = rownames(show_burden_CG_bin4(paste(file,"500kb_with_intron_overlap_mutation.txt",sep = "_"))),
#                    distance = "500kb")
#     f = data.frame(w.burden =show_burden_CG_bin4(paste(file,"1000kb_with_intron_overlap_mutation.txt",sep = "_"))$w.burden,
#                    label = rownames(show_burden_CG_bin4(paste(file,"1000kb_with_intron_overlap_mutation.txt",sep = "_"))),
#                    distance = "1000kb") 
#   }
#   else if(mode == "genebody without intron"){
#     a = data.frame(w.burden =show_burden_CG_bin4(paste(file,"10kb_without_intron_overlap_mutation.txt",sep = "_"))$w.burden,
#                    label = rownames(show_burden_CG_bin4(paste(file,"10kb_without_intron_overlap_mutation.txt",sep = "_"))),
#                    distance = "10kb")
#     b = data.frame(w.burden =show_burden_CG_bin4(paste(file,"25kb_without_intron_overlap_mutation.txt",sep = "_"))$w.burden,
#                    label = rownames(show_burden_CG_bin4(paste(file,"25kb_without_intron_overlap_mutation.txt",sep = "_"))),
#                    distance = "25kb")
#     c = data.frame(w.burden =show_burden_CG_bin4(paste(file,"50kb_without_intron_overlap_mutation.txt",sep = "_"))$w.burden,
#                    label = rownames(show_burden_CG_bin4(paste(file,"50kb_without_intron_overlap_mutation.txt",sep = "_"))),
#                    distance = "50kb")
#     d = data.frame(w.burden =show_burden_CG_bin4(paste(file,"100kb_without_intron_overlap_mutation.txt",sep = "_"))$w.burden,
#                    label = rownames(show_burden_CG_bin4(paste(file,"100kb_without_intron_overlap_mutation.txt",sep = "_"))),
#                    distance = "100kb") 
#     e = data.frame(w.burden =show_burden_CG_bin4(paste(file,"500kb_without_intron_overlap_mutation.txt",sep = "_"))$w.burden,
#                    label = rownames(show_burden_CG_bin4(paste(file,"500kb_without_intron_overlap_mutation.txt",sep = "_"))),
#                    distance = "500kb")
#     f = data.frame(w.burden =show_burden_CG_bin4(paste(file,"1000kb_without_intron_overlap_mutation.txt",sep = "_"))$w.burden,
#                    label = rownames(show_burden_CG_bin4(paste(file,"1000kb_without_intron_overlap_mutation.txt",sep = "_"))),
#                    distance = "1000kb") 
#   }              
#   all_burden = rbind(a,b,c,d,e,f)
#   all_burden$label = factor(all_burden$label, levels = a$label)
#   all_burden$distance = factor(all_burden$distance, levels = c("10kb","25kb","50kb","100kb","500kb","1000kb"))
#   p <- ggplot(all_burden, aes(x = distance, y = w.burden, fill = factor(label))) + geom_bar(position = "dodge", stat="identity")
#   dodge <- position_dodge(width=0.9)
#   p+geom_hline(yintercept=1)
#   
# }
# 
# draw_wtd_burden_for_promoter_utr5 <-function(file){
#   all_burden = data.frame(w.burden = show_burden_CG_bin4(file)$w.burden, label = rownames(show_burden_CG_bin4(file)))
#   all_burden$label = factor(all_burden$label, levels = all_burden$label)
#   p <- ggplot(all_burden, aes(x = label, y = w.burden)) + geom_bar(position = "dodge", stat="identity")
#   dodge <- position_dodge(width=0.9)
#   p+geom_hline(yintercept=1)+theme(axis.text.x = element_text(angle = 90))
# }
# 
# draw_wtd_burden_for_coding <-function(file){
#   all_burden = data.frame(w.burden = show_burden_CG_bin4(file, "coding")$w.burden, label = rownames(show_burden_CG_bin4(file, "coding")))
#   all_burden$label = factor(all_burden$label, levels = all_burden$label)
#   p <- ggplot(all_burden, aes(x = label, y = w.burden)) + geom_bar(position = "dodge", stat="identity")
#   dodge <- position_dodge(width=0.9)
#   p+geom_hline(yintercept=1)+theme(axis.text.x = element_text(angle = 90))
# }