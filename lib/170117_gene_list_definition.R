brainspan = new.env()
load("../other_annotation/brainspan/160206_brainspan_expression.Rdata", envir = brainspan)

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

ExAC_gene = as.data.frame(readRDS("../other_annotation/ConstraintMat.RDS"))
perc.rank <- function(x) trunc(rank(-x))/length(x)
ExAC_gene = data.frame(ExAC_gene, z_mis_pct = perc.rank(ExAC_gene$mis_z), z_lof_pct = perc.rank(ExAC_gene$lof_z))


Nick_ASD_gene = read.delim("../other_annotation/gene_list/Nick_machine_learning_ASD_gene_data.txt",
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)

all_refseq_genes = coding_tada$RefSeqName

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
ID_genes = Pinto_ID_genelist[,1]
SCZ_composite = Purcell_NG_2014_composite_gene[,1]

#relaxed_ASD genes are named neuronpsychiatric genes
relaxed_ASD = c(coding_tada[coding_tada$qvalue.combined<0.5,2],SFRAI_gene_high_confidence[,1],SFRAI_gene_strong_candidate[,1],Pinto_ID_genelist[,1],
                Bayes_postsynaptic_high_confidence_genelist[,1],Purcell_NG_2014_composite_gene[,1],AutismKB_gene[,1],Betancur_gene[,1],
                SFRAI_gene_syndromic[,1], SFRAI_gene_suggestive_evidence[,1])

#new definition of putative ASD 
#relaxed_ASD = c(coding_tada[coding_tada$qvalue.combined<0.3,2],SFRAI_gene_high_confidence[,1],SFRAI_gene_strong_candidate[,1],Pinto_ID_genelist[,1],
#                     AutismKB_gene[,1],Betancur_gene[,1],SFRAI_gene_syndromic[,1],SFRAI_gene_suggestive_evidence[,1], SFRAI_gene_hypothesized_gene[,1])

#nonASD_genes = coding_tada[coding_tada$qvalue.combined > quantile(coding_tada$qvalue.combined, seq(0,1,0.1))[10],]$TadaName
nonASD_genes = coding_tada[coding_tada$qvalue.combined > quantile(coding_tada$qvalue.combined, seq(0,1,0.1))[9],]$TadaName

constraint_union = c(Petrovski_RVIS_top5_pct_gene[,1],Petrovski_haploinsufficiency_genelist[,1],Huang_Hpi_score[Huang_Hpi_score[,2] >= 0.95,1])
constraint_intersection = intersect(Petrovski_RVIS_top5_pct_gene[,1],intersect(Petrovski_haploinsufficiency_genelist[,1],Huang_Hpi_score[Huang_Hpi_score[,2] >= 0.75,1]))
constraint_union2 = c(Petrovski_RVIS_top5_pct_gene[,1],Petrovski_haploinsufficiency_genelist[,1],Huang_Hpi_score[Huang_Hpi_score[,2] >= 0.95,1],constraint_gene_tableS2[,1])

#nonconstraint_union = setdiff(c(Petrovski_RVIS_bottom_5pct_gene[,1],Huang_Hpi_score[Huang_Hpi_score[,2] <= 0.05,1], Petrovski_LoF_control_gene[,1]),constraint_union)

#nonconstraint_union = setdiff(c(Petrovski_RVIS_bottom_20pct_gene[,1],Huang_Hpi_score[Huang_Hpi_score[,2] <= 0.2,1], Petrovski_LoF_control_gene[,1]),constraint_union)

nonconstraint_union = setdiff(c(Petrovski_RVIS_bottom_10pct_gene[,1],Huang_Hpi_score[Huang_Hpi_score[,2] <= 0.1,1], Petrovski_LoF_control_gene[,1]),constraint_union)

#### gene pct for genes that are in the TADA table
gene_with_exp_mean = read.delim("../other_annotation/brainspan/expression_pct_for_TADA_genes.txt", header = TRUE, sep = "\t",stringsAsFactors = FALSE)

#Petrovski_haploinsufficiency_genelist = read.delim("../other_annotation/Petrovski_plosgen_haploinsufficiency_including_all_without_ncscore_genelist.txt",
#                                                   header = FALSE, sep = "\t", stringsAsFactors = FALSE)

synaptomeDB_presynaptics = read.delim("../other_annotation/gene_list/synaptomeDB_presynaptics.genelist", header = FALSE, sep = "\t", stringsAsFactors = FALSE )
synaptomeDB_postsynaptics = read.delim("../other_annotation/gene_list/synaptomeDB_postsynaptics.genelist", header = FALSE, sep = "\t", stringsAsFactors = FALSE )
synaptomeDB_presynaptics_activezone = read.delim("../other_annotation/gene_list/synaptomeDB_presynaptics_active_zone.genelist", header = FALSE, sep = "\t", stringsAsFactors = FALSE )
synaptomeDB_vesicles = read.delim("../other_annotation/gene_list/synaptomeDB_vesicles.genelist", header = FALSE, sep = "\t", stringsAsFactors = FALSE )
GO_synaptic_transmission = read.delim("../other_annotation/gene_list/GO_0007268_synaptic_transmission.genelist",header = FALSE, sep = "\t", stringsAsFactors = FALSE )


### Top 6 pct TADA genes used in relative risk estimation in the base-level model
TADA_top6_pct = read.delim("../other_annotation/gene_list/TADA_SNV_CNV_combined_Feb7_tadaname_null_PPA_p0_0.94_top6_pct_genes.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
TADA_top6_pct = TADA_top6_pct[,1]