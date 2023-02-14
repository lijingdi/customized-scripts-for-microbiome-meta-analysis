### Title: Script used for microbiome meta-analysis ###
### Author: Jingdi Li, github page: https://github.com/lijingdi/customized-scripts-for-microbiome-meta-analysis ###
### data: 2023/2/14 ###

#### Differential analysis on functions and species ####

##### Differential species analysis on per study basis #####

##using Aldex2
#where each row is a different feature and each column represents a sample.
set.seed(520)
#tb <- beta_sub
tb <- beta_sub.sub
#tb <- tip_glom(beta_sub, h = 0.1)
d.x <- aldex.clr(otu_table(tb), conds = sample_data(tb)$Treatment, mc.samples = 128)
d.eff <- aldex.effect(d.x) 
d.tt <- aldex.ttest(d.x)
res.all <- data.frame(d.eff,d.tt)
nrow(res.all[res.all$wi.eBH<0.05,])
res.all$sig = rep(0,nrow(res.all))
res.all[which(res.all$wi.eBH < 0.05),"sig"] = "Significant"
res.all[which(res.all$sig == 0), "sig"] = "Not significant"
res.all_taxa <- merge(tax_table(tb), res.all, by = "row.names")

if ("heat" %in% sample_data(tb)$Treatment) {
  #heat
  res.all_taxa_asv_add <- res.all_taxa %>%
    mutate(Study = StudyID, Type = type, Group = "Heat Exposure", Level = "ASV")
  res.all_taxa_asv <- rbind(res.all_taxa_asv, res.all_taxa_asv_add)
  write.csv(res.all_taxa_asv, "/Users/lijingdi/Downloads/res_all_asv_heat.csv")
} else{
  #cold
  res.all_taxa_asv_add <- res.all_taxa %>%
    mutate(Study = StudyID, Type = type, Group = "Cold Exposure", Level = "ASV")
  res_all_cold_asv <- rbind(res_all_cold_asv, res.all_taxa_asv_add)
  write.csv(res_all_cold_asv, "/Users/lijingdi/Downloads/res_all_asv_cold.csv")
  
}

##class level
tb <- tax_glom(tb, taxrank = "Class")
d.x <- aldex.clr(otu_table(tb), conds = sample_data(tb)$Treatment, mc.samples = 128)
d.eff <- aldex.effect(d.x) 
d.tt <- aldex.ttest(d.x)
res.all <- data.frame(d.eff,d.tt)
nrow(res.all[res.all$wi.eBH<0.05,])
res.all$sig = rep(0,nrow(res.all))
res.all[which(res.all$wi.eBH < 0.05),"sig"] = "Significant"
res.all[which(res.all$sig == 0), "sig"] = "Not significant"
res.all_taxa <- merge(tax_table(tb), res.all, by = "row.names")

if ("heat" %in% sample_data(tb)$Treatment) {
  #heat
  res.all_taxa_class_add <- res.all_taxa %>%
    mutate(Study = StudyID, Type = type, Group = "Heat Exposure", Level = "Class")
  res.all_taxa_class <- rbind(res.all_taxa_class, res.all_taxa_class_add)
  write.csv(res.all_taxa_class, "/Users/lijingdi/Downloads/res_all_class_heat.csv")
} else{
  #cold
  res.all_taxa_class_add <- res.all_taxa %>%
    mutate(Study = StudyID, Type = type, Group = "Cold Exposure", Level = "Class")
  res_all_cold_class <- rbind(res_all_cold_class, res.all_taxa_class_add)
  write.csv(res_all_cold_class, "/Users/lijingdi/Downloads/res_all_class_cold.csv")
}

#check which taxa contribute most to the community difference, are those the same or different from Ancom result?
#coef <- coefficients(adon)["group1",] #not working
#top.coef <- coef[rev(order(abs(coef)))[1:20]]
#par(mar = c(3, 14, 2, 1))
#barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")

#Ancom test within each study
# library(nlme)
# library(tidyverse)
# library(ggplot2)
# library(compositions)
#library(ANCOMBC)
#class level
prepro <- feature_table_pre_process(feature_table = otu_table(tb),
                                    meta_data = sample_data(tb),
                                    sample_var = "Sample.id",
                                    group_var = NULL, out_cut = 0.05,
                                    zero_cut = 0.9, lib_cut = min(sample_sums(tb)),
                                    neg_lb = FALSE)
# #prepro$feature_table  # Preprocessed feature table
# #prepro$meta_data # Preprocessed metadata
# #prepro$structure_zeros # Structural zero info
# #Step 2, ANCOM for differential analysis
# #t_start = Sys.time()
res_ancom <- ANCOM(feature_table = prepro$feature_table,
                   meta_data = prepro$meta_data,
                   struc_zero = prepro$structure_zeros,
                   main_var = "Treatment",
                   p_adj_method = "BH",
                   alpha = 0.05,
                   adj_formula = NULL,
                   rand_formula = NULL)
res_ancom_tb <- res_ancom$out
res_ancom_tb$sig = rep(0,nrow(res_ancom_tb))
res_ancom_tb[res_ancom_tb$detected_0.9=="TRUE" | res_ancom_tb$detected_0.8=="TRUE" | res_ancom_tb$detected_0.7=="TRUE" | res_ancom_tb$detected_0.6=="TRUE", "sig"] <- "Significant"
res_ancom_tb[which(res_ancom_tb$sig == 0), "sig"] = "Not significant"
nrow(res_ancom_tb[res_ancom_tb$sig=="Significant",])

if ("heat" %in% sample_data(tb)$Treatment) {
  #heat
  res_ancom_tb_class_add <- res_ancom_tb %>%
    mutate(Study = StudyID, Type = type, Group = "Heat Exposure", level = "Class")
  res_ancom_tb_class <- rbind(res_ancom_tb_class, res_ancom_tb_class_add)
  write.csv(res_ancom_tb_class, "/Users/lijingdi/Downloads/res_ancom_tb_class.csv")
} else {
  #cold
  res_ancom_tb_class_add <- res_ancom_tb %>%
    mutate(Study = StudyID, Type = type, Group = "Cold Exposure", level = "Class")
  res_ancom_tb_class_cold <- rbind(res_ancom_tb_class_cold, res_ancom_tb_class_add)
  write.csv(res_ancom_tb_class_cold, "/Users/lijingdi/Downloads/res_ancom_tb_class_cold.csv")
}

#ASV level
#if there are too many ASVs, filter out low-abundant ASVs or low-variance OTUs, when > 1200 ASVs
tb <- beta_sub.sub
rel_abun <- transform_sample_counts(tb, function(x) x/sum(x))
remove_low_abun <- filter_taxa(rel_abun, function(x) max(x) > 0.001, TRUE) #remove low-abundant taxa
GPf = filter_taxa(rel_abun, function(x) var(x) > 1e-06, TRUE)
tb <- GPf
tb
prepro <- feature_table_pre_process(feature_table = otu_table(tb),
                                    meta_data = sample_data(tb),
                                    sample_var = "Sample.id",
                                    group_var = NULL, out_cut = 0.05,
                                    zero_cut = 0.9, lib_cut = min(sample_sums(tb)),
                                    neg_lb = FALSE)
res_ancom <- ANCOM(feature_table = prepro$feature_table,
                   meta_data = prepro$meta_data,
                   struc_zero = prepro$structure_zeros,
                   main_var = "Treatment",
                   p_adj_method = "BH",
                   alpha = 0.05,
                   adj_formula = NULL,
                   rand_formula = NULL)
res_ancom_tb <- res_ancom$out
res_ancom_tb$sig = rep(0,nrow(res_ancom_tb))
res_ancom_tb[res_ancom_tb$detected_0.9=="TRUE" | res_ancom_tb$detected_0.8=="TRUE" | res_ancom_tb$detected_0.7=="TRUE" | res_ancom_tb$detected_0.6=="TRUE", "sig"] <- "Significant"
res_ancom_tb[which(res_ancom_tb$sig == 0), "sig"] = "Not significant"
nrow(res_ancom_tb[res_ancom_tb$sig=="Significant",])

if ("heat" %in% sample_data(tb)$Treatment) {
  #heat
  res_ancom_tb_asv_add <- res_ancom_tb%>%
    mutate(Study = StudyID, Type = type, Group = "Heat Exposure", level = "ASV")
  res_ancom_tb_asv <- rbind(res_ancom_tb_asv, res_ancom_tb_asv_add)
  write.csv(res_ancom_tb_asv, "/Users/lijingdi/Downloads/res_ancom_tb_asv.csv")
} else {
  #cold
  res_ancom_tb_asv_add <- res_ancom_tb%>%
    mutate(Study = StudyID, Type = type, Group = "Cold Exposure", level = "ASV")
  res_ancom_tb_asv_cold <- rbind(res_ancom_tb_asv_cold, res_ancom_tb_asv_add)
  write.csv(res_ancom_tb_asv_cold, "/Users/lijingdi/Downloads/res_ancom_tb_asv_cold.csv")
}
###ancombc is much faster than ancom in computing time when there are >1000 features
#t_start <- Sys.time()
#ANOM_BC results would be unstable when the sample size is <5 per group
# ancom_da <- ancombc(phyloseq = tb, formula = "Treatment",
#                     p_adj_method = "BH", zero_cut = 0.99, lib_cut = 0,
#                     group = "Treatment", struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5,
#                     max_iter = 1000, conserve = TRUE, alpha = 0.05, global = TRUE)

# ancom_res_df <- data.frame(
#   Species = row.names(ancom_da$res$beta),
#   beta = unlist(ancom_da$res$beta),
#   se = unlist(ancom_da$res$se),
#   W = unlist(ancom_da$res$W),
#   p_val = unlist(ancom_da$res$p_val),
#   q_val = unlist(ancom_da$res$q_val),
#   diff_abn = unlist(ancom_da$res$diff_abn))
# 
# fdr_ancom <- ancom_res_df %>%
#   dplyr::filter(q_val < 0.05)
# 
# dim(fdr_ancom)
# 
# t_end <- Sys.time()
# t_run <- t_end - t_start
# n_taxa = ifelse(is.null(prepro$structure_zeros), nrow(prepro$feature_table), sum(apply(prepro$structure_zeros, 1, sum) == 0))
# res_df_add <- merge(res_ancom$fig$data, res_ancom$out, by = "taxa_id")
# res_df_add$Number <- n_taxa-1
# res_df_add$Study <- StudyID
# res_df_add$Type <- "HV"
# res_df_add$Group <- "Heat Exposure"
# res_df <- rbind(res_df, res_df_add)
# write.csv(res_df, "/Users/lijingdi/Downloads/res_df.csv")
# 
# 
# #step 3: Volcano plot
# # Number of taxa except structural zeros
# n_taxa = ifelse(is.null(prepro$structure_zeros), nrow(prepro$feature_table), sum(apply(prepro$structure_zeros, 1, sum) == 0))
# # Cutoff values for declaring deferentially abundant taxa
# cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
# names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")
# 
# # plot
# dat_ann = data.frame(x = min(res_ancom$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")
# 
# fig = res_ancom$fig +  
#   geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
#   geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
#             size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
# fig  

####if using ANCOM across study, can using rand_formula to stare random coefficients/slope model
#but do not need for within-study analysis


###somehow I got different results from ancom and ancombc, which is also different from that from ALdex2,
###ANCOM and ANCOMBC fail to control FDR at sample sizes <10 (https://www.nature.com/articles/s41522-020-00160-w)
###ALDEx2 can have high precision but low recall (especially when sample size is small), so ideally sample size should be 20, or at least 10 (https://europepmc.org/article/med/30021534)
###try tip_glom and tax_glom to see if there are signifiant signals

#random forest, suitable for large number of features, and modest number of observations
#so might be better if used on pooled samples rather than within-study smaller sample size and smaller number of features
#code here: https://www.listendata.com/2014/11/random-forest-with-r.html



##filtered data: detected in at least 10% of samples in each dataset, but not advisable if I want to discover 
##rare taxa, also as I am doing within-study analysis, is it necessary, though the computational time is long




##### Network analysis on differential species #####
#sparsity levels
# prev <- apply(X = otu_table(beta_sub),
#               MARGIN = 1,
#               FUN = function(x){sum(x > 0)})
# hist(prev/ncol(otu_table(beta_sub)))
#filter_prev_perc <- 0.4
#beta_sub_filt <- filter_taxa(beta_sub, function(x){sum(x > 0) > filter_prev_perc*ncol(otu_table(beta_sub))}, prune = TRUE)
#ASV level

tb <- beta_sub.sub

pseq_rel <- core(microbiome::transform(tb, "compositional"),
                 detection=0.1/100, prevalence=1/100)
pseq_rel_core <- prune_taxa(taxa_names(pseq_rel), pseq)
#Note that the confidence values were chosen as the edge stabilities under the optimal choice of the tuning parameter selected by StARS for HARMONIES, SPIEC-EASI-Glasso, and SPIEC-EASI-mb
se.mb.amgut <- spiec.easi(tb, method='mb', lambda.min.ratio=1e-2,
                          nlambda=20, pulsar.params=list(rep.num=50))
ig.mb <- adj2igraph(getRefit(se.mb.amgut), vertex.attr=list(name=taxa_names(tb)))
#plot_network(ig2.mb, tb, type='taxa', color = "Phylum")

#vsize <- rowMeans(clr(t(otu_table(tb)), 1))+6
#am.coord <- layout.fruchterman.reingold(ig.mb)
#plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
getOptMerge # edge stability
sebeta <- symBeta(getOptBeta(se.mb.amgut), mode='maxabs') #edge weights
colnames(sebeta) <- rownames(sebeta) <- colnames(t(otu_table(tb)))
hist(summary(sebeta)[,3])
dd.mb <- degree.distribution(ig.mb)
#plot(0:(length(dd.mb)-1), dd.mb, ylim=c(0,0.1), type='b', 
#     ylab="Frequency", xlab="Degree", main="Degree Distributions")
#network properties
betaMat=as.matrix(sebeta)
# We divide by two since an edge is represented by two entries in the matrix.
positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 
total=length(betaMat[betaMat!=0])/2 
#degree, betweenness and closeness
degree_calc_f(ig.mb)
#se.gl.amgut <- spiec.easi(tb, method='glasso', lambda.min.ratio=1e-2,
#                          nlambda=20, pulsar.params=list(rep.num=50))

if ("heat" %in% sample_data(tb)$Treatment) {
  #heat
  beta_sub_treat <- subset_samples(tb, Treatment == "heat")
  beta_sub_treat <- phyloseq::prune_taxa(phyloseq::taxa_sums(beta_sub_treat) > 0,beta_sub_treat)
  beta_sub_control <- subset_samples(tb, Treatment == "control")
  beta_sub_control <- phyloseq::prune_taxa(phyloseq::taxa_sums(beta_sub_control) > 0,beta_sub_control)
  
} else {
  #cold
  beta_sub_treat <- subset_samples(tb, Treatment == "cold")
  beta_sub_treat <- phyloseq::prune_taxa(phyloseq::taxa_sums(beta_sub_treat) > 0,beta_sub_treat)
  beta_sub_control <- subset_samples(tb, Treatment == "control")
  beta_sub_control <- phyloseq::prune_taxa(phyloseq::taxa_sums(beta_sub_control) > 0,beta_sub_control)
  
}
# sparcc.amgut <- SpiecEasi::sparcc(as.matrix(t(otu_table(beta_sub_treat)))) #normalization inside the function
# #get p-value, which is too computational extensive and no way to get through
# sparcc.boot <- SpiecEasi::sparccboot(as.matrix(t(otu_table(beta_sub_treat))), R = 100, ncpus = 4)
# sparcc_p <- pval.sparccboot(sparcc.boot)
# rownames(sparcc.amgut$Cor) = colnames(sparcc.amgut$Cor) = taxa_names(otu_table(beta_sub_treat))
# ## Define arbitrary threshold for SparCC correlation matrix for the graph
# sparcc.graph <- abs(sparcc.amgut$Cor) >= 0.3
# #con_value <- pulsar::natural.connectivity(sparcc.amgut$Cor) #the stability of the network
# diag(sparcc.graph) <- 0
# sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
# ig.sparcc <- adj2igraph(sparcc.graph)
# #vsize    <- rowMeans(clr(as.matrix(t(otu_table(beta_sub_treat))), 1))+6
# #am.coord <- layout.fruchterman.reingold(ig.sparcc)
# #plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")
# ###using functions from mnmnets
# #evaluate the weights on edges
# sparcc_df <- degree_calc_f(ig.sparcc)
# sparcc_df$names <- rownames(sparcc.amgut$Cor)
#Note that the confidence values were chosen as the edge stabilities under the optimal choice of the tuning parameter selected by StARS for HARMONIES, SPIEC-EASI-Glasso, and SPIEC-EASI-mb
se.mb.amgut <- spiec.easi(beta_sub_treat, method='mb', lambda.min.ratio=1e-2,
                          nlambda=20, pulsar.params=list(rep.num=50))
ig.mb <- adj2igraph(getRefit(se.mb.amgut), vertex.attr=list(name=taxa_names(tb)))
#getOptMerge(se.mb.amgut) # edge stability
sebeta <- symBeta(getOptBeta(se.mb.amgut), mode='maxabs') #edge weights
colnames(sebeta) <- rownames(sebeta) <- colnames(t(otu_table(tb)))
#filter edge weights with absolute value >0.2
hist(summary(sebeta)[,3])
dd.mb <- degree.distribution(ig.mb)
#plot(0:(length(dd.mb)-1), dd.mb, ylim=c(0,0.1), type='b', 
#     ylab="Frequency", xlab="Degree", main="Degree Distributions")
#network properties
betaMat=as.matrix(sebeta)
# We divide by two since an edge is represented by two entries in the matrix.
positive=length(betaMat[betaMat>0])/2 
negative=length(betaMat[betaMat<0])/2 
total=length(betaMat[betaMat!=0])/2 
#degree, betweenness and closeness
degree_calc_f(ig.mb)
#ASV level
if ("heat" %in% sample_data(tb)$Treatment) {
  #heat
  # sparcc_df_heat <- sparcc_df %>%
  #   mutate(Study = StudyID, Type = type, Treatment = "heat")
  # sparcc_df_combine <- rbind(sparcc_df_combine, sparcc_df_heat)
} else {
  #cold
  # sparcc_df_cold <- sparcc_df %>%
  #   mutate(Study = StudyID, Type = type, Treatment = "cold")
  # sparcc_df_combine_cold <- rbind(sparcc_df_combine_cold, sparcc_df_cold)
}
#control group
# sparcc.amgut <- SpiecEasi::sparcc(as.matrix(t(otu_table(beta_sub_control))))
# rownames(sparcc.amgut$Cor) = colnames(sparcc.amgut$Cor) = taxa_names(otu_table(beta_sub_control))
# ## Define arbitrary threshold for SparCC correlation matrix for the graph
# sparcc.graph <- abs(sparcc.amgut$Cor) >= 0.3
# diag(sparcc.graph) <- 0
# sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
# ig.sparcc <- adj2igraph(sparcc.graph)
# #vsize    <- rowMeans(clr(as.matrix(t(otu_table(beta_sub_treat))), 1))+6
# #am.coord <- layout.fruchterman.reingold(ig.sparcc)
# #plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")
# ###using functions from mnmnets
# #evaluate the weights on edges
# sparcc_df <- degree_calc_f(ig.sparcc)
# sparcc_df$names <- rownames(sparcc.amgut$Cor)
# sparcc_df_control <- sparcc_df %>%
#   mutate(Study = StudyID, Type = type, Treatment = "control")
if ("heat" %in% sample_data(tb)$Treatment) {
  #heat
  sparcc_df_combine <- rbind(sparcc_df_combine, sparcc_df_control)
  write.csv(sparcc_df_combine, "/Users/lijingdi/Downloads/sparcc_df_combine.csv")
} else {
  #cold
  sparcc_df_combine_cold <- rbind(sparcc_df_combine_cold, sparcc_df_control)
  write.csv(sparcc_df_combine_cold, "/Users/lijingdi/Downloads/sparcc_df_combine_cold.csv")
}


##class level
if ("heat" %in% sample_data(tb)$Treatment) {
  tb <- tax_glom(beta_sub.sub, taxrank = "Class")
  beta_sub_treat <- subset_samples(tb, Treatment == "heat")
  beta_sub_control <- subset_samples(tb, Treatment == "control")
} else {
  #cold
  tb <- tax_glom(beta_sub.sub, taxrank = "Class")
  beta_sub_treat <- subset_samples(tb, Treatment == "cold")
  beta_sub_control <- subset_samples(tb, Treatment == "control")
}
#treat group
# sparcc.amgut <- SpiecEasi::sparcc(as.matrix(t(otu_table(beta_sub_treat)))) #normalization inside the function
# #too computational extensive
# #sparcc.boot <- SpiecEasi::sparccboot(as.matrix(t(otu_table(beta_sub_treat))), R = 100, ncpus = 8)
# #sparcc_p <- pval.sparccboot(sparcc.boot)
# rownames(sparcc.amgut$Cor) = colnames(sparcc.amgut$Cor) = taxa_names(otu_table(beta_sub_treat))
# ## Define arbitrary threshold for SparCC correlation matrix for the graph
# sparcc.graph <- abs(sparcc.amgut$Cor) >= 0.3
# #con_value <- pulsar::natural.connectivity(sparcc.amgut$Cor) #the stability of the network
# diag(sparcc.graph) <- 0
# sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
# ig.sparcc <- adj2igraph(sparcc.graph)
# #vsize    <- rowMeans(clr(as.matrix(t(otu_table(beta_sub_treat))), 1))+6
# #am.coord <- layout.fruchterman.reingold(ig.sparcc)
# #plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")
# ###using functions from mnmnets
# #evaluate the weights on edges
# sparcc_df <- degree_calc_f(ig.sparcc)
# sparcc_df$names <- rownames(sparcc.amgut$Cor)
##class level
if ("heat" %in% sample_data(tb)$Treatment) {
  #heat
  # sparcc_df_heat <- sparcc_df %>%
  #   mutate(Study = StudyID, Type = type, Treatment = "heat")
  # sparcc_df_combine_class <- rbind(sparcc_df_combine_class, sparcc_df_heat)
} else {
  #cold
  # sparcc_df_cold <- sparcc_df %>%
  #   mutate(Study = StudyID, Type = type, Treatment = "cold")
  # sparcc_df_combine_class_cold <- rbind(sparcc_df_combine_class_cold, sparcc_df_cold)
}
#control group
# sparcc.amgut <- SpiecEasi::sparcc(as.matrix(t(otu_table(beta_sub_control))))
# rownames(sparcc.amgut$Cor) = colnames(sparcc.amgut$Cor) = taxa_names(otu_table(beta_sub_control))
# ## Define arbitrary threshold for SparCC correlation matrix for the graph
# sparcc.graph <- abs(sparcc.amgut$Cor) >= 0.3
# diag(sparcc.graph) <- 0
# sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
# ig.sparcc <- adj2igraph(sparcc.graph)
# #vsize    <- rowMeans(clr(as.matrix(t(otu_table(beta_sub_treat))), 1))+6
# #am.coord <- layout.fruchterman.reingold(ig.sparcc)
# #plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")
# ###using functions from mnmnets
# #evaluate the weights on edges
# sparcc_df <- degree_calc_f(ig.sparcc)
# sparcc_df$names <- rownames(sparcc.amgut$Cor)
# sparcc_df_control <- sparcc_df %>%
#   mutate(Study = StudyID, Type = type, Treatment = "control")
##class level
if ("heat" %in% sample_data(tb)$Treatment) {
  # sparcc_df_combine_class <- rbind(sparcc_df_combine_class, sparcc_df_control)
  # write.csv(sparcc_df_combine_class, "/Users/lijingdi/Downloads/sparcc_df_combine_class.csv")
} else {
  #cold
  # sparcc_df_combine_class_cold <- rbind(sparcc_df_combine_class_cold, sparcc_df_control)
  # write.csv(sparcc_df_combine_class_cold, "/Users/lijingdi/Downloads/sparcc_df_combine_class_cold.csv")
  #   
}



##### Functional analysis on per study basis #####
#functional annotation
ko_tb <- read.csv("/Users/lijingdi/Downloads/picrust2_out_pipeline_all/KO_metagenome_out/pred_metagenome_unstrat.tsv", header = 1, row.names = 1, sep = "\t")
ec_tb <- read.csv("/Users/lijingdi/Downloads/picrust2_out_pipeline_all/EC_metagenome_out/pred_metagenome_unstrat.tsv", header = 1, row.names = 1, sep = "\t")
path_tb <- read.csv("/Users/lijingdi/Downloads/picrust2_out_pipeline_all/pathways_out/path_abun_unstrat.tsv", header = 1, row.names = 1, sep = "\t")
meta_beta_subset <- read.csv("/Users/lijingdi/Downloads/meta_beta_subset.txt", header = 1, row.names = 1, sep = "\t")

#posadas
#ko_tb <- read.csv("/Users/lijingdi/Downloads/Posadas2021/picrust2_result/KO_metagenome_out/pred_metagenome_unstrat.tsv", header = 1, row.names = 1, sep = "\t")
#ec_tb <- read.csv("/Users/lijingdi/Downloads/Posadas2021/picrust2_result/EC_metagenome_out/pred_metagenome_unstrat.tsv", header = 1, row.names = 1, sep = "\t")
#path_tb <- read.csv("/Users/lijingdi/Downloads/Posadas2021/picrust2_result/path_abun_unstrat.tsv", header = 1, row.names = 1, sep = "\t")
ko_phy <- phyloseq(otu_table(ko_tb, taxa_are_rows = TRUE), phyloseq::sample_data(meta_beta_subset))
ec_phy <- phyloseq(otu_table(ec_tb, taxa_are_rows = TRUE), phyloseq::sample_data(meta_beta_subset))
path_phy <- phyloseq(otu_table(path_tb, taxa_are_rows = TRUE), phyloseq::sample_data(meta_beta_subset))


ko_sub <- phyloseq::subset_samples(ko_phy, Study == StudyID)
ko_sub <- phyloseq::subset_samples(ko_sub,Sampling_part == "gut_lab")
ko_sub <- phyloseq::prune_taxa(phyloseq::taxa_sums(ko_sub) > 0, ko_sub)
ko_sub

ec_sub <- phyloseq::subset_samples(ec_phy, Study == StudyID)
ec_sub <- phyloseq::subset_samples(ec_sub, Sampling_part == "gut_lab")
ec_sub <- phyloseq::prune_taxa(phyloseq::taxa_sums(ec_sub) > 0, ec_sub)
ec_sub

path_sub <- phyloseq::subset_samples(path_phy, Study == StudyID)
path_sub <- phyloseq::subset_samples(path_sub, Sampling_part == "gut_lab")
path_sub <- phyloseq::prune_taxa(phyloseq::taxa_sums(path_sub) > 0, path_sub)
path_sub


func_sub <- ko_sub
func <- "ko"

func_sub <- ec_sub
func <- "ec"

func_sub <- path_sub
func <- "path"

##wilcox test
path_merge_tb <- merge(t(otu_table(func_sub)), meta_beta_subset, by = "row.names")
row.names(path_merge_tb) <- path_merge_tb$Row.names
path_merge_tb$Row.names <- NULL
modelList<-list()
for (i in 1:length(row.names(otu_table(func_sub)))) {
  fmla <- formula(paste(paste0("`", row.names(otu_table(func_sub))[i], "`"), "~ Treatment"))
  modelList[[i]] <- wilcox.test(fmla, data = path_merge_tb, alternative ="two.sided")$p.value
}
dt_p <- data.frame(Tax = row.names(otu_table(func_sub)), P_value = unlist(modelList))
dt_p$P_adjusted <- p.adjust(dt_p$P_value, method = 'BH')
nrow(dt_p[dt_p$P_adjusted<= 0.05,])
sum_path <- psmelt(func_sub) %>%
  as_tibble %>%
  group_by(OTU, Treatment) %>%
  summarize(Mean = mean(Abundance), SD = sd(Abundance)) %>%
  ungroup
if ("heat" %in% sample_data(func_sub)$Treatment) {
  #save temp file, heat
  write.csv(sum_path, paste0("/Users/lijingdi/Downloads/tempp/", StudyID, "_sum_path_", func, "_", type, ".csv"))
} else {
  #save cold
  write.csv(sum_path, paste0("/Users/lijingdi/Downloads/tempp/", StudyID, "_sum_path_", func, "_", type, "_cold.csv"))
}

#export result
sum_path_diff <- matrix(0, nrow = 0, ncol = 2)
sum_path_diff$OTU <- unique(sum_path$OTU)
sum_path_diff <- as.data.frame(sum_path_diff)

if ("heat" %in% sample_data(func_sub)$Treatment) {
  #heat
  for (i in sum_path_diff$OTU) {
    sum_path_diff[sum_path_diff$OTU==i, "control_minus_heat"] <- sum_path[sum_path$OTU==i & sum_path$Treatment=="control", "Mean"]-sum_path[sum_path$OTU==i & sum_path$Treatment=="heat", "Mean"]
  }
  dt_path_merge_add <- merge(dt_p, sum_path_diff, by.x = "Tax", by.y = "OTU") %>%
    mutate(Significance = ifelse(P_adjusted <= 0.05, "Significant", "Non-significant"), Study=StudyID, Type = type)
  
  dt_path_merge <- read.csv(paste0("/Users/lijingdi/Downloads/dt_path_merge_", func, ".csv"), header = 1, row.names = 1)
  dt_path_merge <- rbind(dt_path_merge, dt_path_merge_add)
  write.csv(dt_path_merge, paste0("/Users/lijingdi/Downloads/dt_path_merge_", func, ".csv"))
} else {
  #cold
  for (i in sum_path_diff$OTU) {
    sum_path_diff[sum_path_diff$OTU==i, "control_minus_cold"] <- sum_path[sum_path$OTU==i & sum_path$Treatment=="control", "Mean"]-sum_path[sum_path$OTU==i & sum_path$Treatment=="cold", "Mean"]
  }
  dt_path_merge_add <- merge(dt_p, sum_path_diff, by.x = "Tax", by.y = "OTU") %>%
    mutate(Significance = ifelse(P_adjusted <= 0.05, "Significant", "Non-significant"), Study=StudyID, Type = type)
  
  dt_path_merge_cold <- read.csv(paste0("/Users/lijingdi/Downloads/dt_path_merge_", func, "_cold.csv"), header = 1, row.names = 1)
  dt_path_merge_cold <- rbind(dt_path_merge_cold, dt_path_merge_add)
  write.csv(dt_path_merge_cold, paste0("/Users/lijingdi/Downloads/dt_path_merge_", func, "_cold.csv"))
}
#write.csv(dt_path_merge, paste0("/Users/lijingdi/Downloads/Posadas2021/dt_path_", func, "_", type, "_merge.csv"))



