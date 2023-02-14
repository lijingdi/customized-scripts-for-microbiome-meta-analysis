### Title: Script used for microbiome meta-analysis ###
### Author: Jingdi Li, github page: https://github.com/lijingdi/customized-scripts-for-microbiome-meta-analysis ###
### data: 2023/2/14 ###


#### Beta diversity and dispersion effect-size calculation ####

##### load data, functions can be loaded from alpha-div script #####
StudyID = "Hartman2020" 

# convert QIIME2 output to a phyloseq object
assign("beta", qza_to_phyloseq_new(features = "/PATH/To/merged_table.qza", taxonomy = "/PATH/To/classification.qza", metadata = "/PATH/To//meta_beta_subset.txt"))
tree <- read.tree("/PATH/To//tree.nwk")
phy_tree(beta) <- tree
beta_sub <- phyloseq::subset_samples(beta, Study == StudyID)
beta_sub <- phyloseq::prune_taxa(phyloseq::taxa_sums(beta_sub) > 0, beta_sub) 
beta_sub

### try all available normalization methods, can skip some 
####unifrac and wunifrac on rarefied data
beta_sub.rare <- phyloseq::rarefy_even_depth(beta_sub, rngseed = 520, sample.size = min(sample_sums(beta_sub)), replace = F)
beta_subr.rel <- phyloseq::transform_sample_counts(beta_sub.rare, function(x){x / sum(x)})
unifrac_dist = phyloseq::distance(beta_subr.rel, method="unifrac", weighted=FALSE)
wunifrac_dist = phyloseq::distance(beta_subr.rel, method="unifrac", weighted=T)
ado_u <- adonis(unifrac_dist ~ sample_data(beta_subr.rel)$Host_genotype + sample_data(beta_subr.rel)$Treatment, permutations = 9999) 
ado_w <- adonis(wunifrac_dist ~ sample_data(beta_subr.rel)$Host_genotype + sample_data(beta_subr.rel)$Treatment, permutations = 9999)
#clr for whole dataset
beta_sub.clr <- microbiome::transform(beta_sub, 'clr')
euclidean_dist <- phyloseq::distance(beta_sub.clr, method="euclidean")
ado_e <- adonis(euclidean_dist ~ sample_data(beta_sub.clr)$Host_genotype + sample_data(beta_sub.clr)$Treatment, permutations = 9999) 
#philr for whole dataset
physeq <- phyloseq::transform_sample_counts(beta_sub, function(x) {x+1})
physeq.philr <- philr(t(otu_table(physeq)), phyloseq::phy_tree(physeq), part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt', return.all = F)
physeq.dist <- dist(physeq.philr, method="euclidean")
ado_p <- adonis(physeq.dist~ sample_data(physeq)$Host_genotype + sample_data(physeq)$Treatment, permutations = 9999) 


####### calculate effect-size for beta-div
eff_u <- MicEco::adonis_OmegaSq(ado_u) #negative omegaSQ are treated as 0
eff_w <- MicEco::adonis_OmegaSq(ado_w)
eff_e <- MicEco::adonis_OmegaSq(ado_e)
eff_p <- MicEco::adonis_OmegaSq(ado_p)
#export the result
ava_df_pre <- matrix(NA, ncol = 14, nrow = 0)
colnames(ava_df_pre) <- c("Study", "Type", "Normalization", "Dissimilarity metric", "Grouping", "F", "Pr(>F)", "R2", "parOmegaSq", "dispr", "dispr_avg_dist_to_meadian_diff", "dispr_lwr", "dispr_upr", "dispr_p_value")
if (exists("eff_p")) {
  ava_df_pre_add <- do.call("rbind", list(
    c(StudyID, "combined", "rarefy", "unweighted", group, eff_u$aov.tab$F.Model[1], eff_u$aov.tab$`Pr(>F)`[1], eff_u$aov.tab$R2[1], eff_u$aov.tab$parOmegaSq[1], row.names(tukey_u$group), tukey_u$group[, "diff"], tukey_u$group[, "lwr"], tukey_u$group[, "upr"], tukey_u$group[, "p adj"]),
    c(StudyID, "combined", "rarefy", "weighted", group, eff_w$aov.tab$F.Model[1], eff_w$aov.tab$`Pr(>F)`[1], eff_w$aov.tab$R2[1], eff_w$aov.tab$parOmegaSq[1], row.names(tukey_w$group), tukey_w$group[, "diff"], tukey_w$group[, "lwr"], tukey_w$group[, "upr"], tukey_w$group[, "p adj"]),
    c(StudyID, "combined", "clr", "euclidean", group, eff_e$aov.tab$F.Model[1], eff_e$aov.tab$`Pr(>F)`[1], eff_e$aov.tab$R2[1], eff_e$aov.tab$parOmegaSq[1], row.names(tukey_e$group), tukey_e$group[, "diff"], tukey_e$group[, "lwr"], tukey_e$group[, "upr"], tukey_e$group[, "p adj"]),
    c(StudyID, "combined", "philr", "euclidean", group, eff_p$aov.tab$F.Model[1], eff_p$aov.tab$`Pr(>F)`[1], eff_p$aov.tab$R2[1], eff_p$aov.tab$parOmegaSq[1], row.names(tukey_p$group), tukey_p$group[, "diff"], tukey_p$group[, "lwr"], tukey_p$group[, "upr"], tukey_p$group[, "p adj"])
  ))} else {
    ava_df_pre_add <- do.call("rbind", list(
      c(StudyID, "combined", "rarefy", "unweighted", group, eff_u$aov.tab$F.Model[1], eff_u$aov.tab$`Pr(>F)`[1], eff_u$aov.tab$R2[1], eff_u$aov.tab$parOmegaSq[1], row.names(tukey_u$group), tukey_u$group[, "diff"], tukey_u$group[, "lwr"], tukey_u$group[, "upr"], tukey_u$group[, "p adj"]),
      c(StudyID, "combined", "rarefy", "weighted", group, eff_w$aov.tab$F.Model[1], eff_w$aov.tab$`Pr(>F)`[1], eff_w$aov.tab$R2[1], eff_w$aov.tab$parOmegaSq[1], row.names(tukey_w$group), tukey_w$group[, "diff"], tukey_w$group[, "lwr"], tukey_w$group[, "upr"], tukey_w$group[, "p adj"]),
      c(StudyID, "combined", "clr", "euclidean", group, eff_e$aov.tab$F.Model[1], eff_e$aov.tab$`Pr(>F)`[1], eff_e$aov.tab$R2[1], eff_e$aov.tab$parOmegaSq[1], row.names(tukey_e$group), tukey_e$group[, "diff"], tukey_e$group[, "lwr"], tukey_e$group[, "upr"], tukey_e$group[, "p adj"]),
      c(StudyID, "combined", "philr", "euclidean", group, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
    ))
  }
colnames(ava_df_pre_add) <- c("Study", "Type", "Normalization", "Dissimilarity metric", "Grouping", "F", "Pr(>F)", "R2", "parOmegaSq", "dispr", "dispr_avg_dist_to_meadian_diff", "dispr_lwr", "dispr_upr", "dispr_p_value")
ava_df_pre_add[ava_df_pre_add[,"parOmegaSq"] < 0, "parOmegaSq"] <- 0
ava_df_pre <- rbind(ava_df_pre, ava_df_pre_add) ### this is the output table for beta-diversity effect sizes


###### calculate effect-size for beta-dispersion ######
library(esc)
f_value <- ava$`F value`[1]
sampled <- as.data.frame(sample_data(beta_sub))
grp1 <- nrow(sampled[sampled$Treatment == "heat", ])
grp2 <- nrow(sampled[sampled$Treatment == "control", ])
esc::esc_f(f = f_value, grp1n = grp1, grp2n = grp2, es.type="g")

