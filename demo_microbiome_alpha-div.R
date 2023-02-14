### Title: Script used for microbiome meta-analysis ###
### Author: Jingdi Li, github page: https://github.com/lijingdi/customized-scripts-for-microbiome-meta-analysis ###
### data: 2023/2/14 ###


#### Alpha diversity effect-size calculation ####
#reference, check for details: https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/ 
#https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121 
###### load library ######
library(phyloseq)
library(ape)
library(tidyr)
library(qiime2R)
library(vegan)
library(ggplot2)
library(picante)
library(dplyr)
library(tibble)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(esc)
library(effectsize)
library(microbiome)

###### qza_to_phyloseq function, revise accordingly ######
qza_to_phyloseq_new<-function(features,tree,taxonomy,metadata, tmp){
  
  if(missing(features) & missing(tree) & missing(taxonomy) & missing(metadata)){
    stop("At least one required artifact is needed (features/tree/taxonomy/) or the metadata.")
  }
  
  if(missing(tmp)){tmp <- tempdir()}
  
  argstring<-""
  
  if(!missing(features)){
    features<-read_qza(features, tmp=tmp)$data
    argstring<-paste(argstring, "otu_table(features, taxa_are_rows=T),")
  }
  
  if(!missing(taxonomy)){
    taxonomy<-read_qza(taxonomy, tmp=tmp)$data
    taxonomy<-parse_taxonomy(taxonomy)
    taxonomy<-as.matrix(taxonomy)
    argstring<-paste(argstring, "tax_table(taxonomy),")
  }
  
  if(!missing(tree)){
    tree<-read_qza(tree, tmp=tmp)$data
    argstring<-paste(argstring, "phy_tree(tree),")
  }
  
  if(!missing(metadata)){
    if(is_q2metadata(metadata)){
      metadata<-read_q2metadata(metadata)
      rownames(metadata)<-metadata$SampleID
      metadata$SampleID<-NULL
    } else{
      metadata<-read.table(metadata, row.names=1, sep='\t', header=TRUE)
    }
    argstring<-paste(argstring, "sample_data(metadata),")
  }
  
  argstring<-gsub(",$","", argstring) #remove trailing ","
  
  physeq<-eval(parse(text=paste0("phyloseq(",argstring,")")))
  
  return(physeq)
}


#Two implications to consider are that 
#(1) sampling with replacement is faster and more memory efficient as currently implemented; and 
#(2), sampling with replacement means that there is a chance that the number of reads for a given 
#OTU in a given sample could be larger than the original count value, as opposed to sampling 
#without replacement where the original count value is the maximum possible. 

###### load output file from QIIME2 ######
meta_combined <- read.table("/PATH/TO/YOUR_METADATA_TABLE.tsv", sep = "\t", header = 1)
row.names(meta_combined) <- meta_combined$Sample.id
StudyID = "Hartman2020"

###### convert QIIME2 output to a phyloseq object ######
assign("beta", qza_to_phyloseq_new(features = "/PATH/To/merged_table.qza", taxonomy = "/PATH/To/classification.qza", metadata = "/PATH/To//meta_beta_subset.txt"))
tree <- read.tree("/PATH/To/tree.nwk")
phy_tree(beta) <- tree

get(StudyID) ### now this is a phyloseq object for dataset Hartman2020
deblur_sub_raw <- subset_samples(get(StudyID), Study == StudyID & Treatment %in% c("heat", "control") & Sampling_point == "end")

###### rarefy ######
#check read count per sample, and determine the threshold of filtering by checking rarefraction curve
sort(sample_sums(deblur_sub_raw))
raremin <- min(sample_sums(deblur_sub_raw)) ### this is the minimum number of read count, if filtering by this number, then no sample will be filtered
#https://fromthebottomoftheheap.net/2015/04/16/drawing-rarefaction-curves-with-custom-colours/ custom color
rarecurve(t(otu_table(deblur_sub_raw)), sample = raremin, step = 50, cex = 0.5) ## this is for visualizing the rarefraction curve

rare_value <- raremin ## this value needs to be adjusted according to rarefraction curve
deblur_sub <- rarefy_even_depth(deblur_sub_raw, rngseed = 520, sample.size = rare_value, replace = F) ## rarefy 
rarecurve(t(otu_table(deblur_sub)), step = 50, cex = 0.5) ## visualizing rarefraction curve after rarefy


##### Alpha diversity analysis #####

#ASV to genus level, not needed if performing analysis on ASV level
#deblur_sub_genus <- tax_glom(deblur_sub, taxrank = "Genus")

####### calculate alpha diversity
deblur_alpha <- data.frame(
  "Observed" = phyloseq::estimate_richness(deblur_sub, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(deblur_sub, measures = "Shannon"),
  "PD" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(deblur_sub)))), tree = phy_tree(beta_sub.sub)) [, 1],
  #  "PD" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(deblur_sub)))), tree = phyloseq::phy_tree(deblur_sub))[, 1],
  "Treatment" = phyloseq::sample_data(deblur_sub)$Treatment)
deblur_alpha$Evenness <- deblur_alpha$Shannon/log(deblur_alpha$Observed)


####### visualizing alpha diversity (can be skipped) #######
alpha_plot <- deblur_alpha %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "PD", "Evenness")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "PD", "Evenness"))) %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Evenness")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Evenness"))) %>%
  ggplot(aes(x = Treatment, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Treatment), height = 0, width = .2) +
  scale_color_manual(values = c("#0072B2","#D55E00")) +
  labs(x = "Treatment", y = "Alpha diversity estimates") +
  facet_wrap(~ metric, scales = "free") +
  theme_bw()  +
  theme(panel.grid.minor.y = element_blank(), legend.position = "none") 
alpha_plot


###### effect size calculation ######
# distribution of alpha-div, normality test within group
#the assumption for t-test is that the residuals are normally distributed.
#if the two populations differ in means, the distribution of the combined groups might be 
#bi-model, even if both are normally distributed on their own
melt_deblur <- 
  deblur_alpha %>%
  melt(., id.vars = "Treatment", measure.vars = c("Observed", "Shannon", "PD", "Evenness"))
#or use `melt(data, id = c("Treatment"))`
head(melt_deblur)
##check for density plot
melt_deblur %>% 
  ggplot(., aes(x = value, color = Treatment)) +
  geom_density() +
  facet_wrap(~variable, scales = "free") 
##check for qqplot
melt_deblur %>% 
  ggqqplot(., x = "value", color = "Treatment") +
  facet_wrap(~variable, scales = "free") 
##shapiro.test
melt_deblur %>% 
  group_by(Treatment, variable) %>%
  group_map(~ shapiro.test(.x$value)) ##returns a list of results, to return a data.frame, use group_modify
#p-value<0.05, reject the null hypothesis of normality

###### separate groups ######
#need data.table v1.9.6
#here replace "heat" with your desired grouping factor
packageVersion("data.table")
if ("heat" %in% melt_deblur$Treatment) {
  alpha_heat <-
    melt_deblur %>%
    group_by(variable, Treatment) %>%
    summarise(mean = mean(value), sd = sd(value), n = n()) %>%
    ungroup() %>%
    mutate(id=paste0(variable, "_", Treatment)) %>%
    as.data.table(.) %>%
    data.table::dcast(.,  Treatment ~ id, value.var=c("mean", "sd", "n")) %>%
    select(-Treatment) %>%
    summarise_all(funs(sum_NA)) %>%
    mutate(Study=StudyID, Type = type)
} else {
  alpha_heat <-
    melt_deblur %>%
    group_by(variable, Treatment) %>%
    summarise(mean = mean(value), sd = sd(value), n = n()) %>%
    ungroup() %>%
    mutate(id=paste0(variable, "_", Treatment)) %>%
    as.data.table(.) %>%
    data.table::dcast(.,  Treatment ~ id, value.var=c("mean", "sd", "n")) %>%
    select(-Treatment) %>%
    dplyr::summarise_all(funs(sum_NA)) %>%
    mutate(Study=StudyID, Type = type)

#### calculate sed.pd, can skip ####  
alpha.ses.pd <- ses.pd(samp = t(as.matrix(as.data.frame(phyloseq::otu_table(deblur_sub)))), #community data matrix has sample name as row.names
                         tree = phyloseq::phy_tree(deblur_sub),
                         null.model = "taxa.labels",
                         runs = 999)
alpha_ses <- alpha.ses.pd %>%
    mutate(Treatment = sample_data(deblur_sub)$Treatment) %>%
    group_by(Treatment) %>%
    summarise(mean = mean(pd.obs.z), sd = sd(pd.obs.z), n = n()) %>%
    ungroup() %>%
    as.data.table(.) %>%
    mutate(id = paste0("ses.pd", "_", Treatment)) %>%
    data.table::dcast(., Treatment ~ id, value.var=c("mean", "sd", "n")) %>%
    select(-Treatment) %>%
    dplyr::summarise_all(funs(sum_NA)) %>%
    mutate(StudyID=StudyID, Type = type)

### combine values for different alpha-div metric
alpha_heat <- rbind(alpha_heat, alpha_ses)

#### calculate effect-size, hedge's g
alpha_heat_ob <- escalc(measure = "SMD", m1i = mean_Observed_control, sd1i = sd_Observed_control, n1i = n_Observed_control, m2i = mean_Observed_heat, sd2i = sd_Observed_heat, n2i = n_Observed_heat, data = alpha_heat)
alpha_heat_sha <- escalc(measure = "SMD", m1i = mean_Shannon_control, sd1i = sd_Shannon_control, n1i = n_Shannon_control, m2i = mean_Shannon_heat, sd2i = sd_Shannon_heat, n2i = n_Shannon_heat, data = alpha_heat)
alpha_heat_pd <- escalc(measure = "SMD", m1i = mean_PD_control, sd1i = sd_PD_control, n1i = n_PD_control, m2i = mean_PD_heat, sd2i = sd_PD_heat, n2i = n_PD_heat, data = alpha_heat)
alpha_heat_even <- escalc(measure = "SMD", m1i = mean_Evenness_control, sd1i = sd_Evenness_control, n1i = n_Evenness_control, m2i = mean_Evenness_heat, sd2i = sd_Evenness_heat, n2i = n_Evenness_heat, data = alpha_heat)
alpha_heat_ses <- escalc(measure = "SMD", m1i = mean_ses.pd_control, sd1i = sd_ses.pd_control, n1i = n_ses.pd_control, m2i = mean_ses.pd_heat, sd2i = sd_ses.pd_heat, n2i = n_ses.pd_heat, data = alpha_ses)


####### seperate different groups that you want to compare, eg. heat vs. control #######
melt_ob_ctr <- melt_deblur %>%
  filter(variable == "Observed" & Treatment == "control") %>%
  summarise(mean_ob_ctr = mean(value), sd_ob_ctr = sd(value), n_ob_ctr = n())

melt_ob_heat <- melt_deblur %>%
  filter(variable == "Observed" & Treatment == "heat") %>%
  summarise(mean_ob_heat = mean(value), sd_ob_heat = sd(value), n_ob_heat = n())

melt_sha_ctr <- melt_deblur %>%
  filter(variable == "Shannon" & Treatment == "control") %>%
  summarise(mean_sha_ctr = mean(value), sd_sha_ctr = sd(value), n_sha_ctr = n())

melt_sha_heat <- melt_deblur %>%
  filter(variable == "Shannon" & Treatment == "heat") %>%
  summarise(mean_sha_heat = mean(value), sd_sha_heat = sd(value), n_sha_heat = n())

melt_pd_ctr <- melt_deblur %>%
  filter(variable == "PD" & Treatment == "control") %>%
  summarise(mean_pd_ctr = mean(value), sd_pd_ctr = sd(value), n_pd_ctr = n())

melt_pd_heat <- melt_deblur %>%
  filter(variable == "PD" & Treatment == "heat") %>%
  summarise(mean_pd_heat = mean(value), sd_pd_heat = sd(value), n_pd_heat = n())

alpha_heat <- cbind(melt_ob_ctr, melt_ob_heat, melt_sha_ctr, melt_sha_heat, melt_pd_ctr, melt_pd_heat)


####### calculate effect-size for alpha-diversity #######
# escalc calculate effect-size, can also add evenness and other metrics you are interested
#default, vtype="LS", which uses the usual large-sample approximation to compute the sampling variance, "UB" provides unbiased estimates of the sampling variance. 
alpha_heat_ob <- escalc(measure = "SMD", m1i = mean_ob_ctr, sd1i = sd_ob_ctr, n1i = n_ob_ctr, m2i = mean_ob_heat, sd2i = sd_ob_heat, n2i = n_ob_heat, data = alpha_heat)
alpha_heat_sha <- escalc(measure = "SMD", m1i = mean_sha_ctr, sd1i = sd_sha_ctr, n1i = n_sha_ctr, m2i = mean_sha_heat, sd2i = sd_sha_heat, n2i = n_sha_heat, data = alpha_heat)
alpha_heat_pd <- escalc(measure = "SMD", m1i = mean_pd_ctr, sd1i = sd_pd_ctr, n1i = n_pd_ctr, m2i = mean_pd_heat, sd2i = sd_pd_heat, n2i = n_pd_heat, data = alpha_heat)

#pool effect size together
alpha_eff<-matrix(NA,ncol=20,nrow=0)
colnames(alpha_eff) <- c("es_ob", "se_ob", "variance_ob", "lower_ci_ob", "upper_ci_ob", "weight_ob", "es_shannon", "se_shannon", "variance_shannon", "lower_ci_shannon", "upper_ci_shannon", "weight_shannon", "es_pd", "se_pd", "variance_pd", "lower_ci_pd", "upper_ci_pd", "weight_pd", "Studyid", "Type")
#export table
write.csv(alpha_eff, "/Users/user/Downloads/alpha_eff.csv")


