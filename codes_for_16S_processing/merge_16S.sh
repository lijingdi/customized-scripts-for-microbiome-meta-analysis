#! /bin/bash

# -----------------------------------------------------------------------------
# Script Name:  merge_16S.sh
# Description:  This script demonstrates how to merge all the output ASV tables and sequences from pipeline_16S scripts
# Author:       Jingdi (Judy) Li
# Date:         1 Jan 2023
# Version:      1.0
# Notes:        Always check QIIME2 website for updates of functions and pipelines here: https://docs.qiime2.org/
# -----------------------------------------------------------------------------



####merge all ASV tables and sequences you got from individual studies, eg. you have two datasets for merging: McKenzie et al. 2017 and Eve et al. 2021
qiime feature-table merge \
 --i-tables ./qiime2_result/McKenzie2017/deblur/denoise/table.qza \
 --i-tables ./qiime2_result/Eve2021/deblur/denoise/table.qza \
 --o-merged-table ./merged_table_V4.qza


qiime feature-table merge-seqs \
 --i-data ./qiime2_result/McKenzie2017/deblur/denoise/representative_sequences.qza \
 --i-data ./qiime2_result/Eve2021/deblur/denoise/representative_sequences.qza \
 --o-merged-data ./merged_rep-seqs_V4.qza


####classify the features/ASVs
qiime feature-classifier classify-sklearn \
            --i-reads ./merged_rep-seqs_V4.qza \
            --i-classifier ./ref_alignments/silva-138-99-nb-classifier.qza \ ###this can be downloaded from QIIME2 website or trained by yourself, eg. it may be better to train a V4 region classifier
            --p-n-jobs 16 \
            --output-dir ./taxonomy_V4 

qiime metadata tabulate \
  --m-input-file ./taxonomy_V4/classification.qza \
  --o-visualization ./taxonomy_V4.qzv
#here can visualize rep-seq-dada2.qzv and choose ~5 ASVs to BLAST and to validate taxonomic classification.

# filter out rare ASVs/potential contaminants and those not classified at kingdom level
qiime taxa filter-table \
            --i-table ./merged_table_V4.qza \
            --i-taxonomy ./taxonomy_V4/classification.qza \
            --p-include d__ \
            --p-exclude mitochondria,chloroplast,eukaryota \
            --o-filtered-table ./table_filt_decomtam_V4.qza


qiime feature-table filter-seqs \
   --i-data ./merged_rep-seqs_V4.qza \
   --i-table ./table_filt_decomtam_V4.qza \
   --o-filtered-data ./merged_rep-seqs_final_V4.qza


# generate a phylogenetic tree
qiime fragment-insertion sepp \
   --i-representative-sequences ./merged_rep-seqs_final_V4.qza \
   --i-reference-database ./sepp-refs-gg-13-8.qza \ ###this was downloaded from QIIME2 website, need to check
   --o-tree ./asvs-tree_V4.qza \
   --o-placements ./insertion-placements_V4.qza \
   --p-threads 16


##export the tree
qiime tools export \
  --input-path ./asvs-tree_V4.qza \
  --output-path ./





