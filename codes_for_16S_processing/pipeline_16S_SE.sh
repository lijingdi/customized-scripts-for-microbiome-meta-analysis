#! /bin/bash

# -----------------------------------------------------------------------------
# Script Name:  pipeline_16S_SE.sh
# Description:  This script demonstrates how to download single-end or pre-joined 16S amplicon raw data from public repository such as NCBI,process the raw data using deblur (trim all sequences to V4 region) and get the final ASV table.
# Author:       Jingdi (Judy) Li
# Date:         4 Dec 2024
# Version:      1.0
# Notes:        This pipeline should be performed on individual study basis, as primer sequences and trim length need to be manually checked for individual study
# Notes:        Always check QIIME2 website for updates of functions and pipelines here: https://docs.qiime2.org/
# -----------------------------------------------------------------------------


## V4 primer
primer_f="GTGCCAGCMGCCGCGGTAA"
primer_r="GACTACHVGGGTWTCTAAT"

#if the sequence is V3-V4, needs to trim the V3 region using seq_end
seq_end="TTACCGCGGCNGCTGGCAC"
#if the sequence is V3-V4, needs to trim the V3 region using seq_end
seq_end="ATTAGAWACCCBDGTAGTCC"

#individual study
PROJECT=
STUDY=


#download the metadata table
wget -O ./$STUDY/"$PROJECT".tsv 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession='$PROJECT'&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,library_name,library_strategy,library_source,library_selection,read_count,base_count,center_name,experiment_title,study_title,study_alias,experiment_alias,run_alias,fastq_ftp,submitted_ftp,sra_ftp,cram_index_galaxy,sample_alias,sample_title&format=tsv&download=true'

#download sequencing data, notice: download_reads.py may need some revision, based on the metadata table you downloaded above
python3 download_reads.py $PROJECT $STUDY


#fastqc, QC and manually check for each dataset
mkdir ./qiime2_result/$STUDY/ &&
mkdir ./qiime2_result/$STUDY/fastqc_out &&
fastqc -t 12 ./16S/$STUDY/raw/*.fastq.gz -o ./qiime2_result/$STUDY/fastqc_out &&
multiqc ./qiime2_result/$STUDY/fastqc_out -o ./qiime2_result/$STUDY/


#####single or merged reads
qiime tools import \
            --type SampleData[SequencesWithQuality] \
            --input-path ./16S/$STUDY/raw/ \
            --output-path ./qiime2_result/$STUDY/reads.qza \
            --input-format CasavaOneEightSingleLanePerSampleDirFmt
#fi

qiime demux summarize \
            --i-data ./qiime2_result/$STUDY/reads.qza \
            --o-visualization ./qiime2_result/$STUDY/reads_summary.qzv


#### trim primers using cutadapt
qiime cutadapt trim-paired \
            --i-demultiplexed-sequences ./qiime2_result/$STUDY/reads.qza \
            --p-cores 12 \
            --p-front $primer_f \
            --p-adapter $adapter_f  \
            --p-discard-untrimmed \
            --p-no-indels \
            --o-trimmed-sequences ./qiime2_result/$STUDY/reads_trimmed.qza

qiime demux summarize \
            --i-data ./qiime2_result/$STUDY/reads_trimmed.qza \
            --o-visualization ./qiime2_result/$STUDY/reads_trimmed_summary.qzv


###trim and deblur
qiime deblur denoise-16S \
            --i-demultiplexed-seqs ./qiime2_result/$STUDY/deblur/reads_trimmed_joined_filt.qza \
            --p-trim-length $trim_length \
            --p-left-trim-len 0 \
            --p-sample-stats \
            --p-jobs-to-start 16 \
            --p-min-reads 10 \
            --output-dir ./qiime2_result/$STUDY/deblur/denoise

# summarize, this is simply for visualizing what you have got: ASV table, sequences
qiime deblur visualize-stats \
   --i-deblur-stats ./qiime2_result/$STUDY/deblur/denoise/stats.qza \
   --o-visualization ./qiime2_result/$STUDY/deblur/denoise/deblur-stats.qzv


qiime feature-table summarize \
   --i-table ./qiime2_result/$STUDY/deblur/denoise/table.qza \
   --o-visualization ./qiime2_result/$STUDY/deblur/denoise/deblur_table_summary.qzv

qiime feature-table tabulate-seqs \
   --i-data ./qiime2_result/$STUDY/deblur/denoise/representative_sequences.qza \
   --o-visualization ./qiime2_result/$STUDY/deblur/denoise/deblur_rep-seqs.qzv


