#! /bin/bash

# -----------------------------------------------------------------------------
# Script Name:  pipeline_16S_PE.sh
# Description:  This script demonstrates how to download paired-end 16S amplicon raw data from public repository such as NCBI,process the raw data using deblur (trim all sequences to V4 region) and get the final ASV table.
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

#individual study, eg. Sun et al. 2020 has a dataset with BioProject ID PRJNA526611
PROJECT="PRJNA526611"
STUDY="Sun2020"



#download the metadata table
wget -O ./$STUDY/"$PROJECT".tsv 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession='$PROJECT'&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,library_name,library_strategy,library_source,library_selection,read_count,base_count,center_name,experiment_title,study_title,study_alias,experiment_alias,run_alias,fastq_ftp,submitted_ftp,sra_ftp,cram_index_galaxy,sample_alias,sample_title&format=tsv&download=true'

#download sequencing data, notice: download_reads.py may need some revision, based on the metadata table you downloaded above
python3 download_reads.py $PROJECT $STUDY

#fastqc for QC, the output should be manually checked 
mkdir ./qiime2_result/$STUDY/ &&
mkdir ./qiime2_result/$STUDY/fastqc_out &&
fastqc -t 12 ./$STUDY/raw/*.fastq.gz -o ./qiime2_result/$STUDY/fastqc_out &&
multiqc ./qiime2_result/$STUDY/fastqc_out -o ./qiime2_result/$STUDY/


####paired-end data
##import data to qiime2
qiime tools import \
            --type SampleData[PairedEndSequencesWithQuality] \
            --input-path ./$STUDY/raw/ \
            --output-path ./qiime2_result/$STUDY/reads.qza \
            --input-format CasavaOneEightSingleLanePerSampleDirFmt


qiime demux summarize \
            --i-data ./qiime2_result/$STUDY/reads.qza \
            --o-visualization ./qiime2_result/$STUDY/reads_summary.qzv


#### trim primers using cutadapt
##### for V3-V4:--p-adapter-r $seq_end \ for V4-V5: --p-adapter-f $seq_end
qiime cutadapt trim-paired \
            --i-demultiplexed-sequences ./qiime2_result/$STUDY/reads.qza \
            --p-cores 12 \
            --p-front-f $primer_f \
            --p-front-r $primer_r \
            --p-discard-untrimmed \
            --p-no-indels \
            --o-trimmed-sequences ./qiime2_result/$STUDY/reads_trimmed_pre.qza \
			--verbose


qiime cutadapt trim-paired \
            --i-demultiplexed-sequences ./qiime2_result/$STUDY/reads_trimmed_pre.qza \
            --p-cores 12 \
	    	--p-adapter-r $seq_end \
            --p-no-indels \
            --o-trimmed-sequences ./qiime2_result/$STUDY/reads_trimmed.qza \
			--verbose


qiime demux summarize \
            --i-data ./qiime2_result/$STUDY/reads_trimmed.qza \
            --o-visualization ./qiime2_result/$STUDY/reads_trimmed_summary.qzv


####join reads using deblur
if [ ! -d ./qiime2_result/$STUDY/deblur ]; then
        mkdir ./qiime2_result/$STUDY/deblur
fi

qiime vsearch join-pairs \
            --i-demultiplexed-seqs ./qiime2_result/$STUDY/reads_trimmed.qza \
            --o-joined-sequences ./qiime2_result/$STUDY/deblur/reads_trimmed_joined.qza \
	         --verbose \
	         --p-allowmergestagger

# summarize the joined reads            
qiime demux summarize \
            --i-data ./qiime2_result/$STUDY/deblur/reads_trimmed_joined.qza \
            --o-visualization ./qiime2_result/$STUDY/deblur/reads_trimmed_joined_summary.qzv

# filter out low-quality reads, here you can also add customized filtering parameter
qiime quality-filter q-score \
            --i-demux ./qiime2_result/$STUDY/deblur/reads_trimmed_joined.qza \
            --o-filter-stats ./qiime2_result/$STUDY/deblur/filt_stats.qza \
            --o-filtered-sequences ./qiime2_result/$STUDY/deblur/reads_trimmed_joined_filt.qza


# summarize, It is a good idea at this point just to verify that there haven't been any substantial losses of reads, before going through the whole ASV process, at either the joining or quality-filtering steps above:
qiime demux summarize \
            --i-data ./qiime2_result/$STUDY/deblur/reads_trimmed_joined_filt.qza \
            --o-visualization ./qiime2_result/$STUDY/deblur/reads_trimmed_joined_filt_summary.qzv


#####here by manually checking the reads_trimmed_joined_filt_summary.qzv file
#####you can decide the trim_length,for V4 region, 250bp is normally fine, so you can set trim_length=250

###trim and deblur, parameters you can adjust: --p-jobs-to-start, --p-min-reads
qiime deblur denoise-16S \
            --i-demultiplexed-seqs ./qiime2_result/$STUDY/deblur/reads_trimmed_joined_filt.qza \
            --p-trim-length $trim_length \
            --p-left-trim-len 0 \
            --p-sample-stats \
            --p-jobs-to-start 16 \
            --p-min-reads 10 \
            --output-dir ./qiime2_result/$STUDY/deblur/denoise

# summarize, this is for visualization of ASV table and sequences output
qiime deblur visualize-stats \
   --i-deblur-stats ./qiime2_result/$STUDY/deblur/denoise/stats.qza \
   --o-visualization ./qiime2_result/$STUDY/deblur/denoise/deblur-stats.qzv


qiime feature-table summarize \
   --i-table ./qiime2_result/$STUDY/deblur/denoise/table.qza \
   --o-visualization ./qiime2_result/$STUDY/deblur/denoise/deblur_table_summary.qzv

qiime feature-table tabulate-seqs \
   --i-data ./qiime2_result/$STUDY/deblur/denoise/representative_sequences.qza \
   --o-visualization ./qiime2_result/$STUDY/deblur/denoise/deblur_rep-seqs.qzv



