#! /bin/bash

PROJECT="PRJNA673151" #here project ID is the NCBI bioproject ID from which we will download raw reads data
STUDY="StudyName" # for each study, the name is different
### need to collect the info about primer pair from the literature
primer_f="ACTCCTACGGGAGGCAGCAG" ### forward primer 
primer_r="GGACTACHVGGGTWTCTAAT" ### reverse primer


# 1. download the metadata table and raw reads
### create a folder for the study if not existed already
if [ ! -d /data/zool-temp-micro/pemb6000/16S/$STUDY ]; then
        mkdir /data/zool-temp-micro/pemb6000/16S/$STUDY && mkdir /data/zool-temp-micro/pemb6000/16S/$STUDY/raw
fi
### download the metadata table for the study
wget -O /data/zool-temp-micro/pemb6000/16S/$STUDY/"$PROJECT".tsv 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession='$PROJECT'&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,library_name,library_strategy,library_source,library_selection,read_count,base_count,center_name,experiment_title,study_title,study_alias,experiment_alias,run_alias,fastq_ftp,submitted_ftp,sra_ftp,cram_index_galaxy,sample_alias,sample_title&format=tsv&download=true'

### download raw data and rename it for fitting qiime2 format, I worte a custom script in python for this
for i in $(find /data/zool-temp-micro/pemb6000/16S/$STUDY/raw/ -type f -size 0c); do rm $i;done &&
python3 download_reads_for_qiime2.py $PROJECT $STUDY 


# 2. processing the raw reads
### 2.1 fastqc, for checking the quality of reads, and use multiqc for visualization
fastqc -t 12 /data/zool-temp-micro/pemb6000/16S/$STUDY/raw/*.fastq.gz -o /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/fastqc_out &&
multiqc /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/fastqc_out -o /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/

##  For paired-end data
### 2.2 import data to qiime2 and trim primers (if there are any)
qiime tools import \
            --type SampleData[PairedEndSequencesWithQuality] \
            --input-path /data/zool-temp-micro/pemb6000/16S/$STUDY/raw/ \
            --output-path /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/reads.qza \
            --input-format CasavaOneEightSingleLanePerSampleDirFmt

qiime demux summarize \
            --i-data /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/reads.qza \
            --o-visualization /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/reads_summary.qzv
#### trim primers using cutadapt, can also use for trimming to the same region such as V4, just need a bit adjustment of codes
qiime cutadapt trim-paired \
            --i-demultiplexed-sequences /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/reads.qza \
            --p-cores 12 \
            --p-front-f $primer_f \
            --p-front-r $primer_r \
            --p-discard-untrimmed \
            --p-no-indels \
            --o-trimmed-sequences /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/reads_trimmed.qza

qiime demux summarize \
            --i-data /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/reads_trimmed.qza \
            --o-visualization /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/reads_trimmed_summary.qzv

### 2.3 join pair-end reads using vsearch
if [ ! -d /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/deblur ]; then
        mkdir /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/deblur
fi

qiime vsearch join-pairs \
            --i-demultiplexed-seqs /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/reads_trimmed.qza \
            --o-joined-sequences /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/deblur/reads_trimmed_joined.qza \
	    --verbose \
	    --p-allowmergestagger
# summarize the joined reads            
qiime demux summarize \
            --i-data /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/deblur/reads_trimmed_joined.qza \
            --o-visualization /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/deblur/reads_trimmed_joined_summary.qzv

# filter out low-quality reads
qiime quality-filter q-score \
            --i-demux /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/deblur/reads_trimmed_joined.qza \
            --o-filter-stats /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/deblur/filt_stats.qza \
            --o-filtered-sequences /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/deblur/reads_trimmed_joined_filt.qza

# summarize, It is a good idea at this point just to verify that there haven't been any substantial losses of reads, before going through the whole ASV process, at either the joining or quality-filtering steps above:
qiime demux summarize \
            --i-data /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/deblur/reads_trimmed_joined_filt.qza \
            --o-visualization /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/deblur/reads_trimmed_joined_filt_summary.qzv



## For single or merged reads 
### 2.2 import data to qiime2 and trim primers (if there are any)
qiime tools import \
            --type SampleData[SequencesWithQuality] \
            --input-path /data/zool-temp-micro/pemb6000/16S/$STUDY/raw/ \
            --output-path /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/reads.qza \
            --input-format CasavaOneEightSingleLanePerSampleDirFmt
#fi

qiime demux summarize \
            --i-data /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/reads.qza \
            --o-visualization /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/reads_summary.qzv

#### 2.3 trim primers using cutadapt
qiime cutadapt trim-paired \
            --i-demultiplexed-sequences /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/reads.qza \
            --p-cores 12 \
            --p-front $primer_f \
            --p-adapter $adapter_f  \
            --p-discard-untrimmed \
            --p-no-indels \
            --o-trimmed-sequences /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/reads_trimmed.qza

qiime demux summarize \
            --i-data /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/reads_trimmed.qza \
            --o-visualization /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/reads_trimmed_summary.qzv

#### 2.4 here we don't need to join reads
if [ ! -d /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/deblur ]; then
        mkdir /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/deblur
fi

cp /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/reads.qza /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/deblur/reads_trimmed_joined.qza

#### filter out low-quality reads
qiime quality-filter q-score \
            --i-demux /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/deblur/reads_trimmed_joined.qza \
            --o-filter-stats /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/deblur/filt_stats.qza \
            --o-filtered-sequences /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/deblur/reads_trimmed_joined_filt.qza

# summarize, It is a good idea at this point just to verify that there haven't been any substantial losses of reads, before going through the whole ASV process, at either the joining or quality-filtering steps above:
qiime demux summarize \
            --i-data /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/deblur/reads_trimmed_joined_filt.qza \
            --o-visualization /data/zool-temp-micro/pemb6000/qiime2_result/$STUDY/deblur/reads_trimmed_joined_filt_summary.qzv


### Now we process merged reads, whether merged from pair-end or single
###trim and deblur, the trim_length is mannualy decided from visualizing reads_trimmed_joined_filt_summary.qzv, trim at the base when quality drops a lot but you don't want to lose too many reads
### 3. deblur for generating ASV tables
trim_length = "XXX"
qiime deblur denoise-16S \
            --i-demultiplexed-seqs /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/deblur/reads_trimmed_joined_filt.qza \
            --p-trim-length $trim_length \
            --p-left-trim-len 0 \
            --p-sample-stats \
            --p-jobs-to-start 16 \
            --p-min-reads 10 \
            --output-dir /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/deblur/denoise

# summarize
qiime deblur visualize-stats \
   --i-deblur-stats /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/deblur/denoise/stats.qza \
   --o-visualization /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/deblur/denoise/deblur-stats.qzv


qiime feature-table summarize \
   --i-table /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/deblur/denoise/table.qza \
   --o-visualization /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/deblur/denoise/deblur_table_summary.qzv

qiime feature-table tabulate-seqs \
   --i-data /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/deblur/denoise/representative_sequences.qza \
   --o-visualization /data/zool-temp-micro/pemb6000/chap1/qiime2_result/$STUDY/deblur/denoise/deblur_rep-seqs.qzv

######################################### Above is for within-study data processing ###############################
###### Now after processing dataset from individual study, we can merge all tables and sequences
###### below is an example of merging 4 datasets from 4 studies, but there can be more

qiime feature-table merge \
 --i-tables /data/zool-temp-micro/pemb6000/chap1/qiime2_result/Study1/deblur/denoise/table.qza \
 --i-tables /data/zool-temp-micro/pemb6000/chap1/qiime2_result/Study1/deblur/denoise/table.qza \
 --i-tables /data/zool-temp-micro/pemb6000/chap1/qiime2_result/Study3/deblur/denoise/table.qza \
 --i-tables /data/zool-temp-micro/pemb6000/chap1/qiime2_result/Study4/deblur/denoise/table.qza \
 --o-merged-table /data/zool-temp-micro/pemb6000/chap1/merged_table.qza


qiime feature-table merge-seqs \
 --i-data /data/zool-temp-micro/pemb6000/chap1/qiime2_result/Study1/deblur/denoise/representative_sequences.qza \
 --i-data /data/zool-temp-micro/pemb6000/chap1/qiime2_result/Study2/deblur/denoise/representative_sequences.qza \
 --i-data /data/zool-temp-micro/pemb6000/chap1/qiime2_result/Study3/deblur/denoise/representative_sequences.qza \
 --i-data /data/zool-temp-micro/pemb6000/chap1/qiime2_result/Study4/deblur/denoise/representative_sequences.qza \
 --o-merged-data /data/zool-temp-micro/pemb6000/chap1/merged_rep-seqs.qza


###### After merging, now we need to annotate the merged dataset

### 1. classify the features
qiime feature-classifier classify-sklearn \
            --i-reads /data/zool-temp-micro/pemb6000/chap1/merged_rep-seqs.qza \
            --i-classifier /data/zool-temp-micro/pemb6000/chap1/ref_alignments/silva-138-99-nb-classifier.qza \
            --p-n-jobs 16 \
            --output-dir /data/zool-temp-micro/pemb6000/chap1/taxonomy  &&

qiime metadata tabulate \
  --m-input-file /data/zool-temp-micro/pemb6000/chap1/taxonomy/classification.qza \
  --o-visualization /data/zool-temp-micro/pemb6000/chap1/taxonomy.qzv
#here can visualize rep-seq-dada2.qzv and choose ~5 ASVs to BLAST and to validate taxonomic classification.

### 2. filter out rare ASVs/potential contaminants and those not classified at kingdom level (or lower level)
qiime taxa filter-table \
            --i-table /data/zool-temp-micro/pemb6000/chap1/merged_table.qza \
            --i-taxonomy /data/zool-temp-micro/pemb6000/chap1/taxonomy/classification.qza \
            --p-include d__ \
            --p-exclude mitochondria,chloroplast,eukaryota \
            --o-filtered-table /data/zool-temp-micro/pemb6000/chap1/table_filt_decomtam.qza

qiime feature-table filter-seqs \
   --i-data /data/zool-temp-micro/pemb6000/chap1/merged_rep-seqs.qza \
   --i-table /data/zool-temp-micro/pemb6000/chap1/table_filt_decomtam.qza \
   --o-filtered-data /data/zool-temp-micro/pemb6000/chap1/merged_rep-seqs_final.qza


### 3. generate a phylogenetic tree
qiime fragment-insertion sepp \
   --i-representative-sequences /data/zool-temp-micro/pemb6000/chap1/merged_rep-seqs_final.qza \
   --i-reference-database /data/zool-temp-micro/pemb6000/chap1/sepp-refs-gg-13-8.qza \
   --o-tree /data/zool-temp-micro/pemb6000/chap1/asvs-tree.qza \
   --o-placements /data/zool-temp-micro/pemb6000/chap1/insertion-placements.qza \
   --p-threads 16

##export the tree
qiime tools export \
  --input-path /data/zool-temp-micro/pemb6000/chap1/asvs-tree.qza \
  --output-path /data/zool-temp-micro/pemb6000/chap1/saiga

