#!/bin/bash
#SBATCH -A zeller
#SBATCH -t 24:00:00
#SBATCH -n 1
#SBATCH --mem 12G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o ./collateResults.out
#SBATCH -e ./collateResults.err

#Define nextflow project folder
Project_Folder=$(pwd)

#define vknight --output_dir
out_dir=Results/

#create temporary folder to Collate symlinks of all files from the same tool in a folder per tool
mkdir -p ${Project_Folder}/tmp_Results

#Create folder to store rds files of the profiler outputs and QC information to download
mkdir -p ${Project_Folder}/vknight_results_collated


#### Collate the individual result files ####
#Collate kraken2 results
mkdir -p ${Project_Folder}/tmp_Results/kraken2
find ${Project_Folder}/${out_dir}/ -name "*kraken2_report.txt" -exec ln -sf {} ${Project_Folder}/tmp_Results/kraken2/ \;

#Collate PathSeq results
mkdir -p ${Project_Folder}/tmp_Results/PathSeq
find ${Project_Folder}/${out_dir}/ -name "*pathseq.scores" -exec ln -sf {} ${Project_Folder}/tmp_Results/PathSeq/ \;

#Collate mOTUs results
mkdir -p ${Project_Folder}/tmp_Results/mOTUs
find ${Project_Folder}/${out_dir}/ -name "*motus.txt" -exec ln -sf {} ${Project_Folder}/tmp_Results/mOTUs/ \;

#Collate libsize
mkdir -p ${Project_Folder}/tmp_Results/libsize
find ${Project_Folder}/${out_dir}/ -name "*libsize.txt" -exec ln -sf {} ${Project_Folder}/tmp_Results/libsize/ \;

#Collate libary layout
mkdir -p ${Project_Folder}/tmp_Results/lib_layout
find ${Project_Folder}/${out_dir}/ -name "*is_paired.txt" -exec ln -sf {} ${Project_Folder}/tmp_Results/lib_layout/ \;

#Collate flagstats results
mkdir -p ${Project_Folder}/tmp_Results/flagstats
find ${Project_Folder}/${out_dir}/ -name "*flagstats.txt" -exec ln -sf {} ${Project_Folder}/tmp_Results/flagstats/ \;

#Collate OTU tables
#mkdir -p ${Project_Folder}/tmp_Results/otu_tables
#find ${Project_Folder}/Results/ -name "*bac_ssu.tsv" -exec ln -sf {} ${Project_Folder}/tmp_Results/otu_tables/ \;

#manually collate and assign OTU tables from .mseq results
mkdir -p ${Project_Folder}/tmp_Results/mapseq
find ${Project_Folder}/work/ -type f -name "*.mseq" -exec ln -sf {} ${Project_Folder}/tmp_Results/mapseq/ \;

#Collate mtags tables
mkdir -p ${Project_Folder}/tmp_Results/mtags_tables
find ${Project_Folder}/${out_dir}/ -name "merged_profile.genus.tsv" -exec ln -sf {} ${Project_Folder}/tmp_Results/mtags_tables/ \;

#Collate read_counter results
mkdir -p ${Project_Folder}/tmp_Results/read_counter
find ${Project_Folder}/${out_dir}/ -name "*read_counter.txt" -exec ln -sf {} ${Project_Folder}/tmp_Results/read_counter/ \;

#Collate the fasta files with extracted ribosomal sequences from mtags_extract
mkdir -p ${Project_Folder}/tmp_Results/mtags_extract_fastq
find ${Project_Folder}/work -name "*bac_ssu.fasta" -exec ln -sf {} ${Project_Folder}/tmp_Results/mtags_extract_fastq/ \;


#Copy the multiqc report to the vknight_results_folder
find ${Project_Folder}/${out_dir}/ -name "multiqc_report.html" -exec cp {} ${Project_Folder}/vknight_results_collated/ \;

#If QC is performed, also create symplings for the QCed FastQ files
mkdir -p ${Project_Folder}/tmp_Results/QCed_FastQ
for f in $(find ${Project_Folder}/work/ -type f -name '*bbduk_stats.txt'); 
do d=$(dirname $f); ln -s $d/*R?.fastq.gz ${Project_Folder}/tmp_Results/QCed_FastQ/; done

#Copy the files for run overview and resource usage into the results folder
find ${Project_Folder} -name "timeline*" -type f -exec cp {} ${Project_Folder}/vknight_results_collated/ \;
find ${Project_Folder} -name "trace*" -type f -exec cp {} ${Project_Folder}/vknight_results_collated/ \;
find ${Project_Folder} -name "report*" -type f -exec cp {} ${Project_Folder}/vknight_results_collated/ \;

#### Compute rds files with count matrices / libsize / lib_layout etc ####
module load R/4.1.0-foss-2021a
Rscript /g/scb/zeller/fspringe/RScripts/ExtractProfiledCounts_210823.R \
	--kraken2_res_path ${Project_Folder}/tmp_Results/kraken2/ \
	--mOTUs_res_path ${Project_Folder}/tmp_Results/mOTUs/ \
	--PathSeq_res_path ${Project_Folder}/tmp_Results/PathSeq/ \
	--mTAGs_res_path ${Project_Folder}/tmp_Results/mtags_tables/ \
	--mapseq_res_path ${Project_Folder}/tmp_Results/mapseq/ \
	--libsize_res_path ${Project_Folder}/tmp_Results/libsize/ \
	--lib_layout_res_path ${Project_Folder}/tmp_Results/lib_layout/ \
	--N_raw_counts_path ${Project_Folder}/Results/raw_counts/ \
	--read_counter_res_path ${Project_Folder}/tmp_Results/read_counter/ \
	--out_folder ${Project_Folder}/vknight_results_collated

#### Send as zipped email to user ####
#tar czf ${Project_Folder}/vknight_results_collated.tar.gz ${Project_Folder}/vknight_results_collated
#echo 'Hi, here is the results folder!' | mailx -s vknight-results -a ${Project_Folder}/vknight_results_collated.tar.gz fabian.springer@embl.de




