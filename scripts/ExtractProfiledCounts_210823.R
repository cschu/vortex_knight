#!/usr/bin/env Rscript

.libPaths(c("/g/scb/zeller/fspringe/Software/R/4.1", .libPaths()))
source("/g/scb/zeller/fspringe/RScripts/functions/functions_read_in_profiled_data_210702.R")

library(optparse)
library(tidyverse)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 

option_list = list(
  make_option(c("--kraken2_res_path"), type="character", default=NULL, 
              help="Path to folder with kraken2 result tables", metavar="character"),
  make_option(c("--mOTUs_res_path"), type="character", default=NULL, 
              help="Path to folder with mOTUs result tables", metavar="character"),
  make_option(c("--PathSeq_res_path"), type="character", default=NULL, 
              help="Path to folder with PathSeq result tables", metavar="character"),
  make_option(c("--mTAGs_res_path"), type="character", default=NULL, 
              help="Path to folder with mTAGs result tables", metavar="character"),
  make_option(c("--mapseq_res_path"), type="character", default=NULL,
              help="Path to folder with mapseq result tables", metavar="character"),
  make_option(c("--libsize_res_path"), type="character", default=NULL, 
              help="Path to folder with libsize results", metavar="character"),
  make_option(c("--lib_layout_res_path"), type="character", default=NULL, 
              help="Path to folder with library_layout results", metavar="character"),
  make_option(c("--flagstats_res_path"), type="character", default=NULL, 
              help="Path to folder with flagstats", metavar="character"),
  make_option(c("--N_raw_counts_path"), type="character", default=NULL, 
              help="Path to folder with number of raw_counts", metavar="character"),
  make_option(c("--read_counter_res_path"), type="character", default=NULL, 
              help="Path to folder with read_counter results", metavar="character"),
  
  
  make_option(c("--out_folder"), type="character", default=NULL, 
              help="output folder path", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

out.folder <- opt$out_folder
message("kraken2")
#kraken2
if(!(is.null(opt$kraken2_res_path))){
  if(length(list.files(opt$kraken2_res_path))>0){
    res.kraken2 <- .f_read_in_files_kraken2(path_to_folder = opt$kraken2_res_path,
                                            tax.level = "genus")
    saveRDS(res.kraken2,paste0(out.folder,"/res_kraken2.rds"))  
  }else{
    message("kraken2 path empty")
  }
}
message("motus")
#mOTUs
if(!(is.null(opt$mOTUs_res_path))){
  if(length(list.files(opt$mOTUs_res_path))>0){
    res.mOTUs <- .f_read_in_files_mOTUs(path_to_folder = opt$mOTUs_res_path)
    saveRDS(res.mOTUs,paste0(out.folder,"/res_mOTUs.rds"))
  } else{
    message("mOTUs path empty")
  }
}
#mTAGs
message("mtags")
if(!(is.null(opt$mTAGs_res_path))){
  if(length(list.files(opt$mTAGs_res_path))>0){
  res.mTAGs <- .f_read_in_files_mTAGsProfile_merged(path_to_folder = opt$mTAGs_res_path,
                                                    tax.level = "genus")
  saveRDS(res.mTAGs,paste0(out.folder,"/res_mTAGs.rds"))
  } else{
    message("mTAGs path empty")
  }
}

#mapseq
message("mapseq")
if(!(is.null(opt$mapseq_res_path))){
  if(length(list.files(opt$mapseq_res_path))>0){
  res.mapseq <- .f_read_in_mapseq_raw(path_to_folder = opt$mapseq_res_path,tax.lvl = "genus")
  saveRDS(res.mapseq,paste0(out.folder,"/res_mapseq.rds"))
  } else{
    message("mapseq path empty")
  }
}

#pathseq
message("pathseq")
if(!(is.null(opt$PathSeq_res_path))){
  if(length(list.files(opt$PathSeq_res_path))>0){
  res.pathseq <- .f_read_in_files_PathSeq(path_to_folder = opt$PathSeq_res_path,
                                                    tax.level = "genus")
  saveRDS(res.pathseq,paste0(out.folder,"/res_pathseq.rds"))
  } else{
    message("PathSeq path empty")
  }
}

#libsize
if(!(is.null(opt$libsize_res_path))){
  if(length(list.files(opt$libsize_res_path))>0){
    res.libsize <- .f_read_in_libsize(path_to_folder = opt$libsize_res_path)
    saveRDS(res.libsize,paste0(out.folder,"/res_libsize.rds"))
  } else{
    message("libsize path empty")
  }
}

#liblayout
if(!(is.null(opt$lib_layout_res_path))){
  if(length(list.files(opt$lib_layout_res_path))>0){
    res.lib_layout <- .f_read_in_lib_layout(path_to_folder = opt$lib_layout_res_path)
    saveRDS(res.lib_layout,paste0(out.folder,"/res_lib_layout.rds"))
  } else{
    message("liblayout path empty")
  }
}

#flagstats --> total numbers of reads before QC
#No needed anymore in vknight build >= 827290457b
if(!(is.null(opt$flagstats_res_path))){
  if(length(list.files(opt$flagstats_res_path))>0){
    res.flagstats <- .f_read_in_flagstats(path_to_folder = opt$flagstats_res_path)
    saveRDS(res.flagstats,paste0(out.folder,"/res_flagstat_total_reads.rds"))
  } else{
    message("flagstats path empty")
  }
}

#raw_counts folder: total number of raw_counts (before QC)
if(!(is.null(opt$N_raw_counts_path))){
  if(file.exists(opt$N_raw_counts_path)){
    if(length(list.files(opt$N_raw_counts_pat))>0){ #important for vknight implementation
      res.N_raw_counts <- .f_read_in_raw_counts_number(path_to_folder = opt$N_raw_counts_path)
      saveRDS(res.N_raw_counts,paste0(out.folder,"/res_libsize_raw_before_QC.rds"))
    } else{
      message("No raw_counts folder present - No QC performed in the current run")
    }
  }
}

message("read_counter")
#read_counter folder (Alessio's tool that is the basis of mOTUs)
if(!(is.null(opt$read_counter_res_path))){
  if(length(list.files(opt$read_counter_res_path))>0){
    #Load dataframe with information on the number of marker genes per species
    marker_gene_lengths.df <- read_tsv("/g/scb/zeller/fspringe/Database/GTDB/GTDB_marker_gene_lengths.tsv")
    min_genes_vec <- seq(2,20)
    res.rc.list <- list()
    for(i in seq(1,length(min_genes_vec))){
      message(paste0("Collating read_counter results with g= ",min_genes_vec[i]))
      tmp.res <- .f_read_in_read_counter(path_to_folder = opt$read_counter_res_path,
                                         min_genes = min_genes_vec[i],
                                         marker_genes.df = marker_gene_lengths.df, 
                                         tax.level = "genus",
                                         norm_to_gene_length = TRUE)
      res.rc.list[[i]] <- tmp.res
    }
    names(res.rc.list) <- paste0("g",min_genes_vec)
    #save res.rc.list
    saveRDS(res.rc.list,paste0(out.folder,"/res_read_counter_norm.rds"))
  } else{
    message("read_counter path empty")
  }
}






