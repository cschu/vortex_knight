#!/usr/bin/env Rscript

# .libPaths(c("/g/scb/zeller/fspringe/Software/R/4.1", .libPaths()))
# source("/g/scb/zeller/fspringe/RScripts/functions/functions_read_in_profiled_data_210702.R")
# library(here)
# source(here::here("functions_read_in_profiled_data_210702.R"))
source("/scratch/schudoma/vknight/scripts/functions_read_in_profiled_data_210702.R")

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
  make_option(c("--libsize_res_path"), type="character", default=NULL, 
              help="Path to folder with libsize results", metavar="character"),
  make_option(c("--lib_layout_res_path"), type="character", default=NULL, 
              help="Path to folder with library_layout results", metavar="character"),
  
  
  make_option(c("--out_folder"), type="character", default=NULL, 
              help="output folder path", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

out.folder <- opt$out_folder

#kraken2
if(!(is.null(opt$kraken2_res_path))){
  res.kraken2 <- .f_read_in_files_kraken2(path_to_folder = opt$kraken2_res_path,
                                          tax.level = "genus")
  saveRDS(res.kraken2,paste0(out.folder,"/res_kraken2.rds"))
}
#mOTUs
if(!(is.null(opt$mOTUs_res_path))){
  res.mOTUs <- .f_read_in_files_mOTUs(path_to_folder = opt$mOTUs_res_path)
  saveRDS(res.mOTUs,paste0(out.folder,"/res_mOTUs.rds"))
}
#mTAGs
if(!(is.null(opt$mTAGs_res_path))){
  res.mTAGs <- .f_read_in_files_mTAGsProfile_merged(path_to_folder = opt$mTAGs_res_path,
                                                    tax.level = "genus")
  saveRDS(res.mTAGs,paste0(out.folder,"/res_mTAGs.rds"))
}
#pathseq
if(!(is.null(opt$PathSeq_res_path))){
  res.pathseq <- .f_read_in_files_PathSeq(path_to_folder = opt$PathSeq_res_path,
                                                    tax.level = "genus")
  saveRDS(res.pathseq,paste0(out.folder,"/res_pathseq.rds"))
  
}

#libsize
if(!(is.null(opt$libsize_res_path))){
  res.libsize <- .f_read_in_libsize(path_to_folder = opt$libsize_res_path)
  saveRDS(res.libsize,paste0(out.folder,"/res_libsize.rds"))
}

#liblayout
if(!(is.null(opt$lib_layout_res_path))){
  res.lib_layout <- .f_read_in_lib_layout(path_to_folder = opt$lib_layout_res_path)
  saveRDS(res.lib_layout,paste0(out.folder,"/res_lib_layout.rds"))
}








