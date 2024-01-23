#!/usr/bin/env Rscript

idtaxa_func <- function(fasta,training_object_path , threshold , strands) {
  
  require(tidyverse)
  require(DECIPHER)
  
 

  if (!(threshold %in% c(40, 50, 60))) {
    cat("Be careful about threshold, default is 60% (very high confidence), 50% (high confidence), 40% (moderate confidence)\n")
  }
  
 

  fas<- fasta
  
  print(fas) # fasta file as input
  
  
  fileseqs <- readDNAStringSet(fas) # or readRNAStringSet
  
  # remove any gaps (if needed)
  fileseqs <- RemoveGaps(fileseqs)
  
  load(training_object_path)
  
  
  ids <- IdTaxa(fileseqs,
                trainingSet,
                strand=strands, # The default of "both" will classify using both orientations(reverse and forward on training seq) and choose the result with highest confidence.
                threshold= as.numeric(threshold), # 60 (cautious) or 50 (sensible)
                processors= NULL) # use all available processors with NULL
  
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
  
  
  taxid <- t(sapply(ids, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
  }))
  
  colnames(taxid) <- ranks
  
  taxid_table<-as.data.frame(taxid) 
  
  return(taxid_table)
  
}

args <- commandArgs(trailingOnly = TRUE)
fasta <- args[1]
training_object_path <- args[2]
threshold <- args[3]
strand <- args[4]

count_table <- idtaxa_func(fasta, training_object_path,threshold, strand)

write.table(count_table,paste0(str_remove(fasta,'.fasta'),'.count.table.csv'),quote = FALSE, row.names = TRUE,sep = '\t')







