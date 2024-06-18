#library(stringr)
library(progress)

.f_read_in_files_mapseq <- function(path_to_folder){
  #Deprecated since vknight uses by default mapseq with --outfmt simepl (-> read in funciton for that is called ".f_read_in_mapseq_raw")
  #Attention: This function refers to the "default" mapseq output and not the mapseq output generated via --outfmt simple
  # reads in mapseq output file; attention: in mapseq output file, the first two lines are commented out (with a "#"); 
  # The first line (showing the processing data and mapseq version) will be ignored
  # The second line, giving the total sum of counts, will be added to the output matrix
  # If there are mapseq output files for forward and reverse reads, the counts will be combined
  
  #consider only counts from bacterial small subunit (bac_ssu)
  file_list <- list.files(path=path_to_folder,pattern = "*_bac_ssu.tsv")
  
  #if forward and reverse files are present: merge them and sum up counts
  count.df <- tibble()
  if(length(file_list)>0){
    for(i in seq(1,length(file_list))){
      if(file.info(paste0(path_to_folder,file_list[i]))$size>0){
        tmp.file <- read.table(paste0(path_to_folder,file_list[i]),sep = "\t",comment.char = "",skip = 1)
        colnames(tmp.file)[-1] <- tmp.file[1,-1]
        tmp.file <- tmp.file[-1,]
        #colnames(tmp.file) <- gsub(colnames(tmp.file), pattern = "*_bac_ssu.mseq",replacement = "")
        #colnames(tmp.file) <- str_extract(colnames(tmp.file), "[^_]+")
        tmp.file <- tmp.file %>% mutate_at(vars(-V1),as.numeric)
        
        count.df <- bind_rows(count.df,tmp.file)
      }else{
        message(paste0(file_list[i]," empty"))
        next
      }
    }
    if(nrow(count.df)>0){
      count.df <- as.data.frame(count.df %>% group_by(V1) %>% summarise_all(sum))
      #create output matrix (assign rownames and remove column with tax names)
      rownames(count.df) <- sub('.*\\;', '', count.df$V1)
      count.df$V1 <- NULL
      count.mat <- as.matrix(count.df)
      rownames(count.mat)[1] <- "Bacteria"
      #remove everything after the first point in the samplenames
      #colnames(count.mat) <- sub(colnames(count.mat),pattern = "\\..*",replacement = "")
      return(count.mat)
    }else{
      return(as.matrix(count.df))
    }
  }else{
    stop("no files found in the given directory")
  }
}

.f_metaphlanize_mapseq <- function(c.combined_df){
  #Helper function for read_in_mapseq_raw
  #takes mapseq output and creates a dataframe with the sum of counts at every taxonomic level
  #also adds the prefixes to the individual tax levels
  
  taxLevel_short <- c("k__","p__","c__","o__","f__","g__","s__")
  taxLevels <- c("kingdom","phylum","class","order","family","genus","species")
  res_df <- tibble(tax=character(0),
                   count=double(0))
  #step 1 create df with only the tax levels and the added prefixes
  tmp <- c.combined_df %>% 
    select(-tax) %>% 
    distinct() %>% 
    as.data.frame()
  j <- 2
  for(j in seq(2,ncol(tmp))){
    tmp[,j] <- paste0(taxLevel_short[j-1],tmp[,j])
  }
  taxReference_df <- 
    tmp %>% 
    mutate(across(.cols = -count,.fns = ~map(., str_replace_all,pattern = ".__NA",replacement = NA_character_))) %>% 
    #unite("tax",kingdom:species,sep = "|",na.rm = T,remove = F) %>% 
    unnest(cols=c(kingdom, phylum, class, order, family, genus, species)) %>% 
    #relocate(tax) %>% 
    as_tibble()
  
  
  #Step2 now do the grouping individually and create res_df with total counts at the given tax level
  res_df <- tibble(tax = character(0),
                   count = double(0))
  for(i in seq(1,length(taxLevels))){
    taxLevel_df <- 
      taxReference_df %>% 
      select(count,!!as.symbol(taxLevels[1]):!!as.symbol(taxLevels[i])) %>% #select increasing number of tax levels
      filter(!is.na(!!as.symbol(taxLevels[i]))) %>% #remove NAs in current tax level (reads are not resolved there)
      unite("tax",-count,sep = "|") %>% 
      group_by(tax) %>% 
      summarise(count = sum(count))
    res_df <- bind_rows(res_df,taxLevel_df)
  }
  return(res_df)
}
.f_correctNonMatchingTaxAssignments_mapseq <- function(mat1,mat2){
  # case1: taxonomic assignment of one read is "deeper" than the other (but all levels are matching):
  # Take deeper annotation
  
  # case2: taxonomic assignment differs between fwd and reverse reads at some level:
  # Take last common ancestor
  
  #this function was mostly written by ChatGPT and significantly increases performance compared to the old version (pre-2023-02-07)
  n_rows <- nrow(mat1)
  consensus_matrix <- matrix(NA, nrow = n_rows, ncol = 7)
  
  pb_1 <- progress_bar$new(total=nrow(consensus_matrix))
  for (i in 1:n_rows) {
    taxonomy1 <- as.character(mat1[i,])
    taxonomy2 <- as.character(mat2[i,])
    j <- 3
    #message(i)
    for (j in 1:7) {
      if (isTRUE(taxonomy1[j] == taxonomy2[j])) {
        consensus_matrix[i, j] <- taxonomy1[j]
      } else if (is.na(taxonomy1[j] == taxonomy2[j])){#if one of the entries is NA, fill with the other entry
        consensus_matrix[i, j] <- ifelse(is.na(taxonomy1[j]), taxonomy2[j], taxonomy1[j])
      } else{
        break #if there os no match and also no NA in one of the two: break
      }
    }
    pb_1$tick()
  }
  #remove rows with only NA values (when fwd and rev read have entirely different taxonomies)
  rownames(consensus_matrix) <- rownames(mat1)
  colnames(consensus_matrix) <- colnames(mat1)
  consensus_matrix <- consensus_matrix[rowSums(is.na(consensus_matrix))<ncol(consensus_matrix),,drop=F]
  return(consensus_matrix)
}

.f_read_in_mapseq_raw <- function(path_to_folder,fname_mode="vknight",
                                  account_for_paired_reads=NULL,
                                  output_style="standard",
                                  tax.lvl=NULL){
  #account for paired reads and tax.lvl not necessary anymore but left optional in order to not break scripts etc
  account_for_paired_reads <- NULL
  tax.lvl <- NULL
  
  
  #fname_mode defines how the filenames were constructed: 
  #vknight: SAMPLE_ID_R[0-9]_bac_ssu.mseq
  #manual: SAMPLE_ID_1.mseq
  #output_style: "standard" --> just export the standardized mapseq output (summarized counts for every taxonomic level)
  #              "metaphlan"--> metaphlanized output version
  #!!! Works only if .mseq file was computed with --outfmt simple !!!#
  
  file_list <- list.files(path_to_folder,pattern = ".mseq")
  if(length(file_list)<1){message("No .mseq files in given directory")
    stop
  }
  
  #Get the filenames that are expected in the .mseq folder
  if(fname_mode == "vknight"){
    sample.names <- unique(str_remove(file_list,pattern = "_R[0-9]_bac_ssu.mseq"))
    fwd_file_ending <- "_R1_bac_ssu.mseq"
    rev_file_ending <- "_R2_bac_ssu.mseq"
  }else if(fname_mode == "manual"){
    sample.names <- unique(str_remove(file_list,pattern = "_[0-9].mseq"))
    fwd_file_ending <- "_1.mseq"
    rev_file_ending <- "_2.mseq"
    # }else if(fname_mode == "dada2"){
    #   sample.names <- unique(str_remove(file_list,pattern = ".mseq"))
    #   fwd_file_ending <- ".mseq"
    #   #rev_file_ending <- "_2.mseq"
    # }
  }else{
    message("Incorrect fname_mode")
    stop()
  }
  message(paste0("Importing ",length(sample.names)," samples"))
  count_df <- tibble(tax = character(0))
  pb <- progress_bar$new(total=length(sample.names))
  i <- 1
  for(i in seq(1,length(sample.names))){
    #message(i)
    c.sample <- sample.names[i]
    
    #load forward and reverse reads
    c.fwd <- tryCatch(
      {data.table::fread(file = paste0(path_to_folder,c.sample,fwd_file_ending),skip = 1,header = T)},
      error=function(e){
        c.fwd <- data.frame(matrix(ncol = 15,nrow = 0))
      }
    )
    c.rev <- tryCatch(
      {data.table::fread(file = paste0(path_to_folder,c.sample,rev_file_ending),skip = 1,header = T)},
      error=function(e){
        c.rev <- data.frame(matrix(ncol = 15,nrow = 0))
      }
    )
    colnames(c.fwd) <- as.character(seq(1,ncol(c.fwd)))
    colnames(c.rev) <- as.character(seq(1,ncol(c.rev)))
    
    #check if both, forward and reverse reads are present. If not, treat sample as unpaired
    if(nrow(c.fwd)==0 & nrow(c.rev)==0){#skip if both read files are empty
      next
      pb$tick()
    } 
    if(nrow(c.fwd)>0 & nrow(c.rev)>0){
      sample.type <- "paired"
    }else{
      sample.type <- "single"
    }
    
    #Compare taxonomic assignments of forward and reverse reads
    c.combined <- rbind(c.fwd[,c(1,14)] %>% add_column(read = "fwd"),
                        c.rev[,c(1,14)] %>% add_column(read = "rev"))
    colnames(c.combined) <- c("read_id","tax","read")
    
    ### split the taxonomy, add the taxon-prefixes and re-merge it
    tax.split <- data.table::tstrsplit(x = c.combined$tax,";")
    #make list with 7 entries to match taxonomic tree
    if(length(tax.split)<7){
      N <- 7-length(tax.split)
      NA_list <- vector("list", N)
      # Populate the list with vectors of NAs (must be NA characters since otherwise they would be converted to logicals by data.table)
      for (z in 1:N) {
        NA_list[[z]] <- rep(as.character(NA), nrow(c.combined))
      }
      tax.split <- c(tax.split,NA_list)
    }
    c.combined[,c('kingdom', 'phylum', 'class','order', 'family', 'genus', 'species')] <- tax.split
    
    # If sample is paired-end do assign the taxonomie based on the agreement between read pairs
    if(sample.type == "paired"){
      fwd_mat <- c.combined %>% filter(read == "fwd") %>% select(-read) %>% column_to_rownames("read_id")
      rev_mat <- c.combined %>% filter(read == "rev") %>% select(-read) %>% column_to_rownames("read_id")
      
      # Get intersection of read names and compare assigned taxonomies
      readNames_intersect <- intersect(rownames(fwd_mat),rownames(rev_mat))
      fwd_readsIntersect_mat <- fwd_mat[readNames_intersect,]
      rev_readsIntersect_mat <- rev_mat[readNames_intersect,]
      
      # for the taxonomies that match: keep them as is
      matching_reads_mat <- fwd_readsIntersect_mat[fwd_readsIntersect_mat$tax == rev_readsIntersect_mat$tax,
                                                   colnames(fwd_readsIntersect_mat)!="tax"]
      # If there are reads with an assigned taxonomy that are only found in either the fwd or rev: add them to the matching_mat (since this is the only information there is)
      orphans_mat <- rbind(fwd_mat[!(rownames(fwd_mat)%in%readNames_intersect),],
                           rev_mat[!(rownames(rev_mat)%in%readNames_intersect),])
      matching_reads_mat <- rbind(matching_reads_mat,orphans_mat[,colnames(orphans_mat)!="tax"])
      
      #get some statistics about matching/non-matching
      nMatching <- nrow(matching_reads_mat)
      nNonMatching <- length(unique(c.combined$read_id))-nMatching
      fracNonMatching <- round(nNonMatching / (nMatching+nNonMatching) * 100,0)
      message("\n",nNonMatching," reads (",fracNonMatching,"%) have non-identical taxonomic annotation")
      
      # for the reads that don't have matching taxonomy annotations: process them accordingly
      fwdNonMatching_mat <- fwd_readsIntersect_mat[fwd_readsIntersect_mat$tax != rev_readsIntersect_mat$tax,
                                                   colnames(fwd_readsIntersect_mat)!="tax"]
      revNonMatching_mat <- rev_readsIntersect_mat[fwd_readsIntersect_mat$tax != rev_readsIntersect_mat$tax,
                                                   colnames(rev_readsIntersect_mat)!="tax"]
      # If there are non matching taxonomies between fwd and rev reads: harmonize them
      if(nrow(fwdNonMatching_mat)>0){
        nonMatching_reads_fixed_mat <- .f_correctNonMatchingTaxAssignments_mapseq(mat1 = fwdNonMatching_mat,
                                                                                  mat2 = revNonMatching_mat)
      }else{nonMatching_reads_fixed_mat <- data.frame(matrix(NA,nrow = 0,ncol = ncol(matching_reads_mat),dimnames = list(NULL,colnames(matching_reads_mat))))}
      
      #combine the read matrices and proceed and unite the taxonomy in order to do an initial grouping before assigning the tax-prefixes (saves time)
      combined_reads_df <- rbind(matching_reads_mat,
                                 nonMatching_reads_fixed_mat) %>% 
        rownames_to_column("read_id") %>% 
        as_tibble() %>% 
        unite(tmpTax,kingdom:species,remove = F,na.rm = T,sep = "|") %>% 
        group_by(tmpTax) %>% 
        mutate(count = n()) %>% 
        ungroup() %>% 
        select(-read_id) %>% 
        distinct() %>% 
        relocate(tmpTax,count)
    }else{#if sample is single-end, simply take the provided taxonomic annotation of the read (do a grouping and summary based on the taxonomy in order to save computing time)
      combined_reads_df <- 
        c.combined %>% 
        as_tibble() %>% 
        group_by(tax) %>% 
        mutate(count = n()) %>% 
        ungroup() %>% 
        select(-read_id,-read) %>% 
        distinct() %>% 
        rename(tmpTax = tax) %>% 
        relocate(tmpTax,count)
      
    }
    
    # Add taxonomic prefixes
    taxLevel_short <- c("k__","p__","c__","o__","f__","g__","s__")
    taxLevels <- c("kingdom","phylum","class","order","family","genus","species")
    res_df <- tibble(tax=character(0),
                     count=double(0))
    #create df with only the tax levels and the added prefixes
    tmp <- combined_reads_df %>% 
      select(-tmpTax) %>% 
      distinct() %>% 
      as.data.frame()
    for(j in seq(2,ncol(tmp))){
      tmp[,j] <- paste0(taxLevel_short[j-1],tmp[,j])
    }
    taxReference_df <- 
      tmp %>% 
      mutate(across(.cols = -count,.fns = ~map(., str_replace_all,pattern = ".__NA",replacement = NA_character_))) %>% 
      unnest(cols=c(kingdom, phylum, class, order, family, genus, species)) %>% 
      as_tibble()
    
    #Convert to the desired output style matrix
    if(output_style == "metaphlan"){
      # Metaphlanize the mapseq output (give total sum of counts at every taxonomic combination)
      res_df <- .f_metaphlanize_mapseq(taxReference_df = taxReference_df)
      
    }else if(output_style == "standard"){
      res_df <- 
        taxReference_df %>% 
        unite("tax",kingdom:species,sep = "|",na.rm = T,remove = F) %>% 
        select(tax,count)
    }
    
    
    ### merge with count_df
    count_df <- 
      count_df %>% 
      full_join(.,res_df %>% rename(!!as.symbol(c.sample) := count),
                by="tax")
    
    #remove c.combined (to prevent for whatever reason that a failed sample leads to duplication of the previous one)
    rm(res_df,c.combined,taxReference_df,tmp)
    pb$tick()
  }
  ### Create output matrix
  count_mat <- count_df %>% 
    as_tibble() %>% 
    column_to_rownames("tax") %>% 
    replace(is.na(.),0) %>% 
    as.matrix()
  
  return(count_mat)
}
i
.f_read_in_IDtaxa <- function(path_to_folder,output_style="standard") {
  # Collate all IDtaxa files into one count matrix with full taxonomy as rownmes.
  # Follow the same strategy for building consensus taxonomy from rev and fwd read as with MAPseq: 
  # For every read pair, check: 
  # If taxonomy is identical (same resolution, same assignment) -> Take as is
  # If taxonomy diverges at some point (e.g. different species assignment) -> Take last common ancestor
  # If taxonomy is identical but different resolution (e.g. fwd until species, rev only genus) -> take higher resolution


  file_list <- list.files(path_to_folder, pattern = "*IDTaxa.tsv")
  if (length(file_list) < 1) {
    message("No IDTaxa.tsv files in given directory")
    stop
  }
  sample.names <- unique(str_remove(file_list,pattern = "_R[12]_IDTaxa.tsv"))
  
  #* Import fwd and rev file and combine them ----
  # Due to laziness this is copied from the MAPseq import function - sorry for the bad coding practice in the old code.
  message(paste0("Importing ",length(sample.names)," samples"))
  count_df <- tibble(tax = character(0))
  pb <- progress_bar$new(total=length(sample.names))
  i <- 30
  for(i in seq(1,length(sample.names))){
    #message(i)
    c.sample <- sample.names[i]
    
    #load forward and reverse reads
    c.fwd <- tryCatch(
      {data.table::fread(file = paste0(path_to_folder,c.sample,"_R1_IDTaxa.tsv"),sep="\t",skip = 0,header = F)},
      error=function(e){
        c.fwd <- data.frame(matrix(ncol = 3,nrow = 0)) #IDTaxa raw output has 3 columns
      }
    )
    c.rev <- tryCatch(
      {data.table::fread(file = paste0(path_to_folder,c.sample,"_R2_IDTaxa.tsv"),sep="\t",skip = 0,header = F)},
      error=function(e){
        c.rev <- data.frame(matrix(ncol = 3,nrow = 0))
      }
    )
    colnames(c.fwd) <- as.character(seq(1,ncol(c.fwd)))
    colnames(c.rev) <- as.character(seq(1,ncol(c.rev)))
    
    #check if both, forward and reverse reads are present. If not, treat sample as unpaired
    if(nrow(c.fwd)==0 & nrow(c.rev)==0){#skip if both read files are empty
      next
      pb$tick()
    } 
    if(nrow(c.fwd)>0 & nrow(c.rev)>0){
      sample.type <- "paired"
    }else{
      sample.type <- "single"
    }
    c.fwd %>% as_tibble()
    #Compare taxonomic assignments of forward and reverse reads
    c.combined <- rbind(c.fwd[,c(2,3)] %>% add_column(read = "fwd"),
                        c.rev[,c(2,3)] %>% add_column(read = "rev"))
    colnames(c.combined) <- c("read_id","tax","read")

    #* Part II: split the taxonomy, add the taxon-prefixes and re-merge it ----
    tax.split <- data.table::tstrsplit(x = c.combined$tax,";")
    #make list with 7 entries to match taxonomic tree
    if(length(tax.split)<7){
      N <- 7-length(tax.split)
      NA_list <- vector("list", N)
      # Populate the list with vectors of NAs (must be NA characters since otherwise they would be converted to logicals by data.table)
      for (z in 1:N) {
        NA_list[[z]] <- rep(as.character(NA), nrow(c.combined))
      }
      tax.split <- c(tax.split,NA_list)
    }
    c.combined[,c('kingdom', 'phylum', 'class','order', 'family', 'genus', 'species')] <- tax.split
    
    # If sample is paired-end do assign the taxonomie based on the agreement between read pairs
    if(sample.type == "paired"){
      fwd_mat <- c.combined %>% filter(read == "fwd") %>% select(-read) %>% column_to_rownames("read_id")
      rev_mat <- c.combined %>% filter(read == "rev") %>% select(-read) %>% column_to_rownames("read_id")
      
      # Get intersection of read names and compare assigned taxonomies
      readNames_intersect <- intersect(rownames(fwd_mat),rownames(rev_mat))
      fwd_readsIntersect_mat <- fwd_mat[readNames_intersect,]
      rev_readsIntersect_mat <- rev_mat[readNames_intersect,]
      
      # for the taxonomies that match: keep them as is
      matching_reads_mat <- fwd_readsIntersect_mat[fwd_readsIntersect_mat$tax == rev_readsIntersect_mat$tax,
                                                   colnames(fwd_readsIntersect_mat)!="tax"]
      # If there are reads with an assigned taxonomy that are only found in either the fwd or rev: add them to the matching_mat (since this is the only information there is)
      orphans_mat <- rbind(fwd_mat[!(rownames(fwd_mat)%in%readNames_intersect),],
                           rev_mat[!(rownames(rev_mat)%in%readNames_intersect),])
      matching_reads_mat <- rbind(matching_reads_mat,orphans_mat[,colnames(orphans_mat)!="tax"])
      
      #get some statistics about matching/non-matching
      nMatching <- nrow(matching_reads_mat)
      nNonMatching <- length(unique(c.combined$read_id))-nMatching
      fracNonMatching <- round(nNonMatching / (nMatching+nNonMatching) * 100,0)
      message("\n",nNonMatching," reads (",fracNonMatching,"%) have non-identical taxonomic annotation")
      
      # for the reads that don't have matching taxonomy annotations: process them accordingly
      fwdNonMatching_mat <- fwd_readsIntersect_mat[fwd_readsIntersect_mat$tax != rev_readsIntersect_mat$tax,
                                                   colnames(fwd_readsIntersect_mat)!="tax"]
      revNonMatching_mat <- rev_readsIntersect_mat[fwd_readsIntersect_mat$tax != rev_readsIntersect_mat$tax,
                                                   colnames(rev_readsIntersect_mat)!="tax"]
      # If there are non matching taxonomies between fwd and rev reads: harmonize them
      if(nrow(fwdNonMatching_mat)>0){
        nonMatching_reads_fixed_mat <- .f_correctNonMatchingTaxAssignments_mapseq(mat1 = fwdNonMatching_mat,
                                                                                  mat2 = revNonMatching_mat)
      }else{nonMatching_reads_fixed_mat <- data.frame(matrix(NA,nrow = 0,ncol = ncol(matching_reads_mat),dimnames = list(NULL,colnames(matching_reads_mat))))}
      
      #combine the read matrices and proceed and unite the taxonomy in order to do an initial grouping before assigning the tax-prefixes (saves time)
      combined_reads_df <- rbind(matching_reads_mat,
                                 nonMatching_reads_fixed_mat) %>% 
        rownames_to_column("read_id") %>% 
        as_tibble() %>% 
        unite(tmpTax,kingdom:species,remove = F,na.rm = T,sep = "|") %>% 
        group_by(tmpTax) %>% 
        mutate(count = n()) %>% 
        ungroup() %>% 
        select(-read_id) %>% 
        distinct() %>% 
        relocate(tmpTax,count)
    }else{#if sample is single-end, simply take the provided taxonomic annotation of the read (do a grouping and summary based on the taxonomy in order to save computing time)
      combined_reads_df <- 
        c.combined %>% 
        as_tibble() %>% 
        group_by(tax) %>% 
        mutate(count = n()) %>% 
        ungroup() %>% 
        select(-read_id,-read) %>% 
        distinct() %>% 
        rename(tmpTax = tax) %>% 
        relocate(tmpTax,count)
      
    }
    #* Part III: Add taxonomic prefixes and export ----
    # Add taxonomic prefixes
    taxLevel_short <- c("k__","p__","c__","o__","f__","g__","s__")
    taxLevels <- c("kingdom","phylum","class","order","family","genus","species")
    res_df <- tibble(tax=character(0),
                     count=double(0))
    #create df with only the tax levels and the added prefixes
    tmp <- combined_reads_df %>% 
      select(-tmpTax) %>% 
      distinct() %>% 
      as.data.frame()
    for(j in seq(2,ncol(tmp))){
      tmp[,j] <- paste0(taxLevel_short[j-1],tmp[,j])
    }
    taxReference_df <- 
      tmp %>% 
      mutate(across(.cols = -count,.fns = ~map(., str_replace_all,pattern = ".__NA",replacement = NA_character_))) %>% 
      unnest(cols=c(kingdom, phylum, class, order, family, genus, species)) %>% 
      as_tibble()
    
    #Convert to the desired output style matrix
    if(output_style == "metaphlan"){
      # Metaphlanize the mapseq output (give total sum of counts at every taxonomic combination)
      res_df <- .f_metaphlanize_mapseq(taxReference_df = taxReference_df)
      
    }else if(output_style == "standard"){
      res_df <- 
        taxReference_df %>% 
        unite("tax",kingdom:species,sep = "|",na.rm = T,remove = F) %>% 
        select(tax,count)
    }
    
    
    ### merge with count_df
    count_df <- 
      count_df %>% 
      full_join(.,res_df %>% rename(!!as.symbol(c.sample) := count),
                by="tax")
    
    #remove c.combined (to prevent for whatever reason that a failed sample leads to duplication of the previous one)
    rm(res_df,c.combined,taxReference_df,tmp)
    pb$tick()
  }
  ### Create output matrix
  count_mat <- count_df %>% 
    as_tibble() %>% 
    column_to_rownames("tax") %>% 
    replace(is.na(.),0) %>% 
    as.matrix()
  
  return(count_mat)
}


.f_read_in_files_kraken2_TableOutput <- function(path_to_folder,tax.level){
  ### Read in kraken2 result files and return matrix with counts per bacteria and sample
  
  #convert tax.level to kraken2 compatible taxonoic symbol in order to subset the DF
  tax.sym <- toupper(str_extract(tax.level,pattern = "^[A-z]"))
  
  var.names <- c("pct.total","counts.sum","counts.only.here","tax.symbol","tax.ID","tax.name")
  file_list <- list.files(path=path_to_folder)  
  ### initialize output df
  counts.df <- tibble(tax.name = character(0))
  ### iterate over every file and select counts at the selected tax level
  pb <- progress_bar$new(total=length(file_list))
  total.counts.mapped <- tibble(Sample_ID = !!gsub(x = file_list,pattern = ".txt",replacement = ""),tot.counts.mapped = double(length(file_list)))
  i <- 1
  for(i in seq(1,length(file_list))){
    #message(i)
    ### read in file
    #chec if file is empty
    if(file.info(paste0(path_to_folder,file_list[i]))$size==0){
      message(paste0("\nSample ",file_list[i]," is empty - skipping"))
      pb$tick()
      next
    }
    #c.f <- read_tsv(paste0(path_to_folder,file_list[i]),col_names = var.names,col_types = cols())
    c.f <-read.csv(paste0(path_to_folder,file_list[i]),header = FALSE,sep = "\t",comment.char = "",check.names = FALSE,col.names = var.names) %>% 
      mutate(tax.name = trimws(.$tax.name),
             rowNames = seq(1,nrow(.)))
    
    
    c.total.counts.mapped <- c.f %>% filter(tax.name == "root") %>% pull(counts.sum)
    if(!(is_empty(c.total.counts.mapped))){
      total.counts.mapped[i,2] <- c.total.counts.mapped
    }  
    ### extract all bacterial counts
    domains_vec <- 
      c.f %>% 
      filter(tax.symbol == "D") %>% 
      arrange(rowNames) %>% 
      select(tax.name,rowNames) %>% 
      deframe
    rowNamesVec <- as.numeric(domains_vec)
    
    #get start and end positions of individual domains (bacteria and archaea; viruses are not relevant)
    bac.start <- as.numeric(domains_vec["Bacteria"])
    bac.end <- rowNamesVec[which(rowNamesVec==bac.start)+1]-1
    archaea.start <- as.numeric(domains_vec["Archaea"])
    archaea.end <- rowNamesVec[which(rowNamesVec==archaea.start)+1]-1
    
    if(is.na(bac.start)){
      message(paste0("\nNo bacteria profiled in sample ",file_list[i], " - skipping"))
      pb$tick()
      next
    }
    
    #check if the respective "end" number is empty -> if so, it must be changed to the last row of the file
    if(!(is.na(bac.start))){
      if(is.na((bac.end))){
        bac.end <- nrow(c.f)
      } 
    }
    if(!(is.na(archaea.start))){
      if(is.na(archaea.end)){
        archaea.end <- nrow(c.f)
      }
    }
    #subset data too keep only bacterial and archaeal reads
    if(!(is.na(archaea.start))){#get bacterial and archaeal reads
      bac.reads <- c.f[c(bac.start:bac.end,
                         archaea.start:archaea.end),]
    }else{#keep only bacterial reads
      bac.reads <- c.f[c(bac.start:bac.end),]
    }
    
    ### select counts and tax names
    c.counts.df <- 
      bind_rows(
        bac.reads %>% filter(tax.symbol == tax.sym) %>% select(tax.name,counts.sum), #select counts at genus level
        bac.reads %>% filter(tax.symbol == "D") %>% select(tax.name,counts.sum) %>% #select "Bacterial" and "Archeal" counts, sum them togehter as "Bacteria" and merge with counts data
          mutate(tax.name = "Bacteria") %>% 
          group_by(tax.name) %>% 
          summarise(counts.sum = sum(counts.sum))) %>% 
      arrange(-counts.sum) %>% 
    rename(!!gsub(x = file_list[i],pattern = ".txt",replacement = "") := counts.sum)
    
    # Add to data from other samples
    counts.df <- suppressMessages(full_join(counts.df,c.counts.df,by="tax.name"))
    pb$tick()
    #print(i)
  }
  ### convert to output matrix
  counts.df <- counts.df %>% mutate_at(vars(-tax.name),as.numeric)
  counts.mat <- as.matrix(counts.df[,-1])
  rownames(counts.mat) <- (counts.df$tax.name)
  counts.mat[is.na(counts.mat)] <- 0
  #remove everything after the first point in the samplenames
  colnames(counts.mat) <- sub(colnames(counts.mat),pattern = ".kraken2_report",replacement = "")
  
  res.list <- list(counts.mat,total.counts.mapped)
  names(res.list) <- c("counts.mat","total.counts.mapped")
  #return(counts.mat)
  return(res.list)
}

.f_read_in_files_kraken2 <- function(path_to_folder,tax.level=NULL){
  ### Read in kraken2 result files generated using the --use-mpa-style option (full taxonomy)
  
  file_list <- list.files(path=path_to_folder)  
  
  #tax.level not necessary anymore but still kept as potential input option in order to not break the other scripts
  tax.level = NULL
  
  ### initialize output df
  count_df <- tibble(tax = character(0))
  ### iterate over every file and select counts at the selected tax level
  pb <- progress_bar$new(total=length(file_list))
  i <- 1
  for(i in seq(1,length(file_list))){
    #message(i)
    ### read in file
    #chec if file is empty
    if(file.info(paste0(path_to_folder,file_list[i]))$size==0){
      message(paste0("\nSample ",file_list[i]," is empty - skipping"))
      pb$tick()
      next
    }
    c.f <-read.csv(paste0(path_to_folder,file_list[i]),header = FALSE,sep = "\t",comment.char = "",check.names = FALSE) %>% 
      rename(tax = V1,count = V2) %>% 
      as_tibble() %>% 
      filter(str_detect(tax,"d__Bacteria|d__Archaea")) %>% ###keep only Archaea and Fungi
      rename(!!as.symbol(str_remove(file_list[i],".kraken2.*")) := count)
    
    # join with count_df
    count_df <- 
      full_join(count_df,
                c.f,
                by="tax")
    pb$tick()
    #print(i)
  }
  
  ### convert to output matrix
  count_mat <- 
    count_df %>% 
    replace(is.na(.),0) %>% 
    column_to_rownames("tax") %>% 
    as.matrix()
  
  return(count_mat)
}

.f_read_in_files_PathSeq <- function(path_to_folder,tax.level){
  ### Read PathSeq output files, select tax level (eg genus) and return 2 matrices:
  #1) matrix with pathseq scores values
  #2) with count values of !unambiguously! mapped reads
  
  
  file_list <- list.files(path=path_to_folder,pattern = "\\.txt$")  
  if(is_empty(file_list)){
    file_list <- list.files(path=path_to_folder,pattern = "\\.scores$")  
  }
  
  
  
  ### initialize output df
  score.df <- tibble(tax.name = character(0),tax_id=double(0))
  counts.df <- tibble(tax.name = character(0),tax_id=double(0))
  
  ### iterate over every file and select counts at the selected tax level
  empty.counter <- 0
  pb <- progress_bar$new(total=length(file_list))
  i <- 1
  
  for(i in seq(1,length(file_list))){
    ## read in file
    #c.f <-read_tsv(paste0(path_to_folder,file_list[i]),col_types = cols())
    c.f <-read.csv(paste0(path_to_folder,file_list[i]),header = TRUE,sep = "\t",comment.char = "",check.names = FALSE) %>% 
      mutate(name = str_replace(name,pattern = "\\_",replacement = " "))
    
    #check if file is empty
    if(nrow(c.f) == 0){
      message(paste0("\nSample ",file_list[i]," is empty - skipping"))
      pb$tick()
      empty.counter <- empty.counter+1
      next
    }
    
    ### select counts and tax names
    #Edit 22-09-01: When Archaeal reads also should be considered: 
    #Select "Archaea" and "Bacteria", sum up "score" and "unambigous" of the "superkingdom" and then normalize every tax.lvl counts against this number
    
    c.combined.df <- 
      c.f %>% 
      filter(type == "superkingdom",
             kingdom %in% c("Bacteria","Archaea")) %>% 
      mutate(name = "Bacteria") %>% 
      group_by(name) %>% 
      summarise(score = sum(score),
                unambiguous = sum(unambiguous)) %>% 
      add_column(tax_id = 0000) %>%  #just to have a tax_id
      bind_rows(.,c.f %>% filter(kingdom %in% c("Bacteria","Archaea"),
                                 type == tax.level)) %>% 
      mutate(score_normalized = score/score[1]*100)
    
    # Extract the pathseq scores 
    c.score.df <-
      c.combined.df %>% 
      select(name,tax_id,score) %>%
      rename(!!gsub(x = file_list[i],pattern = ".pathseq.txt",replacement = "") := score,
             tax.name = name)
    
    # also unambiguous counts
    c.counts.df <-
      c.combined.df %>% 
      select(name,tax_id,unambiguous) %>%
      rename(!!gsub(x = file_list[i],pattern = ".pathseq.txt",replacement = "") := unambiguous,
             tax.name = name)
    
    
    #Join with data from other samples
    score.df <- full_join(score.df,c.score.df,by=c("tax.name","tax_id"))
    counts.df <- full_join(counts.df,c.counts.df,by=c("tax.name","tax_id"))
    
    pb$tick()
  }
  
  #'In the PathSeq DB of Dohlman et al there seems to be a naming issue:
  #'Both,tax_ID 1937007 and tax_ID 1980680 correspond to the genus "Ileibacterium". 
  #'According to NCBI taxonomy, tax_ID 1980680 refers to "Alterileibacterium
  #'--> correct the naming manually
  score.df <- score.df %>% mutate(tax.name = case_when(tax_id == 1980680 ~ "Alterileibacterium",
                                                       TRUE ~ tax.name)) %>% 
    select(-tax_id)
  counts.df <- counts.df %>% mutate(tax.name = case_when(tax_id == 1980680 ~ "Alterileibacterium",
                                                         TRUE ~ tax.name)) %>% 
    select(-tax_id)
  
  
  ### convert to output matrices and return list
  score.df <- score.df %>% mutate_at(vars(-tax.name),as.numeric)
  score.mat <- as.matrix(score.df[,-1])
  rownames(score.mat) <- (score.df$tax.name)
  score.mat[is.na(score.mat)] <- 0
  #remove all-zero rows
  score.mat <- score.mat[rowSums(score.mat)>0,,drop=F]
  #remove everything after the first point in the samplenames
  colnames(score.mat) <- sub(colnames(score.mat),pattern = ".pathseq.scores",replacement = "")
  
  counts.df <- counts.df %>% mutate_at(vars(-tax.name),as.numeric)
  counts.mat <- as.matrix(counts.df[,-1])
  rownames(counts.mat) <- (counts.df$tax.name)
  counts.mat[is.na(counts.mat)] <- 0
  counts.mat <- counts.mat[rowSums(counts.mat)>0,,drop=F]
  #remove everything after the first point in the samplenames
  colnames(counts.mat) <- sub(colnames(counts.mat),pattern = ".pathseq.scores",replacement = "")
  
  
  res.list <- list(score.mat,counts.mat)
  names(res.list) <- c("score.mat","counts_unambiguous.mat")
  
  return(res.list)
}

.f_read_in_files_mOTUs <- function(path_to_folder){
  #manually merges mOTUs raw counts from varius runs and returns count matrix with full taxonomy
  #return matrix with full length taxonomy (kingdom --> mOTUs)
  
  file_list <- list.files(path=path_to_folder,pattern = ".txt")
  message("Collating mOTUs results for ",length(file_list), " samples")
  
  count_df <- tibble(tax=character(0))
  pb <- progress_bar$new(total=length(file_list))
  i <- 1
  for(i in seq(1,length(file_list))){
    if(file.info(paste0(path_to_folder,file_list[i]))$size==0){
      message(paste0("Sample ",file_list[i]," is empty - skipping"))
      next
    }
    tmp.file <- read_tsv(paste0(path_to_folder,file_list[i]),skip = 2,show_col_types = F)
    colnames(tmp.file)[1] <- "mOTU"
    colnames(tmp.file)[2] <- "taxonomy"
    colnames(tmp.file)[3] <- str_remove(file_list[i],".txt")
    
    #merge "mOTU" column with the taxonomy column and merge by rownames
    tmp.file <- suppressMessages(tmp.file %>% 
                                   mutate_at(vars(-taxonomy,-mOTU),as.numeric) %>% 
                                   unite("tax",taxonomy,mOTU,sep = "|"))
    
    count_df <- full_join(count_df,tmp.file,by="tax")  
    
    pb$tick()
  }
  #rename the "unassigned" fraction properly
  count_df <- count_df %>% 
    mutate(tax = case_when(tax == "unassigned|unassigned"~"unassigned",TRUE~tax))
  #### convert to output matrix ####
  count_mat <- count_df %>% 
    column_to_rownames("tax") %>% 
    as.matrix()
  #remove all-zero rows
  count_mat <- count_mat[rowSums(count_mat)>0,,drop=FALSE]
  #clean sample names
  colnames(count_mat) <- str_remove(colnames(count_mat),pattern = ".motus")
  
  return(count_mat)
}
.f_read_in_files_mTAGsProfile_merged <- function(path_to_folder = mTAGs_profile_folder,tax.level = "genus"){
  # reads in merged mTAGs profile output file
  # So far: works only on genus level! 
  
  
  
  file_list <- list.files(path_to_folder)
  file_list <- file_list[str_detect(file_list,tax.level)]
  
  if(length(file_list)>1){
    stop("More than one input files with the given tax. level")
  }
  if(!(tax.level == "genus")){
    stop("Select <genus> as tax.level")
  }
  
  count.df <- tibble(genus = character())
  i <- 1
  
  tmp.file <- read.table(paste0(path_to_folder,file_list[i]),sep = "\t",comment.char = "",header = T,check.names = FALSE)
  #colnames(tmp.file) <- c("taxonomy","count")
  colnames(tmp.file)
  tmp.counts <- tmp.file %>% select(!!colnames(tmp.file)[1]) %>% 
    separate(!!colnames(tmp.file)[1], into=c('root', 'domain', 'phylum', 
                                             'class', 'order', 'famiy', 'genus'), 
             sep=';', fill='right')
  tmp.counts <- as_tibble(apply(tmp.counts,MARGIN = 2, function(x){
    sub(".*__",x = x,replacement = "")
  }))
  
  ###Add the counts of the samples
  tmp.counts <- cbind(tmp.counts,tmp.file[,-1])
  
  
  tmp.genus <- tmp.counts %>% filter(root == "Root",domain == "Bacteria") %>% 
    select(genus,!!colnames(tmp.file)[2:ncol(tmp.file)]) %>% 
    group_by(genus) %>% 
    summarise_at(vars(!!(colnames(.)[2:ncol(.)])),
                 funs(sum)) 
  
  #compute total bact counts (sum of unknown + genus resolved)
  tot.bact.counts <- as.data.frame(t(as.matrix(colSums(tmp.genus[,-1],na.rm = TRUE))))
  tot.bact.counts$genus <- "Bacteria"
  tot.bact.counts <- tot.bact.counts %>% relocate(genus)
  count.df <- bind_rows(tot.bact.counts,tmp.genus) %>% 
    filter(!(genus %in% c("unknown","uncultured")))
  
  #create output matrix (assign rownames and remove column with tax names)
  rownames(count.df) <- count.df$genus
  count.df$genus <- NULL
  count.mat <- as.matrix(count.df)
  count.mat[is.na(count.mat)] <- 0
  #remove everything after the first point in the samplenames
  colnames(count.mat) <- sub(colnames(count.mat),pattern = ".bins",replacement = "")
  return(count.mat)
}

.f_read_in_libsize <- function(path_to_folder){
  file_list <- list.files(path_to_folder)
  
  df.libsize <- tibble()
  for(i in seq(1,length(file_list))){
    tmp.file <- read.table(paste0(path_to_folder,file_list[i]),sep = "\t",comment.char = "",header = F)
    tmp.file$Sample_ID <- sub(file_list[i],pattern = ".libsize.txt",replacement = "")
    colnames(tmp.file)[1] <- "libsize_counts"
    df.libsize <- bind_rows(df.libsize,tmp.file)
  }
  df.libsize <- df.libsize %>% relocate(Sample_ID)
  return(df.libsize)
}

.f_read_in_lib_layout <- function(path_to_folder){
  file_list <- list.files(path_to_folder)
  df.lib_layout <- tibble()
  for(i in seq(1,length(file_list))){
    #message(i)
    if(file.info(paste0(path_to_folder,file_list[i]))$size==0){
      message(paste0("Sample ",file_list[i]," is empty - skipping"))
      next
    }else{
      tmp.file <- read.table(paste0(path_to_folder,file_list[i]),sep = "\t",comment.char = "",header = F)
      tmp.file$Sample_ID <- sub(file_list[i],pattern = ".is_paired.txt",replacement = "")
      colnames(tmp.file)[1] <- "lib_layout"
      df.lib_layout <- bind_rows(df.lib_layout,tmp.file)
    }
  }
  df.lib_layout <- df.lib_layout %>% relocate(Sample_ID)
  return(df.lib_layout)
}

.f_read_in_flagstats <- function(path_to_folder){
  #No needed anymore in vknight build >= 827290457b
  #read in the total number of reads from the flagstats output -> Corresponds to ALL reads in the FastQ input file (before QC)
  file_list <- list.files(path_to_folder,pattern = "*.flagstats.txt")
  df.flagstats <- tibble()
  for(i in seq(1,length(file_list))){
    tmp.file <- read.table(paste0(path_to_folder,file_list[i]),sep = "\t",comment.char = "",header = F)
    #Extract total number of reads (1st number in 1st row of the flagstats output)
    tmp.total.reads <- as.numeric(gsub("\\ +.*", "",tmp.file[1,1]))
    names(tmp.total.reads) <- sub(file_list[i],pattern = ".flagstats.txt",replacement = "")
    
    tmp.res <- enframe(tmp.total.reads,name = "ID",value = "total.libsize")
    
    df.flagstats <- bind_rows(df.flagstats,tmp.res)
  }
  return(df.flagstats)
}

.f_read_in_raw_counts_number <- function(path_to_folder){
  #'take raw_counts computed by vknight (>build 827290457b) when preprocessing the files with multiQC
  #'exports a tibble with Sample_ID and raw_read_count
  file_list <- list.files(path_to_folder)
  
  df.raw_counts <- tibble(Sample_ID=character(),libsize_raw_counts=double())
  for(i in seq(1,length(file_list))){
    tmp.file <- read.table(paste0(path_to_folder,file_list[i]),sep = "\t") %>% 
      separate(V1,sep = " ",into = c("paired_info","V1")) %>% 
      mutate(paired_info = as.numeric(paired_info))
    tmp.raw_counts <- tibble(Sample_ID = str_remove(file_list[i],pattern = ".txt"),
                             libsize_raw_counts = tmp.file$paired_info * tmp.file$V2)
    
    df.raw_counts <- bind_rows(df.raw_counts,tmp.raw_counts)
  }
  return(df.raw_counts)
}

.f_read_in_read_counter <- function(path_to_folder,min_genes,marker_genes.df,tax.level=NULL,norm_to_gene_length=TRUE){
  #Creates count matrix from folder of read_counter output files considering all species that were hit with at least <min_genes> marker genes
  #normalization: Assumes read_counter was run with "-y insert.raw_counts"; if norm_to_gene_length == TRUE, normalization is performed via:
  #1) sum up counts per species over all marker genes
  #2) divide by total length of marker genes for the given species
  #3) multiply by average length of marker genes over all species in the DB
  
  #deprecated but left as optional in order to not break the scripts
  tax.level <- NULL
  
  stopifnot(is_tibble(marker_genes.df))
  
  file_list <- list.files(path=path_to_folder)
  message(paste0(length(file_list)," files in the given folder"))
  ### initialize output df
  counts_df <- tibble(tax = character(0))
  ### iterate over every file and select counts at the selected tax level
  pb <- progress_bar$new(total=length(file_list))
  for(i in seq(1,length(file_list))){
    #message(i)
    sample.name <- file_list[i]
    tmp.file <- tryCatch(
      {data.table::fread(file = paste0(path_to_folder,sample.name),sep = "\t",skip = 1)},
      error=function(e){
        c.fwd <- data.frame(matrix(ncol = 14,nrow = 0))
      }
    )
    if(!(nrow(tmp.file)>0 & ncol(tmp.file)>0)){
      message(paste0("Sample ",sample.name," is empty - skipping"))
      pb$tick()
      next
    }
    
    file.name <- str_remove_all(sample.name,pattern = "\\.read_counter.txt|\\.txt")
    c.counts <- suppressMessages(.f_read_counter_to_count_mat(read_counter_out.file = tmp.file,
                                                              min_genes = min_genes,
                                                              marker_genes.df = marker_genes.df,
                                                              tax.level = tax.level,
                                                              norm_to_gene_length=norm_to_gene_length))
    counts_df <- suppressMessages(full_join(counts_df,c.counts %>% rename(!!as.symbol(file.name) := avg)))
    pb$tick()
  }
  
  #create output matrix
  counts_mat <- 
    counts_df %>% 
    as.data.frame() %>% 
    column_to_rownames("tax") %>% 
    replace(is.na(.),0) %>% 
    as.matrix()
    
  return(counts_mat)
}

.f_read_counter_to_count_mat <- function(read_counter_out.file,file.name=NULL,min_genes,marker_genes.df,tax.level="genus",norm_to_gene_length=TRUE){
  
  #deprecated
  file.name <- NULL
  
  #takes output file from read_counter tool and returns the median counts aggregated by the selected taxonomic rank
  c.df <-
    read_counter_out.file %>%
    #mutate(V1 = gsub(pattern = "\\-.*",replacement = "",x = V1)) %>%
    separate(V1,into=c('kingdom', 'phylum', 'class',
                       'order', 'family', 'genus', 'species'),
             sep=";",fill = "right") %>%
    separate(species,into=c('species','gene'),
             sep="[-](?=[^-]+$)",fill = "right") %>%#splits only at the last occurence of the character "-"
    rename(read_count  = V2)
  
  #get bacterial counts (only valid if read_counter is run with mode "-y insert.raw_counts")
  totalBactCounts <- sum(c.df$read_count)
  
  #Determine how many marker genes were identified at species level to threshold species
  c.N_genes <- (c.df %>% select(species,gene) %>%  group_by(species) %>% summarise(n=n()))
  
  #select species and compute average read_count by species (over all marker genes that are in the database for the given species)
  c.by_species <-
    c.df %>%
    select(species,gene,read_count) %>%
    filter(species %in% (c.N_genes %>% filter(n>=min_genes) %>% pull(species)))
  
  # #return NULL if minimal marker genes is not met by any taxon
  # if(nrow(c.by_species)<1){
  #   message("At selected tax.level no taxa meets the min.marker.gene criteria")
  #   return(NULL)
  # } 
  
  if(isTRUE(norm_to_gene_length)){
    #average read_counts over all marker genes and also consider number of marker genes with zero counts
    avg.gene_length <- mean(marker_genes.df$sum_length_bySpecies)
    c.by_species <- 
      c.by_species %>% 
      group_by(species) %>% 
      summarise(sum = sum(read_count)) %>% 
      left_join(.,marker_genes.df) %>%
      mutate(avg = (sum/sum_length_bySpecies)*avg.gene_length)
    
  }else{
    #take median
    c.by_species <- c.by_species %>% 
      group_by(species) %>% 
      summarise(avg = median(read_count))
  }
  
  #select all species that passed the filter and return the entire taxonomic annotation
  c.tax.counts <-
    c.df %>%
    select(-read_count,-gene) %>% 
    inner_join(.,c.by_species %>% select(species,avg)) %>%
    distinct() %>%
    unite("tax",kingdom:species,sep = "|")
  
  #add "unassigned" bacterial counts (assigned to species not passing the marker gene threshold)
  unassigned_df <- tibble(tax = "unassigned",avg = totalBactCounts - sum(c.tax.counts$avg))
  
  res_df <- 
    bind_rows(c.tax.counts,unassigned_df)
  
  return(res_df)
}

# .f_ASV_to_tax <- function(lvl,asv_table,mseq.output){
#   #Function that converts an asv-count matrix (e.g. from a dada2 result) to a count matrix at selected tax-level (lvl), 
#   #by taking the taxonomic infroamtion from the mseq.output file (ASVs must be in column 1 and tax-tree must be in column 14 (-> default when running mapseq with --simple)))
#   
#   stopifnot(lvl %in% c('kingdom', 'phylum', 'class', 
#                        'order', 'family', 'genus', 'species'))
#   tax <- .f_resolve_ASVs(mseq.output = mseq.output)
#   tax
#   tax <- tax %>% 
#     filter(ASV %in% rownames(asv_table))
#   if (any(colSums(asv_table) == 0)){
#     asv_table <- asv_table[,colSums(asv_table) > 0]
#   }
#   
#   groups <- tax %>% 
#     filter(!is.na(!!as.symbol(lvl))) %>% 
#     pull(!!as.symbol(lvl)) %>% 
#     unique
#   groups
#   mat.new <- matrix(0, nrow=length(groups), ncol=ncol(asv_table),
#                     dimnames = list(groups, colnames(asv_table)))
#   # Selects all ASVs that correspond to the current group (e.g current genera) and sums the counts
#   for (g in groups){
#     mat.new[g,] <- colSums(asv_table[
#       tax %>% 
#         filter(!!as.symbol(lvl)==g) %>% 
#         pull(ASV), ,drop=FALSE])
#   }
#   #Add counts of samples that are not resolved at given tax_level
#   #This is represented by the difference of colSums
#   not_resolved <- as.matrix(t(colSums(asv_table)-colSums(mat.new)))
#   rownames(not_resolved) <- "not_resolved"
#   mat.new <- rbind(not_resolved,mat.new)
#   
#   #remove the underscore with the tax_level identifier (e.g. "g__")
#   rownames(mat.new) <- str_remove(rownames(mat.new),pattern = "[a-zA-Z]__")
#   
#   print(summary(colSums(mat.new)))
#   return(mat.new)
# }

.f_map_asv_to_tax <- function(asv_tsv,path_to_ASV.mseq,tax_lvl){
  
  
  #Read in the raw mapseq output and create tsv file
  c.fwd <- tryCatch(
    {data.table::fread(file = path_to_ASV.mseq,skip = 1,header = T)},
    error=function(e){
      c.fwd <- data.frame(matrix(ncol = 15,nrow = 0))
    }
  )
  colnames(c.fwd) <- as.character(seq(1,ncol(c.fwd)))
  c.combined <- c.fwd[,c(1,14)]
  colnames(c.combined) <- c("ASV","tax_info")
  
  tax.split <- data.table::tstrsplit(x = c.combined$tax_info,";")
  #make list with 7 entries to match taxonomic tree
  if(length(tax.split)<7){tax.split[(length(tax.split)+1):7] <- NA}#fill with NAs up to species level
  c.combined[,c('kingdom', 'phylum', 'class','order', 'family', 'genus', 'species')] <- tax.split
  
  #now take the asv_table, merge with the taxonomy and calculate relative abundance at taxonomic level
  tax_counts_df <- 
    asv_tsv %>% 
    rownames_to_column("ASV") %>% 
    as_tibble() %>% 
    gather(-ASV,key = "Sample_ID",value = "count") %>% 
    inner_join(.,c.combined %>% select(-tax_info)) %>% 
    select(Sample_ID,count,!!as.symbol(tax_lvl)) %>% 
    group_by(Sample_ID,!!as.symbol(tax_lvl)) %>% 
    summarise(count = sum(count)) %>% 
    ungroup() %>% 
    mutate(!!as.symbol(tax_lvl) := str_remove(!!as.symbol(tax_lvl),pattern = "[a-zA-Z]__")) %>% 
    mutate(!!as.symbol(tax_lvl) := case_when(is.na(!!as.symbol(tax_lvl))~"not_resolved",
                                             TRUE~!!as.symbol(tax_lvl))) %>% 
    arrange(-count) %>% 
    group_by(Sample_ID) %>% 
    mutate(count_rel = count/sum(count)) %>% 
    ungroup()
  
  #create count matrix  
  counts_mat <- 
    tax_counts_df %>% 
    select(Sample_ID,!!as.symbol(tax_lvl),count) %>% 
    pivot_wider(names_from = Sample_ID,values_from = count,values_fill = 0) %>% 
    column_to_rownames(tax_lvl) %>% 
    as.matrix()
  
  #add bacterial counts
  bacMat <- (as.matrix(t(colSums(counts_mat))))
  rownames(bacMat) <- "Bacteria"
  counts_mat <- rbind(bacMat,counts_mat)
  
  #create relCount matrix
  countsRel_mat <- 
    tax_counts_df %>% 
    select(Sample_ID,!!as.symbol(tax_lvl),count_rel) %>% 
    pivot_wider(names_from = Sample_ID,values_from = count_rel,values_fill = 0) %>% 
    column_to_rownames(tax_lvl) %>% 
    as.matrix()
  
  #so far: only export count matrix (do the read in later with vknight read in functions)
  return(counts_mat)
  
}








