library(stringr)
library(progress)

.f_read_in_files_mapseq <- function(path_to_folder){
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

.f_read_in_mapseq_raw <- function(path_to_folder,tax.lvl="genus",fname_mode="vknight",account_for_paired_reads=TRUE){
  #fname_mode defines how the filenames were constructed: 
  #vknight: SAMPLE_ID_R[0-9]_bac_ssu.mseq
  #manual: SAMPLE_ID_1.mseq
  #takes mapseq output files as input and returns a matrix with total counts at the selected taxonomic_level
  #if reads are paired, the counts are averaged between the forward and the reverse reads (fwd = 4x fuso, rev = 10x fuso --> 7x fuso) [if account_for_paired_reads == TRUE]
  #!!! Works only if .mseq file was computed with --outfmt simple !!!#
  
  flist <- list.files(path_to_folder,pattern = ".mseq")
  if(length(flist)<1){message("No .mseq files in given directory")
    stop
  }
  
  #Get the filenames that are expected in the .mseq folder
  if(fname_mode == "vknight"){
    sample.names <- unique(str_remove(flist,pattern = "_R[0-9]_bac_ssu.mseq"))
    fwd_file_ending <- "_R1_bac_ssu.mseq"
    rev_file_ending <- "_R2_bac_ssu.mseq"
  }else if(fname_mode == "manual"){
    sample.names <- unique(str_remove(flist,pattern = "_[0-9].mseq"))
    fwd_file_ending <- "_1.mseq"
    rev_file_ending <- "_2.mseq"
    # }else if(fname_mode == "dada2"){
    #   sample.names <- unique(str_remove(flist,pattern = ".mseq"))
    #   fwd_file_ending <- ".mseq"
    #   #rev_file_ending <- "_2.mseq"
    # }
  }else{
    message("Incorrect fname_mode")
    stop()
  }
  message(paste0("Importing ",length(sample.names)," samples"))
  res.df <- tibble(!!tax.lvl := character())
  pb <- progress_bar$new(total=length(sample.names))
  i <- 306
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
    
    c.combined <- rbind(c.fwd[,c(1,14)],
                        c.rev[,c(1,14)])
    colnames(c.combined) <- c("read_id","tax_info")
    
    tax.split <- data.table::tstrsplit(x = c.combined$tax_info,";")
    #make list with 7 entries to match taxonomic tree
    if(length(tax.split)<7){tax.split[(length(tax.split)+1):7] <- NA}#fill with NAs up to species level
    c.combined[,c('kingdom', 'phylum', 'class','order', 'family', 'genus', 'species')] <- tax.split
  
    if(all(is.na(c.combined[[tax.lvl]]))){
      next
      pb$tick()
    }#skip if no counts at tax level are present
    #compute absolute abudnance
    tax.counts <- c.combined %>%
      select(-tax_info) %>% 
      filter(str_detect(kingdom,pattern = "Bacteria")) %>%
      group_by(!!as.symbol(tax.lvl)) %>%
      mutate(!!as.symbol(tax.lvl) := str_remove(!!as.symbol(tax.lvl),pattern = "[a-zA-Z]__")) %>% ##remove "g__" from tax.lvl column name
      mutate(!!as.symbol(tax.lvl) := case_when(!!as.symbol(tax.lvl) == "" ~ NA_character_,
                                               TRUE~!!as.symbol(tax.lvl))) %>% 
      summarize(counts = length(read_id)) %>%
      arrange(desc(counts)) %>% 
      mutate(!!tax.lvl := case_when(!(is.na(!!as.symbol(tax.lvl))) ~ !!as.symbol(tax.lvl),
                                    TRUE ~ "not_resolved"))
    rm(c.combined,c.fwd,c.rev)
    #add bacterial counts
    bact.counts <- tibble(genus = "Bacteria",counts = sum(tax.counts$counts))
    tax.counts <- bind_rows(bact.counts,tax.counts)
    
    #If reads are paired, divide counts by 2 (assuming )
    if(isTRUE(account_for_paired_reads) & sample.type=="paired"){
      tax.counts <- tax.counts %>% mutate(counts = round(counts/2,0))
    }
    colnames(tax.counts)[2] <- c.sample
    
    res.df <- suppressMessages(full_join(res.df,tax.counts,by=tax.lvl))
    pb$tick()
    
  }
  
  #create output matrix
  res.mat <- 
    res.df %>% 
    as.data.frame() %>% 
    column_to_rownames(tax.lvl) %>% 
    replace(is.na(.),0) %>% 
    as.matrix()
  
  return(res.mat)
}


.f_read_in_files_kraken2 <- function(path_to_folder,tax.level){
  ### Read in kraken2 result files and return matrix with counts per bacteria and sample
  
  if(tax.level == "genus"){
    tax.sym <- "G"
  }else if(tax.level == "species"){
    tax.sym <- "S"
  }else{
    stop("Please enter a tax level (genus or species)")
  }
  
  var.names <- c("pct.total","counts.sum","counts.only.here","tax.symbol","tax.ID","tax.name")
  file_list <- list.files(path=path_to_folder)  
  ### initialize output df
  counts.df <- tibble(tax.name = character(0))
  ### iterate over every file and select counts at the selected tax level
  pb <- progress_bar$new(total=length(file_list))
  total.counts.mapped <- tibble(Sample_ID = !!gsub(x = file_list,pattern = ".txt",replacement = ""),tot.counts.mapped = double(length(file_list)))
  i <- 81
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
    c.f <-read.csv(paste0(path_to_folder,file_list[i]),header = FALSE,sep = "\t",comment.char = "",check.names = FALSE,col.names = var.names)
    c.f$tax.name <- str_replace_all(c.f$tax.name,pattern = " ",replacement = "")
    c.total.counts.mapped <- c.f %>% filter(tax.name == "root") %>% pull(counts.sum)
    if(!(is_empty(c.total.counts.mapped))){
      total.counts.mapped[i,2] <- c.total.counts.mapped
    }  
    ### extract all bacterial counts
    bac.start <- which(c.f[,6] == "Bacteria")
    if(is_empty(bac.start)){
      message(paste0("\nNo bacteria profiled in sample ",file_list[i], " - skipping"))
      pb$tick()
      next
    }
    
    virus.start <- which(c.f[,6] == "Viruses") #get row in table at which "viruses" start
    archea.start <- bac.end <- which(c.f[,6] == "Archaea") #get row in table at which "archea" start
    
    bac.end <- min(virus.start,archea.start)-1 #get last row with entries for bacteria
    if(!(is.finite(bac.end))){
      bac.end <- nrow(c.f)
    }
    #subset to keep only the entries mapped to bacteria
    bac.reads <- c.f[(bac.start:bac.end),]  
    
    ### select counts and tax names
    c.counts.df <- bac.reads %>% filter(tax.symbol == tax.sym | tax.symbol == "D") %>% select(tax.name,counts.sum) %>% 
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


.f_read_in_files_PathSeq <- function(path_to_folder,tax.level){
  ### Read PathSeq output files, select tax level (eg genus) and return 2 matrices:
  #1) matrix with score_normalized values
  #2) with count values of !unambiguously! mapped reads
  
  
  file_list <- list.files(path=path_to_folder,pattern = "\\.txt$")  
  if(is_empty(file_list)){
    file_list <- list.files(path=path_to_folder,pattern = "\\.scores$")  
  }
  if(tax.level == "genus"){
    tax.sym <- "genus"
  }else if(tax.level == "species"){
    tax.sym <- "species"
  }else{
    stop("Please enter a tax level (genus or species)")
  }
  
  ### initialize output df
  score.df <- tibble(tax.name = character(0),tax_id=double(0))
  counts.df <- tibble(tax.name = character(0),tax_id=double(0))

  ### iterate over every file and select counts at the selected tax level
  empty.counter <- 0
  pb <- progress_bar$new(total=length(file_list))
  for(i in seq(1,length(file_list))){
    ## read in file
    #c.f <-read_tsv(paste0(path_to_folder,file_list[i]),col_types = cols())
    c.f <-read.csv(paste0(path_to_folder,file_list[i]),header = TRUE,sep = "\t",comment.char = "",check.names = FALSE)
    
    #check if file is empty
    if(nrow(c.f) == 0){
      message(paste0("\nSample ",file_list[i]," is empty - skipping"))
      pb$tick()
      empty.counter <- empty.counter+1
      next
    }
    
    ### select counts and tax names
    c.score.df <-
      c.f %>% filter(type == tax.sym | type == "superkingdom",
                     kingdom == "Bacteria") %>%
      select(name,tax_id,score_normalized) %>%
      rename(!!gsub(x = file_list[i],pattern = ".pathseq.txt",replacement = "") := score_normalized,
             tax.name = name)

    # also unambiguous counts
    c.counts.df <-
      c.f %>% filter(type == tax.sym | type == "superkingdom",
                     kingdom == "Bacteria") %>%
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
  score.mat <- score.mat[rowSums(score.mat)>0,]
  #remove everything after the first point in the samplenames
  colnames(score.mat) <- sub(colnames(score.mat),pattern = ".pathseq.scores",replacement = "")
  
  counts.df <- counts.df %>% mutate_at(vars(-tax.name),as.numeric)
  counts.mat <- as.matrix(counts.df[,-1])
  rownames(counts.mat) <- (counts.df$tax.name)
  counts.mat[is.na(counts.mat)] <- 0
  counts.mat <- counts.mat[rowSums(counts.mat)>0,]
  #remove everything after the first point in the samplenames
  colnames(counts.mat) <- sub(colnames(counts.mat),pattern = ".pathseq.scores",replacement = "")
  
  
  res.list <- list(score.mat,counts.mat)
  names(res.list) <- c("score_normalized.mat","counts_unambiguous.mat")
  
  return(res.list)
}

.f_read_in_files_mOTUs <- function(path_to_folder){
  # reads in output from mOTUs profiler
  #attention: no information about the total sum of bacterial reads is provided by mOTUs
  
  #consider only counts from bacterial small subunit (bac_ssu)
  file_list <- list.files(path=path_to_folder,pattern = ".txt")
  
  #if forward and reverse files are present: merge them and sum up counts
  count.df <- tibble(taxonomy=character(0))
  pb <- progress_bar$new(total=length(file_list))
  i <- 1
  for(i in seq(1,length(file_list))){
    if(file.info(paste0(path_to_folder,file_list[i]))$size==0){
      message(paste0("Sample ",file_list[i]," is empty - skipping"))
      next
    }
    tmp.file <- read.table(paste0(path_to_folder,file_list[i]),sep = "\t",comment.char = "",skip = 2,header = T)
    colnames(tmp.file)[1] <- "taxonomy"
    colnames(tmp.file)[-1] <- str_remove(file_list[i],".txt")
    tmp.file <- tmp.file %>% mutate_at(vars(-taxonomy),as.numeric)
    
    count.df <- full_join(count.df,tmp.file,by="taxonomy")  
    pb$tick()
  }
  #add first row (usually with total bacterial counts, but mOTUs does not provide this information in the output) 
  #-> however, add "total counts on genus level", to have a first row (consistent with the other tools)
  
  #add first row with total bacterial counts: in the case of mOTUs: represents the colSums, since "unassigned" bacteria are also included
  #therefore: add colSums and remove the row with "unassigned"
  
  tmp <- as.data.frame(t(colSums(count.df[,-1],na.rm = TRUE))) %>% 
    add_column(taxonomy = "Bacteria") %>% relocate(taxonomy)
  #remove the row containing "unassigned" bacterial counts --> this is already represented in "Bacteria"
  count.df <- bind_rows(tmp,count.df %>% filter(taxonomy != "unassigned"))
  
  ### convert to output matrix
  count.df <- count.df %>% mutate_at(vars(-taxonomy),as.numeric)
  count.mat <- as.matrix(count.df[,-1])
  rownames(count.mat) <- (count.df$taxonomy)
  count.mat[is.na(count.mat)] <- 0
  
  #remove all-zero rows
  count.mat <- count.mat[rowSums(count.mat)>0,,drop=FALSE]
  #clean sample names
  colnames(count.mat) <- sub(colnames(count.mat),pattern = ".motus",replacement = "")
  return(count.mat)
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
    tmp.file$ID <- sub(file_list[i],pattern = ".libsize.txt",replacement = "")
    colnames(tmp.file)[1] <- "libsize"
    df.libsize <- bind_rows(df.libsize,tmp.file)
  }
  df.libsize <- df.libsize %>% relocate(ID)
  return(df.libsize)
}

.f_read_in_lib_layout <- function(path_to_folder){
  file_list <- list.files(path_to_folder)
  df.lib_layout <- tibble()
  for(i in seq(1,length(file_list))){
    tmp.file <- read.table(paste0(path_to_folder,file_list[i]),sep = "\t",comment.char = "",header = F)
    tmp.file$ID <- sub(file_list[i],pattern = ".is_paired.txt",replacement = "")
    colnames(tmp.file)[1] <- "lib_layout"
    df.lib_layout <- bind_rows(df.lib_layout,tmp.file)
  }
  df.lib_layout <- df.lib_layout %>% relocate(ID)
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
  
  df.raw_counts <- tibble(Sample_ID=character(),raw_counts=double())
  for(i in seq(1,length(file_list))){
    tmp.file <- read.table(paste0(path_to_folder,file_list[i]),sep = "\t") %>% 
      separate(V1,sep = " ",into = c("paired_info","V1")) %>% 
      mutate(paired_info = as.numeric(paired_info))
    tmp.raw_counts <- tibble(Sample_ID = str_remove(file_list[i],pattern = ".txt"),
                             raw_counts = tmp.file$paired_info * tmp.file$V2)
    
    df.raw_counts <- bind_rows(df.raw_counts,tmp.raw_counts)
  }
  return(df.raw_counts)
}

.f_read_in_read_counter <- function(path_to_folder,min_genes,marker_genes.df,tax.level="genus",norm_to_gene_length=TRUE){
  #Creates count matrix from folder of read_counter output files considering all species that were hit with at least <min_genes> marker genes
  #normalization: Assumes read_counter was run with "-y insert.raw_counts"; if norm_to_gene_length == TRUE, normalization is performed via:
  #1) sum up counts per species over all marker genes
  #2) divide by total length of marker genes for the given species
  #3) multiply by average length of marker genes over all species in the DB
  
  stopifnot(is_tibble(marker_genes.df))
  
  file_list <- list.files(path=path_to_folder)
  message(paste0(length(file_list)," files in the given folder"))
  ### initialize output df
  counts.df <- tibble(genus = character(0))
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
                                                              file.name = file.name,
                                                              min_genes = min_genes,
                                                              marker_genes.df = marker_genes.df,
                                                              tax.level = tax.level,
                                                              norm_to_gene_length=norm_to_gene_length))
    counts.df <- suppressMessages(full_join(counts.df,c.counts))
    pb$tick()
  }
  
  #create output matrix
  counts.mat <- 
    counts.df %>% 
    mutate(!!as.symbol(tax.level) := str_remove(!!as.symbol(tax.level),pattern = ".__")) %>% 
    as.data.frame() %>% 
    column_to_rownames(tax.level) %>% 
    replace(is.na(.),0) %>% 
    as.matrix()
  
  return(counts.mat)
}

.f_read_counter_to_count_mat <- function(read_counter_out.file,file.name,min_genes,marker_genes.df,tax.level="genus",norm_to_gene_length=TRUE){
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
  bact.counts <- tibble(!!as.symbol(tax.level) := "Bacteria", !!as.symbol(file.name) := sum(c.df$read_count))
  
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
  
  #Select the taxonomic level of interest and sum up all the (normalized)counts up to the selected level
  c.tax.counts <-
    c.df %>%
    select(!!as.symbol(tax.level),species) %>%
    inner_join(.,c.by_species %>% select(species,avg)) %>%
    distinct() %>%
    group_by(!!as.symbol(tax.level)) %>%
    #summarise(!!file.name := round(sum(avg),0)) %>%  #Do not round normalized counts since this could lead to false-zeros
    summarise(!!file.name := sum(avg)) %>% 
    bind_rows(bact.counts,.)
  
  
  
  return(c.tax.counts)
}

.f_ASV_to_tax <- function(lvl,asv_table,mseq.output){
  #Function that converts an asv-count matrix (e.g. from a dada2 result) to a count matrix at selected tax-level (lvl), 
  #by taking the taxonomic infroamtion from the mseq.output file (ASVs must be in column 1 and tax-tree must be in column 14 (-> default when running mapseq with --simple)))
  
  stopifnot(lvl %in% c('kingdom', 'phylum', 'class', 
                       'order', 'family', 'genus', 'species'))
  tax <- .f_resolve_ASVs(mseq.output = mseq.output)
  tax
  tax <- tax %>% 
    filter(ASV %in% rownames(asv_table))
  if (any(colSums(asv_table) == 0)){
    asv_table <- asv_table[,colSums(asv_table) > 0]
  }
  
  groups <- tax %>% 
    filter(!is.na(!!as.symbol(lvl))) %>% 
    pull(!!as.symbol(lvl)) %>% 
    unique
  groups
  mat.new <- matrix(0, nrow=length(groups), ncol=ncol(asv_table),
                    dimnames = list(groups, colnames(asv_table)))
  # Selects all ASVs that correspond to the current group (e.g current genera) and sums the counts
  for (g in groups){
    mat.new[g,] <- colSums(asv_table[
      tax %>% 
        filter(!!as.symbol(lvl)==g) %>% 
        pull(ASV), ,drop=FALSE])
  }
  #Add counts of samples that are not resolved at given tax_level
  #This is represented by the difference of colSums
  not_resolved <- as.matrix(t(colSums(asv_table)-colSums(mat.new)))
  rownames(not_resolved) <- "not_resolved"
  mat.new <- rbind(not_resolved,mat.new)
  
  #remove the underscore with the tax_level identifier (e.g. "g__")
  rownames(mat.new) <- str_remove(rownames(mat.new),pattern = "[a-zA-Z]__")
  
  print(summary(colSums(mat.new)))
  return(mat.new)
}








