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
      tmp.file <- read.table(paste0(path_to_folder,file_list[i]),sep = "\t",comment.char = "",skip = 1)
      colnames(tmp.file)[-1] <- tmp.file[2,-1]  
      tmp.file <- tmp.file[-2,]
      #colnames(tmp.file) <- gsub(colnames(tmp.file), pattern = "*_bac_ssu.mseq",replacement = "")
      colnames(tmp.file) <- str_extract(colnames(tmp.file), "[^_]+")
      tmp.file <- tmp.file %>% mutate_at(vars(-V1),as.numeric)
      
      count.df <- bind_rows(count.df,tmp.file)  
    }
    count.df <- as.data.frame(count.df %>% group_by(V1) %>% summarise_all(sum))
  }else{
    stop("no files found in the given directory")
  }
  
  #create output matrix (assign rownames and remove column with tax names)
  rownames(count.df) <- sub('.*\\;', '', count.df$V1)
  count.df$V1 <- NULL
  count.mat <- as.matrix(count.df)
  rownames(count.mat)[1] <- "Bacteria"
  #remove everything after the first point in the samplenames
  colnames(count.mat) <- sub(colnames(count.mat),pattern = "\\..*",replacement = "")
  return(count.mat)
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
  i <- 1
  for(i in seq(1,length(file_list))){
    
    ### read in file
    #chec if file is empty
    if(file.info(paste0(path_to_folder,file_list[i]))$size==0){
      message(paste0("Sample ",file_list[i]," is empty - skipping"))
      pb$tick()
      next
    }
    #c.f <- read_tsv(paste0(path_to_folder,file_list[i]),col_names = var.names,col_types = cols())
    c.f <-read.csv(paste0(path_to_folder,file_list[i]),header = FALSE,sep = "\t",comment.char = "",check.names = FALSE,col.names = var.names)
    c.f$tax.name <- str_replace_all(c.f$tax.name,pattern = " ",replacement = "")
    c.total.counts.mapped <- c.f %>% filter(tax.name == "root") %>% pull(counts.sum)
    total.counts.mapped[i,2] <- c.total.counts.mapped
    ### extract all bacterial counts
    bac.start <- which(c.f[,6] == "Bacteria")
     bac.end <- which(c.f[,6] == "Viruses")-1
    bac.reads <- c.f[(bac.start:bac.end),]  
    
    ### select counts and tax names
    c.counts.df <- bac.reads %>% filter(tax.symbol == tax.sym | tax.symbol == "D") %>% select(tax.name,counts.sum) %>% 
      rename(!!gsub(x = file_list[i],pattern = ".txt",replacement = "") := counts.sum)
    
    # Add to data from other samples
    counts.df <- suppressMessages(full_join(counts.df,c.counts.df,by="tax.name"))
    pb$tick()
  }
  ### convert to output matrix
  counts.df <- counts.df %>% mutate_at(vars(-tax.name),as.numeric)
  counts.mat <- as.matrix(counts.df[,-1])
  rownames(counts.mat) <- (counts.df$tax.name)
  counts.mat[is.na(counts.mat)] <- 0
  #remove everything after the first point in the samplenames
  colnames(counts.mat) <- sub(colnames(counts.mat),pattern = "\\..*",replacement = "")
  
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
  score.df <- tibble(tax.name = character(0))
  counts.df <- tibble(tax.name = character(0))

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
      select(name,score_normalized) %>%
      rename(!!gsub(x = file_list[i],pattern = ".pathseq.txt",replacement = "") := score_normalized,
             tax.name = name)

    # also unambiguous counts
    c.counts.df <-
      c.f %>% filter(type == tax.sym | type == "superkingdom",
                     kingdom == "Bacteria") %>%
      select(name,unambiguous) %>%
      rename(!!gsub(x = file_list[i],pattern = ".pathseq.txt",replacement = "") := unambiguous,
             tax.name = name)

    #Join with data from other samples
    score.df <- full_join(score.df,c.score.df,by="tax.name")
    counts.df <- full_join(counts.df,c.counts.df,by="tax.name")
    
    pb$tick()
  }
  
  ### convert to output matrices and return list
  score.df <- score.df %>% mutate_at(vars(-tax.name),as.numeric)
  score.mat <- as.matrix(score.df[,-1])
  rownames(score.mat) <- (score.df$tax.name)
  score.mat[is.na(score.mat)] <- 0
  #remove all-zero rows
  score.mat <- score.mat[rowSums(score.mat)>0,]
  #remove everything after the first point in the samplenames
  colnames(score.mat) <- sub(colnames(score.mat),pattern = "\\..*",replacement = "")
  
  counts.df <- counts.df %>% mutate_at(vars(-tax.name),as.numeric)
  counts.mat <- as.matrix(counts.df[,-1])
  rownames(counts.mat) <- (counts.df$tax.name)
  counts.mat[is.na(counts.mat)] <- 0
  counts.mat <- counts.mat[rowSums(counts.mat)>0,]
  #remove everything after the first point in the samplenames
  colnames(counts.mat) <- sub(colnames(counts.mat),pattern = "\\..*",replacement = "")
  
  
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
  count.mat <- count.mat[rowSums(count.mat)>0,]
  #remove everything after the first point in the samplenames
  colnames(count.mat) <- sub(colnames(count.mat),pattern = "\\..*",replacement = "")
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
  colnames(count.mat) <- sub(colnames(count.mat),pattern = "\\..*",replacement = "")
  return(count.mat)
}


.f_read_in_libsize <- function(path_to_folder){
  file_list <- list.files(path_to_folder)
  
  df.libsize <- tibble()
  for(i in seq(1,length(file_list))){
    tmp.file <- read.table(paste0(path_to_folder,file_list[i]),sep = "\t",comment.char = "",header = F)
    tmp.file$ID <- sub(file_list[i],pattern = "\\..*",replacement = "")
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
    tmp.file$ID <- sub(file_list[i],pattern = "\\..*",replacement = "")
    colnames(tmp.file)[1] <- "lib_layout"
    df.lib_layout <- bind_rows(df.lib_layout,tmp.file)
  }
  df.lib_layout <- df.lib_layout %>% relocate(ID)
  return(df.lib_layout)
}




