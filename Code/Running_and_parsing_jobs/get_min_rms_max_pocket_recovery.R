# This scripts load data from the parsed files (FlexPepDock results, RMSD-s calculated with PyMol and binding pocket recovery results.
# It concatenates them, allows for filtering for seed number, model number, recycles and drop-out (training) modes. 
# Input are directories per set, whether to process FPD, PyMol and pocket results and the possible filtering values.

library(plyr)
library(dplyr)
library(stringr)
library(optparse)
library(rlist)
library(reshape2)

option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default=NULL, 
              help="input directory on dataset level", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL, 
              help="output file name", metavar="character"),
  make_option(c("-c", "--column_list"), type="character", default=NULL, 
              help="list of columns to summarize", metavar="character"),
  make_option(c("--fpd"), default=F, action='store_true',
              help="read in FlexPepDock rescoring", metavar="character"),
  make_option(c("--pymol"), default=F, action='store_true',
              help="read in metrics calculated by pymol", metavar="character"),
  make_option(c("--pocket"), default=F, action='store_true',
              help="read in metrics calculated for overlapping interfaces", 
              metavar="character"),
  make_option(c("-f", '--factors'), type='character', default=NULL,
              help="factors to be extracted from filename (e.g recycle, seed, etc.), joined by comma", 
              metavar="character"),
  make_option(c("-s", '--seed'), type="character", default=NULL,
              help="only gather from models with this seed", metavar="character"),
  make_option(c("-m", '--model'), type="character", default=NULL,
              help="only gather from models with this model number", metavar="character"),
  make_option(c("-r", '--recycle'), type="character", default=NULL,
              help="only gather from models with this number of recycles", metavar="character"),
  make_option(c("-t", '--training'), default=FALSE, action='store_true',
              help="look for training/production in filename as additional factor?", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)

# NA dataset in the output means that you should add the name of the dataset here
map_datasets = setNames(c("Motif","Non-motif", 'AutoPeptiDB', 'Non-motif', 'Non-motif', 
                          "Motif","Non-motif", 'improved_dataset', 'approved'), 
                        c("Motif","NonMotif", 'non_redundant', 'run_full_peptide', 
                          "Non-motif", "pfpd_motif","pfpd_nonmotif", 'improved_dataset', 'approved'))

# function for reading results of FlexPepDcok
read_fpd_csv <- function(filename){
  pdb_id <- str_extract(filename, '([0-9][0-9a-zA-Z]{3})_[A-Za-z0-9]{2}')
  pfpd_csv <- read.csv(filename, header=T, stringsAsFactors=F, sep='', skip=1)
  pfpd_csv$pdb_id <- pdb_id
  
  return(pfpd_csv)
}


root <- '' # Root directory, above the directory per dataset
dataset_dirs <- Sys.glob(paste0(root, opt$input_dir))
factors <- unlist(stringr::str_split(opt$factors, ','))

# which columns should matter when merging the different tables
columns_to_merge <- c('pdb_id', 'dataset', 'model', 'rank')
columns_to_merge <- list.append(columns_to_merge, factors)

# which factors should be used to group by for gathering the results
columns_to_group_by <- c('pdb_id', 'dataset')
columns_to_group_by <- list.append(columns_to_group_by, factors)

if (opt$training == TRUE){
  print('training added')
  columns_to_merge <- list.append(columns_to_merge, 'training')
  columns_to_group_by <- list.append(columns_to_group_by, 'training')
}

# load data calculated by PyMOL
if (opt$pymol){
  all_data_list <- list()
  for (dir in dataset_dirs){
    dataset <- basename(dir)
    print(dir)
    if (dataset != 'binders_non_binders' & dir.exists(dir) & !startsWith(dataset, "interpep")){
      
      csvs <- Sys.glob(paste0(dir, '/*.csv'))
      data <- ldply(csvs, read.csv, header=TRUE, stringsAsFactors=F)
      
      if (nrow(data) > 0){
        data$dataset <- dataset
        all_data_list <- list.append(all_data_list, data)
      }
    }  
  }  
  
  all_data_pymol <- bind_rows(all_data_list)

  # annotate each model
  all_data_pymol$dataset = map_datasets[all_data_pymol$dataset]
  all_data_pymol$model <- str_extract(all_data_pymol$model_name, "model_[0-9]")
  all_data_pymol$rank <- str_extract(all_data_pymol$model_name, 'rank_[0-9]')  
  all_data_pymol$training <- grepl('training', all_data_pymol$model_name, fixed = TRUE)
  
  for (factor in factors){
    all_data_pymol[,factor] <- str_extract(all_data_pymol$model_name, paste0(factor, "_[0-9]"))
  }
  
}


# load data calculated by FlexPepDock
if (opt$fpd){
  all_data_list <- list()
  for (dir in dataset_dirs){
    dataset <- basename(dir)
    
      csvs <- Sys.glob(paste0(dir, '/*/*.score.sc'))
      data <- ldply(csvs, read_fpd_csv)
      
      if (nrow(data) > 0){
        data$dataset <- dataset
        all_data_list <- list.append(all_data_list, data)
    }  
  }  

  all_data_fpd <- bind_rows(all_data_list)
  
  # annotate each model  
  all_data_fpd$dataset = map_datasets[all_data_fpd$dataset]
  all_data_fpd$model <- str_extract(all_data_fpd$description, "model_[0-9]")
  all_data_fpd$rank <- str_extract(all_data_fpd$description, 'rank_[0-9]')
  all_data_fpd$training <- grepl('training', all_data_fpd$description, fixed = TRUE)

  for (factor in factors){
    all_data_fpd[,factor] <- str_extract(all_data_fpd$description, paste0(factor, "_[0-9]"))
  }
  
  # select longest consecutive stretch
  bestrms <- all_data_fpd %>% select(pdb_id, dataset, matches("bestRMS"))
  bestrms$bestRMS_0mer_all <- TRUE # add 0, to sum those that has no Xmer in these distances
  melt_aln_data <- melt(bestrms, id.vars = c('dataset', 'pdb_id'))
  
  melt_aln_data$in1 <- melt_aln_data$value <= 1 # is it in 1 A?
  melt_aln_data$in2.5 <- melt_aln_data$value <= 2.5 # is it in 2 A?
  melt_aln_data$in4 <- melt_aln_data$value <= 4 # is it in 4 A?
  melt_aln_data$aln_len <- as.integer(stringr::str_match(melt_aln_data$variable, '[0-9]+'))

  aln_data_1 <- data.frame(melt_aln_data %>% 
                             filter(in1) %>% 
                             group_by(dataset, pdb_id) %>% 
                             summarize_at(vars(max_in1A = aln_len), max))
  
  aln_data_2.5 <- data.frame(melt_aln_data %>% 
                               filter(in2.5) %>% 
                               group_by(dataset, pdb_id) %>% 
                               summarize_at(vars(max_in2.5A = aln_len), max))
  
  aln_data_4 <- data.frame(melt_aln_data %>% 
                             filter(in4) %>% 
                             group_by(dataset, pdb_id) %>% 
                             summarize_at(vars(max_in4A = aln_len), max))
  
  aln_data <- Reduce(function(x, y) merge(x, y, all=TRUE, by=c('dataset', 'pdb_id')), 
         list(aln_data_1, aln_data_2.5, aln_data_4))
  
}

# Load pocket recovery calculated by PyMOL
if (opt$pocket){
  all_data_list <- list()
  for (dir in dataset_dirs){
    dataset <- basename(dir)
    
    if(dataset != 'binders_non_binders' & dir.exists(dir) & !startsWith(dataset, "interpep")){
      
      csvs <- Sys.glob(paste0(dir, '/*/*_overlapping_interface.csv'))
      data <- ldply(csvs, read.csv, header=TRUE, stringsAsFactors=F)
      
      if (nrow(data) > 0){
        data$dataset <- dataset
        all_data_list <- list.append(all_data_list, data)
      }
    }  
  }  
  
  all_data_pocket <- bind_rows(all_data_list)
  colnames(all_data_pocket)[2] <- 'pdb_id'
  
  # annotate each model  
  all_data_pocket$dataset <- map_datasets[all_data_pocket$dataset]
  all_data_pocket$pdb_id <- str_extract(all_data_pocket$model_name, "[0-9][0-9a-z][0-9a-z][0-9a-z]_[0-9a-zA-Z][0-9a-zA-Z]")
  all_data_pocket$model <- str_extract(all_data_pocket$model_name, "model_[0-9]")
  all_data_pocket$rank <- str_extract(all_data_pocket$model_name, 'rank_[0-9]')
  all_data_pocket$training <- grepl('training', all_data_pocket$model_name, fixed = TRUE)
  
  for ( factor in factors ){
    all_data_pocket[,factor] <- str_extract(all_data_pocket$model_name, paste0(factor, "_[0-9]"))
  }
  
  if ( !is.null(opt$seed) ){
    all_data_pocket <- all_data_pocket %>% filter(across(seed, ~ grepl(paste0("seed_", opt$seed), .)))
  }
  
  if ( !is.null(opt$recycle) ){
    all_data_pocket <- all_data_pocket %>% filter(across(recycle, ~ grepl(paste0("recycle_", opt$recycle), .)))
  }
}

# Merge datasets, depending on what is the input
if ( opt$fpd & opt$pymol & opt$pocket ){
  all_data <- Reduce(function(x, y) merge(x, y, all=TRUE, 
                                          by=unlist(c(columns_to_merge, 'training'))), 
                     list(all_data_fpd, all_data_pymol, all_data_pocket))
}else if( opt$fpd & opt$pymol ){
  all_data <- Reduce(function(x, y) merge(x, y, all=TRUE, 
                                          by=unlist(c(columns_to_merge, 'training'))), 
                     list(all_data_fpd, all_data_pymol))
}else if( opt$fpd ){
  all_data <- all_data_fpd
}else if( opt$pymol ){
  all_data <- all_data_pymol
}else if( opt$pocket ){
  all_data <- all_data_pocket
}

# Filtering by given criteria. 
# Filter for 'production' runs if training is not a grouping variable
if (!is.null(opt$seed)){ # filter by seed
  all_data <- all_data %>% filter(seed==paste0("seed_", opt$seed))
}

if (!is.null(opt$recycle)){ # filter by recycle
  all_data <- all_data %>% filter(recycle==paste0("recycle_", opt$recycle))
}

if(!opt$training){ # added for filtering by training/production
  all_data <- all_data %>% filter(training==FALSE) %>% select(-training)
}

if(!is.null(opt$model)){ # added for filtering by model number
  print(str_split(opt$model,','))
  models <- paste0('model_', as.vector(str_split(opt$model,',')[[1]]))
  all_data <- all_data %>% filter(model %in% models)
}

write.csv(all_data, str_replace(opt$outfile, 'min_', 'all_')) 

columns_to_group_by <- unlist(columns_to_group_by)
column_list <- as.vector(unlist(str_split(opt$column_list, ',')))

# Also return training column if GROUPING happened by it
if (opt$training == TRUE){
  columns_to_return <- unlist(c(columns_to_merge, 'training', column_list))
}else{
  columns_to_return <- unlist(c(columns_to_merge, column_list))
}

# Summarize minimum values
all_data_min <- data.frame(all_data %>% 
                   group_by(!!!syms(columns_to_group_by)) %>% 
                   summarize_at(vars(any_of(column_list)), min))

# Summarize maximum values for pocket RMSDs
if(opt$pocket){
  all_data_max <- data.frame(all_data_pocket %>% 
                               group_by(!!!syms(columns_to_group_by)) %>% 
                               summarize_at('common_residues_percent', max))
  all_data_out <- merge(all_data_min %>% select(-common_residues_percent), 
                        all_data_max, by=c('pdb_id', 'dataset', 'seed', 'recycle'), all.x=T)
}else{
  all_data_out <- all_data_min
}

# merge bestXmer data if FPD was read in
if(opt$fpd){
  all_data_out <- merge(all_data_out, aln_data, by=c('pdb_id', 'dataset'), all.x=T)
}

# write output
write.csv(all_data_out, opt$outfile)

