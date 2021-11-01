##########################################################
# load common variables like color, mapping, etc and libraries
getwd()
source('common_variables.R')
##########################################################

# load data
motif_linker_data <- read.csv('../../Data/min_linker_env_motif_new_columns.csv', stringsAsFactors = F)
nonmotif_linker_data <- read.csv('../../Data/min_linker_env_nonmotif_new_columns.csv', stringsAsFactors = F)
approved_linker_data <- read.csv('../../Data/min_linker_env_improved_dataset_new_columns.csv', stringsAsFactors = F)


motif_sep_data <- read.csv('../../Data/min_sep_chains_env_motif_new_columns.csv', stringsAsFactors = F)
nonmotif_sep_data <- read.csv('../../Data/min_sep_chains_env_nonmotif_new_columns.csv', stringsAsFactors = F)
approved_sep_data <- read.csv('../../Data/min_sep_chains_env_improved_dataset_new_columns.csv', stringsAsFactors = F)


all_fig_1_linker_data <- do.call("rbind", list(motif_linker_data, nonmotif_linker_data, approved_linker_data))
all_fig_1_linker_data$link_type <- 'linker'

all_fig_1_sep_data <- do.call("rbind", list(motif_sep_data, nonmotif_sep_data, approved_sep_data))
all_fig_1_sep_data$link_type <- 'sep_chains'

all_rms_data <- rbind(all_fig_1_linker_data, all_fig_1_sep_data) %>% filter(pdb_id!='1lvm_AE')

final_pfpd_list <- unique(substr(all_rms_data %>%
                                   filter(dataset!='approved') %>% pull(pdb_id), 1, 4))

# gatherin minium data points
#all_rms_data <- all_rms_data %>% select(-common_residues_percent)

# calculate binned frequency
list_of_dfs <- list()
for (dataset in unique(all_rms_data$dataset)){
  for (link_type in unique(all_rms_data$link_type)){
    print(paste(dataset, link_type))
    
    small_data <- all_rms_data %>% filter(dataset==!!dataset) %>% filter(link_type==!!link_type)
    
    for (col in colnames_to_eval[!grepl('common', colnames_to_eval)]){
      column_to_eval <- small_data %>% pull(!!col)
      frequencies <- data.frame(calculate_binned_freq(column_to_eval, breaks, length(column_to_eval)))
      colnames(frequencies) <- c('frequency')
      frequencies$rms_type <- col
      frequencies$dataset <- dataset
      frequencies$link_type <- link_type
      frequencies$breaks <- breaks
      
      list_of_dfs <- list.append(list_of_dfs, frequencies)
    }
  }
}

all_frequencies <- do.call("rbind", list_of_dfs)
all_frequencies$rms_type <- rms_mapping[all_frequencies$rms_type]
all_frequencies$rms_atoms <- 'All atom'
all_frequencies$rms_atoms[grepl('_bb$', all_frequencies$rms_type)] <- 'Backbone'
all_frequencies$rms_type <- gsub('_bb', '', all_frequencies$rms_type)
all_frequencies$rms_atoms <- factor(all_frequencies$rms_atoms, 
                                    levels=c('Backbone','All atom'))

capri_all_frequencies <- all_frequencies[grepl('CAPRI', all_frequencies$rms_type),]
all_frequencies <- all_frequencies[!grepl('CAPRI', all_frequencies$rms_type),]

# combining linking types - this will go to all the plots
all_rms_data_combined <- data.frame(all_rms_data %>% 
                                    group_by(pdb_id, dataset) %>% 
                                    summarize_at(vars(all_of(colnames_to_eval)), min))

list_of_dfs_combined <- list()
for (dataset in unique(all_rms_data_combined$dataset)){
  small_data <- all_rms_data_combined %>% filter(dataset==!!dataset)
  print(dataset)
  
  for (col in colnames_to_eval[!grepl('common', colnames_to_eval)]){
    column_to_eval <- small_data %>% pull(!!col)
    frequencies <- data.frame(calculate_binned_freq(column_to_eval, breaks, length(column_to_eval)))
    colnames(frequencies) <- c('frequency')
    frequencies$rms_type <- col
    frequencies$dataset <- dataset
    frequencies$breaks <- breaks
    
    list_of_dfs_combined <- list.append(list_of_dfs_combined, frequencies)
  }
}

all_frequencies_combined <- do.call("rbind", list_of_dfs_combined)
all_frequencies_combined$rms_type <- rms_mapping[all_frequencies_combined$rms_type]
all_frequencies_combined$rms_atoms <- 'All atom'
all_frequencies_combined$rms_atoms[grepl('_bb$', all_frequencies_combined$rms_type)] <- 'Backbone'
all_frequencies_combined$rms_type <- gsub('_bb', '', all_frequencies_combined$rms_type)
all_frequencies_combined$rms_atoms <- factor(all_frequencies_combined$rms_atoms, 
                                             levels=c('Backbone','All atom'))

all_frequencies_combined$method <- 'AlphaFold'
all_frequencies_combined$link_type <- 'combined'
all_frequencies_combined$exp_type <- 'final'

capri_all_frequencies_combined <- all_frequencies_combined[grepl('CAPRI', all_frequencies_combined$rms_type),]
all_frequencies_combined <- all_frequencies_combined[!grepl('CAPRI', all_frequencies_combined$rms_type),]

capri_columns <- colnames(all_rms_data_combined)[grepl('CAPRI', colnames(all_rms_data_combined))]
capri_all_rms_combined <- all_rms_data_combined  %>% select(pdb_id, dataset, !!capri_columns)
all_rms_data_combined <- all_rms_data_combined %>% select(-c(!!capri_columns))
all_rms_data_combined$method <- 'AlphaFold'

##########################################################
# PFPD data
Sys.glob('*')
pfpd <- read.csv('../../Data/min_values_pfpd.csv') %>% select(-X)
pfpd <- pfpd %>% filter(pdb %in% !!final_pfpd_list)

list_of_dfs <- list()
for (dataset in unique(pfpd$dataset)){
  print(paste(dataset))
  
  small_data <- pfpd %>% filter(dataset==!!dataset)
  
  for (col in c('rmsBB_if', 'rmsALL_if')){
    column_to_eval <- small_data %>% pull(!!col)
    frequencies_pfdp <- data.frame(calculate_binned_freq(column_to_eval, breaks, length(column_to_eval)))
    colnames(frequencies_pfdp) <- c('frequency')
    frequencies_pfdp$rms_type <- col
    frequencies_pfdp$dataset <- dataset
    frequencies_pfdp$breaks <- breaks
    
    list_of_dfs <- list.append(list_of_dfs, frequencies_pfdp)
  }
}

all_frequencies_pfpd <- do.call("rbind", list_of_dfs)
all_frequencies_pfpd$rms_type <- rms_mapping[all_frequencies_pfpd$rms_type]
all_frequencies_pfpd$rms_atoms <- 'All atom'
all_frequencies_pfpd$rms_atoms[grepl('_bb$', all_frequencies_pfpd$rms_type)] <- 'Backbone'
all_frequencies_pfpd$rms_type <- gsub('_bb', '', all_frequencies_pfpd$rms_type)
all_frequencies_pfpd$rms_atoms <- factor(all_frequencies_pfpd$rms_atoms, 
                                         levels=c('Backbone','All atom'))
all_frequencies_pfpd$method <- 'PFPD'

##########################################################
# Merge PFPD and AF data
all_frequencies_pfpd$link_type <- 'combined'
all_frequencies_pfpd$exp_type <- 'final'
all_frequencies_combined_pfpd_af <- rbind(all_frequencies_combined, all_frequencies_pfpd)

