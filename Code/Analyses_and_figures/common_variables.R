library(dplyr)
library(plyr)
library(ggplot2)
library(PupillometryR)
library(reshape2)
library(purrr)
library(stringr)
library(optparse)
library(rlist)
library(scales)

# colors
# dark blue: #04617B
# light blue: #009DD9
# red: #C00000
# ref fill: #ff7272
# yellow: #ffbf00
# yellow fill: #ffe69d
# green: #A5C249

### DATA ###
colnames_to_eval <- c('rms_before_ref_align_seq_rec_bb', 'rms_before_ref_align_seq_rec',
                      'rms_before_ref_align_seq_pep_bb', 'rms_before_ref_align_seq_pep', 
                      'rmsBB_if', 'rmsALL_if', 
                      'common_residues_percent',
                      'rmsBB_CAPRI_if', 'rmsALL_CAPRI_if', 
                      'rmsBB')

colnames_to_eval_pretty <- c('Receptor structure_bb', 'Receptor structure',
                             'Peptide structure_bb', 'Peptide structure',
                             "Peptide interface_bb",  
                             "Peptide interface", 
                             'Binding residue recovery',
                             'CAPRI interface_bb', 'CAPRI interface',
                             'Peptide_bb')

rms_mapping <- setNames(colnames_to_eval_pretty, colnames_to_eval)

secstr <- c('pep_is_helix', 'pep_is_strand', 'pep_is_combined', 'pep_is_coil',
            'rec_is_helix', 'rec_is_strand', 'rec_is_combined')

secstr_pretty <- c('Helix', 'Strand',
                   'Combined', 'Coil',
                   'Helix', 'Strand',
                   'Combined')

secstr_mapping <- setNames(secstr_pretty, secstr)

breaks <- seq(0,20,0.5)

### FUNCTIONS ###
calculate_binned_freq <- function(column, breaks, all_count, direction='forward'){
  results <- list()
  for (br in breaks){
    if (direction == 'forward'){
      count <- sum(na.omit(column) <= br) / all_count
    }else{
      count <- sum(na.omit(column) >= br) / all_count
    }
    results <- append(results, count)
  }
  return(unlist(results))
}


### THEMES ### 
new_theme <- theme_minimal()+ theme(aspect.ratio=1,
                                    legend.position="none",
                                    panel.background = element_rect(colour = "darkgray", size=1),
                                    panel.grid.major = element_line(colour="gray85", size=0.5),
                                    strip.text.x = element_text(size = 12),
                                    strip.text.y = element_text(size = 12, angle=-90),
                                    strip.background = element_rect(fill="white", color='white'))
theme_set(new_theme)


raincloud_theme = theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(vjust = 0.5),
  legend.title = element_text(size=16),
  legend.text = element_text(size=16),
  legend.position = "right",
  plot.title = element_text(lineheight=.8, face="bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
  panel.background = element_rect(colour = "darkgray", size=1),
  strip.text.x = element_text(size = 16),
  strip.text.y = element_text(size = 16, angle=-90),
  strip.background = element_rect(fill="white", color='white')
)

colorscale = c("Motif" = "#009DD9", "Non-motif" = '#C00000', 'approved' = "#ffbf00", 'NA' = 'grey30', "motif" = "#009DD9", "nonmotif" = '#C00000')
color_scale_for_raincloud = c("Motif" = '#04617B', "Non-motif" = '#C00000', 'approved' = "#FFBF00")
fill_scale_for_raincloud = c("Motif" = '#009DD9', "Non-motif" = '#ff7272', 'approved' = "#ffe69d")

