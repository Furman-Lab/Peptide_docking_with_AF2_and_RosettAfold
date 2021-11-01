if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
  BiocManager::install("ComplexHeatmap")

library("data.table")
library("ComplexHeatmap")
library(circlize)

# args = commandArgs(trailingOnly=TRUE)
process_row_names <- function(row_name_vec, type_to_use){
  split_vec = strsplit(row_name_vec,split="_")
  print("split_ved")
  print(split_vec)
  
  parsed_vec = "A"
  for(i in 1:length(split_vec)){
    if(type_to_use == "best_lddt"){

      tmp_name = paste0(split_vec[[i]][1:2],collapse = "_")
    }else if(type_to_use == "network_params"){
      tmp_name = paste0(split_vec[[i]][3:4],collapse = "_")
    }
   parsed_vec = c(parsed_vec, tmp_name) 
  }
  return (parsed_vec[-1])
}

bind_with_second_run_type <- function(first_mat, i_prot, metric){
  list_of_files = list.files()[grep(metric, list.files())]
  list_of_files = list_of_files[-grep("comparison",list_of_files)]
  list_of_files = list_of_files[grep(".tsv",list_of_files)]
  list_of_files = list_of_files[grep("sep",list_of_files)]
  list_of_files = list_of_files[grep(i_prot,list_of_files)]

  i_tab = fread(list_of_files[1])
  row_names = i_tab$V1

  i_mat = as.matrix(i_tab[,-1])
  rownames(i_mat) = paste0("sep_",row_names)
  order_of_rows = paste0("sep_",c("model_1","model_2","model_3","model_4","model_5"))

  i_mat = i_mat[order_of_rows,]
  
  print("i_mat")
  print(i_mat)
  print("first_mat")
  print(first_mat)
  
  bound_mat = rbind(first_mat, i_mat)
  return(bound_mat)
}
# 


draw_by_residue_heatmaps <- function(set, metric, value_vec, color_vec, path, how_to_order){
  setwd(path)
  heatmap_list = list()
  
  motif_tab = fread("/vol/ek/share/peptide_docking_with_afold2_and_rosettAfold/motif_analyses/first_submission/motif_peptides_with_motif_annotation.tsv")
  list_of_files = list.files()[grep(metric, list.files())]
  list_of_files = list_of_files[-grep("comparison",list_of_files)]
  list_of_files = list_of_files[grep(".tsv",list_of_files)]
  list_of_files = list_of_files[grep("linker",list_of_files)]
  if(set == "non_motif"){
    list_of_files = list_of_files[-grep("1lvm",list_of_files)]
  }
  # 
  # print("list of files")
  # print(list_of_files)
  # 
  for (i_file in list_of_files){
    # i_file="1jwg_by_residue_rms.tsv"
    i_tab = fread(i_file)
    i_prot = unlist(strsplit(i_file,split="_",fixed=F))[2]
    print(i_prot)
    
    i_motif = motif_tab$motif[which(motif_tab$pdb_id==i_prot)]
    
    motif_bool = unlist(strsplit(i_motif,split="")) %in% LETTERS
    
    row_names = i_tab$V1

    i_mat = as.matrix(i_tab[,-1])
    type_to_sort_rownames_by = how_to_order
    
    ## use this only when running on the first submission data, where pdbs had names with "model_#_model_#"
    ## for params and for rank
    ##
    # processed_rownames = process_row_names(row_names, type_to_sort_rownames_by)
    # rownames(i_mat) = processed_rownames
    ##
    
    rownames(i_mat) = paste0("linker_",row_names)
    order_of_rows = paste0("linker_",c("model_1","model_2","model_3","model_4","model_5"))
    # colnames(i_mat)<-paste("",colnames(i_mat))
    
    i_mat = i_mat[order_of_rows,]
    i_mat = bind_with_second_run_type(i_mat, i_prot, metric)
    
    if(metric == "lddt"){
      i_mat = i_mat/100
    }
    
    order_of_rows_combined = c(order_of_rows, paste0("sep_",c("model_1","model_2","model_3","model_4","model_5")))
    
    i_ha_name <- paste0(i_prot, "_ha")
    
    print(i_mat)
    if(set=="motif"){
      tmp_ha = HeatmapAnnotation(motif =motif_bool, col = list(motif = c("FALSE"="gray70","TRUE"="navy")), show_annotation_name = FALSE, show_legend = FALSE)#,
                               #heatmap_legend_param = list(direction = "horizontal"))
    }else{
      tmp_ha=NULL
    }
    
    colanno = HeatmapAnnotation(cn = anno_text(colnames(i_mat), just = "bottom", location = unit(0, "npc")))
                           #ha_height = max_text_height(colnames(mat)))
    
    col_fun = colorRamp2(value_vec, color_vec)
    col_fun(seq(-3, 3))
    
    i_heatmap_name <- paste0(i_prot, "_heatmap")
    tmp_ht = Heatmap(i_mat, name = i_prot, row_order = order_of_rows_combined, cluster_columns = FALSE, top_annotation = tmp_ha,
                     cluster_rows = FALSE,
                     column_names_rot = 0,
                     #bottom_annotation = colanno,
                     
                     show_heatmap_legend = FALSE,
                     border_gp = gpar(col = "black", lty = 1),
                     show_column_names = TRUE,
                     column_title = i_prot,
                     heatmap_legend_param = list(direction = "horizontal"),
                     #row_names_gp = gpar(fontsize = 10),
                     col=col_fun, show_row_names =FALSE) 
    
    
    
    if(exists("all_ht")){
      all_ht = all_ht+tmp_ht
    }else{
      all_ht = tmp_ht
    }
    # draw(tmp_ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
    # heatmap_list = heatmap_list + get(i_heatmap_name)
    
    
    
  }
  draw(all_ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", row_split = c("linker","linker","linker","linker","linker","sep","sep","sep","sep","sep"))
  return(all_ht)
  
  
}

motif_path="/vol/ek/share/peptide_docking_with_afold2_and_rosettAfold/motif_analyses/review/motif_by_residue_rms"
non_motif_path="/vol/ek/share/peptide_docking_with_afold2_and_rosettAfold/motif_analyses/review/non_motif_by_residue_rms"

how_to_order = "network_params"
type_row_split = c("linker","linker","linker","linker","linker","sep","sep","sep","sep","sep")


motif_lddt_all_ht = draw_by_residue_heatmaps("motif", "lddt", c(0, 0.5,  1), c("coral4","tan3","white"), motif_path, how_to_order)
non_motif_lddt_all_ht = draw_by_residue_heatmaps("non_motif", "lddt", c(0, 0.5,  1), c("coral4","tan3","white"), non_motif_path, how_to_order)
 
motif_rms_all_ht = draw_by_residue_heatmaps("motif", "rms", c(-1, 0, 2.5,  5), c("black","ivory","palegreen","darkgreen"), motif_path, how_to_order)
non_motif_rms_all_ht = draw_by_residue_heatmaps("non_motif", "rms",  c(-1, 0, 2.5,  5), c("black","ivory","palegreen","darkgreen"), non_motif_path, how_to_order)

 
png(file=paste0("review_",how_to_order,"_ord_motif_rms_heatmaps.png"),width = 5000, height = 1000, res = 300)
draw(motif_rms_all_ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom" , auto_adjust=FALSE, row_split = type_row_split)
dev.off()

png(file=paste0("review_",how_to_order,"_ord_motif_lddt_heatmaps.png"),width = 5000, height = 1000, res = 300)
draw(motif_lddt_all_ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom" , auto_adjust=FALSE, row_split = type_row_split)
dev.off()


png(file=paste0("review_",how_to_order,"_ord_non_motif_lddt_heatmaps.png"),width = 5000, height = 1000, res = 300)
draw(non_motif_lddt_all_ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom" , auto_adjust=FALSE, row_split = type_row_split)
dev.off()


png(file=paste0("review_",how_to_order,"_ord_non_motif_rms_heatmaps.png"),width = 5000, height = 1000, res = 300)
draw(non_motif_rms_all_ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom" , auto_adjust=FALSE, row_split = type_row_split)
dev.off()



## color backup
## motif_rms_all_ht = draw_by_residue_heatmaps("motif", "rms", c(-1, 0, 2.5,  5), c("black","ivory","goldenrod1","red3"), motif_path, how_to_order)
## non_motif_lddt_all_ht = draw_by_residue_heatmaps("non_motif", "lddt", c(0, 0.5,  1), c("blue4","lightskyblue","white"), non_motif_path, how_to_order)

