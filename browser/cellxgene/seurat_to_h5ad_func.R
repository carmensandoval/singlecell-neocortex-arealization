seurat_to_h5ad <- function(seurat_object_path = NULL, 
                           seurat_object = NULL, 
                           file_path_out) {
  
  if(!is.null(seurat_object_path)) {
    seurat_object <- read_rds(seurat_object_path)
  }
  
  SeuratDisk::SaveH5Seurat(object = seurat_object,  # Save H5Seurat adds .h5Seurat to filename
                           filename = file_path_out, 
                           overwrite = TRUE)
  
  SeuratDisk::Convert(source = file_path_out, 
                      dest = "h5ad", 
                      overwrite = TRUE)
  sms("Convert h5ad")
}

# seurat_to_h5ad(seurat_object = "~/cse-phd/second-trimester/neocortex/data/210322_neocortex_full_v3.2.3_subset_0.01.rds", 
#               file_path = "~/cse-phd/second-trimester/browser/browser-dbx/data/210323_neocortex_v3.2.3_subset.05.h5seurat")

# /////////////////////////////////////////////////////////////////////////////////////////////////

# For neocortex
cleanMetadataNeocortex <- function(seurat_object){
  
  neocortex_annot <- read_tsv("../neocortex/data/tbls/83d19_Neocortex_allindividuals_combinedclusters_v1.txt")
  gw16_annot <- read_tsv("../neocortex/data/tbls/gw16neo_clusteridentity.txt")
  
  subtypes_annot <- 
    read_tsv("../neocortex/data/tbls/Neocortex_subset1_clustermarkers_combo2_RGandNeuronannotations.csv") %>% 
    mutate_all( .funs = . %>% str_replace("^neuron", "excitatory_neuron") %>%
                  str_replace("rg", "radial_glia") %>% 
                  str_replace("-", "_"))
  
  
  annot <- rbind(neocortex_annot %>% select(cell_id = cell.name, 
                                            cluster_label = combined.cluster.2),
                 gw16_annot %>% unite(cluster_label, cell_type, cluster_id, sep = "_", remove = FALSE) %>%
                   select(cell_id = cell_name, 
                          cluster_label)) %>%
    mutate(cluster_label = str_remove_all(cluster_label, "_combo2") %>% tolower %>%
             str_replace("inteneuron", "interneuron") %>%
             str_replace("endo", "endothelial") %>%
             str_replace("rg|radial glia", "radial_glia") %>%
             str_replace("^neuron", "excitatory_neuron") %>%
             str_replace("^0|other|mixed", "outlier"),
           cell_type = cluster_label %>% str_remove_all("_[:digit:]+") %>%
             str_replace("^cr$", "cajal_retzius") 
    ) %>% left_join(subtypes_annot %>% select(-cell_type), by = "cluster_label") %>%
    mutate(cell_type = case_when(!is.na(annotation ) ~ annotation,
                                 TRUE ~ cell_type),
           # ******
           lamina = str_replace_all(lamina, c("All" = "all", 
                                              "VZSVZ" = "VZ_SVZ"))) %>% 
    rename(area = "cortical_area")
    select(-annotation, -ncount_rna, -nfeature_rna)
    # *****
  
  # 457,965 rows
  metadata_tmp <- seurat_object@meta.data %>% rownames_to_column("cell_id")

  metadata <- left_join(x = metadata_tmp,
                        y = annot, 
                        by = "cell_id") %>% 
    set_colnames(colnames(.) %>% tolower %>% str_replace("[.]", "_")) %>% 
    select(-orig_ident, -v1, -consensuscelltype, -area_sub) %>%
    select(cell_id, age, individual, brain_region = structure, cortical_area = area, 
           lamina, cell_type, cluster_label, everything()) %>%
    mutate_at(c("brain_region", "cortical_area"), .funs = tolower) %>%
    set_rownames(.$cell_id)
    
  cell_types_rgx <- metadata$cell_type %>% unique %>% glue::collapse("|") %>% 
    str_remove("outlier|") %>% paste0('(', ., ')') 
  
  metadata %<>% mutate(individual = case_when(str_detect(cell_id, "gw19_2") ~"GW19_2",
                              TRUE ~ individual),
       cluster_label = case_when(str_detect(cell_id, "GW16") ~ str_replace(string = cluster_label, 
                                                                           pattern = cell_types_rgx, 
                                                                           replacement = "\\1_16"),
                                 TRUE ~ cluster_label),
       
       cell_type = case_when(str_detect(cell_id, "GW16") ~ str_replace(string = cell_type,
                                                                       pattern = "(^[:alpha:]+_[:alpha:]+)_[:alpha:]+",
                                                                       replacement = "\\1"),
                              TRUE ~ cell_type))
  
  return(metadata)
}

# //////////////////////////////////

addReductionRows <- function(reduc_item, new_cells) {
  
  cell_embeddings <- reduc_item@cell.embeddings
  
  new_rows <- matrix(dimnames = list(row_names = new_cells, #gw16_cells, 
                                     col_names = colnames(cell_embeddings)), 
                     nrow = length(new_cells), 
                     ncol = length(colnames(cell_embeddings)), 
                     data = NaN)
  
  cell_embeddings <- rbind(cell_embeddings, new_rows)
  reduc_item@cell.embeddings <- cell_embeddings
  return(reduc_item)
  
}


# ////////////////////////////////////////////////

cleanMetadataWholeBrain <- function(seurat_object){
  
    tmp <- seurat_object@meta.data %>% 
    # rownames_to_column("cell_id") %>% 
    set_colnames(colnames(.) %>% tolower %>% str_replace_all("[.]", "_")) %>%
      
      select(cell_id = cell_name, age, individual, brain_region = structure, area, lamina, 
             cell_type = cell_type_v2, cluster_label = clusterv2, 
             ngene, numi, percent_mito, percent_ribo) %>%

    mutate_at(.vars = c("cluster_label", "cell_type"), 
              .funs = ~ tolower %>% str_replace_all(string = ., 
                                                    pattern = c("([:digit:]+)[:alpha:]+$" = "\\1",
                                                             "^neuron" = "excitatory_neuron", 
                                                             "radialglia|rg" = "radial_glia", 
                                                             "vascular|endo" = "endothelial"))) %>%
      
      mutate(area = case_when(area == "V1" & !(brain_region == "neocortex") ~ brain_region,
                              TRUE ~ area)) %>%

      set_rownames(.$"cell_id")
   
  return(tmp)
}