# title: "buildconstellation"
# date: 2021-03-03
# output: html_document
# params:
  # cells.cl.df: NULL # s.obj@meta.data
  # groups.col: NULL  # character, name of column to group by
  # rd.dat: list()
  # colors.use : c()
  # cl: NULL
  # cl.numeric: NULL
  # n.pcs: 20,
  # run.knn: TRUE
  # k: 15,
  # frac.th: 0.05,
  # cl.df: NULL
  # out.dir: '../out'

buildConstellation <- function(args_list, 
                               groups.col, 
                               node_label,
                               colors.use,
                               subclades, 
                               run.knn, n.pcs, k, frac.th, 
                               out.dir, 
                               name = '', 
                               plot.parts = FALSE, 
                               dodge_value = 3, 
                               max_size = 25, label.size = 2, 
                               exaggeration = exaggeration) {

source_rmd('~/cse-phd/second-trimester/neocortex/constellation/R/scrattch.hicat_fxns.Rmd')
  
cells.cl.df <- args_list$cells_cl_df
rd.dat <- args_list$reductions


# 1. Build dataframes for constellation plots.

# 1.1 cl_and_cl_numeric
# cl and cl numeric can also be passed as a param to the Rmd.

if(!exists('args_list$cl') & !exists('args_list$cl.numeric')) {
  
    message('Getting cl and cl.numeric')
  cell_names <- cells.cl.df$cell_name
    # Cluster, or area-celltype combination.
    cl <- cells.cl.df[[groups.col]] %>% as.factor %>% magrittr::set_names(value = cell_names)

    cl.numeric <- as.numeric(cl) %>% magrittr::set_names(value = cell_names)
    
} else { message('Grabbing cl and cl.numeric from cells.cl.df')
        cl <- args_list$cl
        cl.numeric <- args_list$cl_numeric }

message(paste0('cl: ', head(cl), '\n',
               'cl.numeric: ', head(cl.numeric)))

## 2 cl.df

cl.df <- get_cl_df(cl)

cl.df$clade <- str_split_fixed(cl.df$cluster_label, "_", 2)[ ,1] %>% tolower 
cl.df$area <- str_split_fixed(cl.df$cluster_label, "_", 2)[ ,2] %>% tolower 

# Add clade_id, clade_color to cl.df
# cl.df <- cl.df %>% left_join(clade.cols)


cl.df <- cl.df %>% left_join(colors.use, by = c("clade" = "celltype")) 

# %>%  rename(cluster_color = "colors.use")
# rm(group.cols)

# cells.cl.df: Add cluster_id column from cl.df; remove unused columns. 
cells.cl.df <- cells.cl.df %>% dplyr::select(-cluster_label) %>% dplyr::rename(cluster_label = groups.col) %>%
                        left_join(
                         # %>% select(cell.name, groups.col, combined.cluster.2),
                         cl.df, by = "cluster_label") %>%
                         
                         # Requires cells.cl.df (metadata) to have column being used for groups
                         # named 'cluster_label' to match with cl_df during join.
                 mutate(cluster_id = as.factor(cluster_id))

message('Updated cells.cl.df')

# ----------------------------------------------
## 4 Find cluster centroids from UMAP coordinates
## rd.cl.center

rd.cl.center <- get_RD_cl_center(rd.dat = rd.dat$umap, cl)
message("Got rd.cl.center")

# update-rd.cl.center

rd.cl.center %<>% 
  as.data.frame %>% 
  set_names(c("x", "y")) %>%
  add_column(cl = cl.df$cluster_id, .before = "x") %>%
  # add_column preserves rownames.
  # but moving rownames to column cluster_label anyway bc of left_join below.
  # Needs to be cl (not cl_id) or else you get error:
  # Error in `$<-.data.frame`(`*tmp*`, "edge.frac.within", value = numeric(0)) : 
  # replacement has 0 rows, data has 26 
  rownames_to_column("cluster_label")

message("Updated rd.cl.center")

## 5 Join `cl.df` and `rd.cl.center` into `cl.center.df` for input into `get_KNN_graph`.
## cl.center.df
cl.center.df <- left_join(rd.cl.center, cl.df,
                          by = c("cluster_label")) 

## Update cl.center.df for radial glia and neuron subtypes (deep layer, upper layer, vRG, oRG) ------------------------

if(subclades) {

  updateClCenterDf <- function(cl.center.df) {
    
    message('Updating cl.center.df for subclades')
    
    areas <- c('pfc', 'motor', 'somatosensory', 'parietal', 'temporal', 'v1')
    areas_rgx <- paste0('_', areas) %>% glue_collapse('|')
    
    cl.center.df %<>% mutate_if(.predicate = is.character, tolower) %>%
      # mutate(cluster_label = str_extract(cluster_label, areas_rgx) %>% str_remove('_'),
      #       subclade = case_when(! .$clade == 'ipc' ~ str_remove(.$area, areas_rgx)),
      # ) %>%
      # unite(clade, 'clade', 'subclade') %>% 
      # mutate(clade = str_remove_all(clade, '_NA| |late')) %>%
      # select(-cluster_color) %>%
      left_join(cluster_colors, by = c('clade' = 'celltype'))
    
    message(cl.center.df)
    return(cl.center.df)
  }
  
  cl.center.df <- updateClCenterDf(cl.center.df)
}


# 6 Get knn and cluster counts
# Calls `knn.cl` in scrattch.hicat_fxns.Rmd
# knn.result
if(run.knn == TRUE) {
  
  message('Running RANN::nn2')

knn.result <- RANN::nn2(data = rd.dat$pca[, 1:n.pcs], k = k)

} else { knn.result <- knn.result }

knn.cl <- get_knn_graph(knn.result = knn.result,
                        rd.dat = rd.dat$umap, 
                        cl.df =  cl.df, 
                        cl = cl.numeric,
                        cl.numeric = cl.numeric,
                        knn.outlier.th = 2, 
                        outlier.frac.t = 0.5)

# rm(rd.dat, ncx.clusters)

# -----------------------------------------------------------------------------
# 2. Make constellation plot


# knn_cl_df_filter
message('Calculating knn.cl.df.filter')
# Keep only cells whose $frac >= 0.05.
# frac = fraction of cells in cluster with nearest neighbors in a different cluster.
# Defined in `get_knn_graph`: 
# knn.cl.df$frac = knn.cl.df$Freq / knn.cl.df$cl.from.total
# 10% : 213 edges

filterKNN <- function(knn.cl, frac.th = 0.1) {
knn.cl.df.filter <- knn.cl$knn.cl.df %>% dplyr::filter(frac >= frac.th) %>% 
  mutate(cl.from = as.numeric(cl.from), cl.to = as.numeric(cl.to))
}

knn.cl.df.filter <- filterKNN(knn.cl = knn.cl)

# cl.to, cl.from numeric or factor?
# Need to be numeric for getting the rows where cl.to == cl.from (knn.cl.df$same)





## plot-constellation

# To plot only edges between ExN lineage clusters:
# knn.cl.df %<>% filter_at(vars(cl.from.label, cl.to.label), 
#                        all_vars(str_detect(., "RG|IPC|Neuron|OPC|Dividing")))
st <- format(Sys.time(), "%Y%m%d.%H%M_")

message('Printing constellation plot.')

cl.plot <- plot_constellation(knn.cl.df = knn.cl.df.filter, 
                              cl.center.df = cl.center.df, 
                              out.dir = out.dir,
                              node.label = node_label,  # name of column in cl.center.df with cluster/node names
                              exxageration = exaggeration, curved = TRUE, 
                              plot.parts = plot.parts, plot.hull = NULL, 
                              plot.height = 15, plot.width = 8,
                              node.dodge = TRUE, dodge_value = dodge_value,
                              label.size = label.size, max_size = max_size)

results <- lst(cl, cl.numeric, cl.df, cl.center.df, 
               knn.result, knn.cl, knn.cl.df.filter, cl.plot)

message('Writing rds...')

write_rds(x = results, 
          path = file.path(out.dir, paste0(st, 'constellation_result_', name, '.rds')))

return(results)

}
