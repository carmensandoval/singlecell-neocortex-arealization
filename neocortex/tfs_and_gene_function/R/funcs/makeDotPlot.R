# 2020-11-25
# function: makeDotPlot

# Takes in a list of genes and a Seurat object,
# and retrieves the average gene expression for each group of cells, 
# as specified in group_by (and optionally in split_by)

# Paramteters:
# genes_list can be a vector of genes or a list of vectors.

makeDotPlot <- function(# Params for Seurat::DotPlot
                          seurat_object,  
                          genes_list,   
                          idents_use, 
                          group_by, 
                          split_by,  # cell_type
                          scale = FALSE,
                        # Params for scaling and plotting
                          mode = c('genes_by_area', 'paired_areas') ,
                                # paired-areas refers to genes by age in side by side areas.
                          scale_by = 'gene',
                        # Params for plotting and printing
                          color_scale = 'viridis',
                          max_size = 500/n_genes,
                          group, 
                          area, 
                          title,
                          quantiles = c(.05 , .95)) {

  # 1. Make Seurat DotPlot -------------------------
  
  dotplot <- Seurat::DotPlot(object = seurat_object,
                           idents = idents_use,
                           features = genes_list,
                           group.by = group_by,
                           split.by = split_by,
                           scale = scale,
                           cols =  "RdYlBu",
                           dot.min = 0.01)
        
  # 2. Use dotplot$data to get average gene expression for each group of cells -------------------

      x <- dotplot$data %>%
          separate(id, c('id_1', 'id_2'), remove = FALSE)
  
      print(x)
      ## Clean up x$data for manipulations below (scaling, ggplot)

      # For normal DotPlot  
        # id_1 = area
        # id_2 = cell_type (only one type)
        
        # For or PFC/V1 (type: paired)
        # id_1 = age
        # id_2 = area
      
        # A tibble: 114 x 7
        # avg.exp pct.exp features.plot id       id_1  id_2  avg.exp.scaled
        # <dbl>   <dbl> <fct>         <fct>    <fct> <chr>          <dbl>
        # 1   0.652    23.7 NR2F1         motor_rg motor rg           -0.428 
        # 2   5.71     82.5 NFIA          motor_rg motor rg            1.74  
        # 3   0.609    20.9 ZBTB18        motor_rg motor rg           -0.497 
        # 4   3.86     71.5 TCF4          motor_rg motor rg            1.35  
        # 5   3.08     54.2 SOX11         motor_rg motor rg            1.13  
        # 6   1.03     34.3 NFIX          motor_rg motor rg            0.0249
        # 7   4.29     74.5 NFIB          motor_rg motor rg            1.46  
        # 8   2.67     54.9 HMGB2         motor_rg motor rg            0.982 
        # 9   0.478    19.6 ZNF462        motor_rg motor rg           -0.739 
        # 10   0.344   16.4 ZIC2         motor_rg motor rg           -1.07  
        #  … with 104 more rows
      
      # type: normal
          if(any(x$id_1 == 'pfc')) { 
            x %<>% mutate(id_1 = factor(id_1, levels = cortical_areas)) 
          }
        
     # # type: paired
     #    if(any(str_detect(x$id_1, '18'))) {
     #      x %<>% mutate(id_1 = as.numeric(id_1)) 
     #    }
         
       #  if(any(x$id_2 == 'pfc')) { 
       #    x %<>% mutate(id_1 = factor(id_1, levels = levels(seurat_object@meta.data$area))) 
       #  }
         
        # if(any(str_detect(x$id_1, '[:digit:]+'))) {
        #   x %<>% mutate(id_1 = as.numeric(id_1)) 
        # }
           
  
  # 3. SCALING ----------------------------------------------------
      message('Scaling. \n ')

x %<>% group_by(features.plot) %>% 
  mutate(avg_exp_quant = scales::squish(avg.exp, quantile(avg.exp, quantiles))) %>% 
  ungroup

      # Scale average expression values by individual.
      if(scale_by == 'gene') {
        
        cat('Scaling each gene across each group (column-wise).\n')
        x %<>% group_by(features.plot) %>%
          mutate(scaled_exp = scale(avg_exp_quant, center = TRUE)) %>%
          ungroup
      } 
      
      # mode == 'paired'
      if(scale_by == 'id_2') { 
       cat('Scaling each gene within each id_2.\n')
       # Scale average expression values by individual.
       x %<>% group_by(id_2, features.plot)  %>%
         mutate(scaled_exp = scale(avg.exp, center = TRUE)) %>%
         ungroup
     } 
      # TODO fix this so the three options can coexist (ifelse?)
       #   else {
       #       cat('Scaling each gene within each id_1.\n')
       #       # Scale average expression values by individual.
       #       x %<>% group_by(id_1, features.plot)  %>%
       #         mutate(scaled_exp = scale(avg.exp, center = TRUE)) %>%
       #         ungroup
       #  }
       message("Scaled\n")
       print(x)

    # END SCALING
     
  # 4. Use x$data to make own ggplot dotplot. ---------------------------------------
  
  
  message('n_genes: \n'); cat(n_genes <- n_distinct(x$features.plot))
          
  message('Max size: \n'); cat(max_size <- 500 / n_genes)

      # PFC / V1 DotPlots across ages. (side-by-side areas) ----------------------------------

      
        if(mode == 'paired-areas') {
          
          # id_1 = area/region; 1 = pfc ; 6 = v1
          # id_2 = age
          
          x %<>% dfilter(id %>% str_detect('pfc|v1') & 
                         !(id_2 %>% str_detect('14|25')))

          x %<>% unite(gene_region, 'features.plot', 'id_1', remove = FALSE)

          x %<>% bind_rows(data.frame(features.plot = x$features.plot, 
                                 gene_region = paste0(x$features.plot, '_NA'), 
                                 scaled_exp = NA)) %>%
             mutate(x_label = case_when(gene_region %>% str_detect('_1') ~ as.character(features.plot),
                                                                    TRUE ~ ''))
          
            ## Prepare for clustering
              ## Group rows by gene and age, then get the diff in expression between pfc and v1.
              ## id_1 = area/region; 1 = pfc ; 6 = v1
          
            x %<>% group_by(features.plot, id_2) %>% 
               mutate(diff = scaled_exp[id_1 == 1] - scaled_exp[id_1 == 6]) %>% ungroup %>%
              
              dfilter(id_1 == 1) %>% 
                select(features.plot, id_2, diff) %>% 
                spread(key = id_2, value = diff)
        }
  
     # General DotPlots: genes in columns, areas in rows.  -----------------------------------
   
         if(mode == 'genes_by_area'){
               
           ## Prepare for clustering
              # Spread
           
            # id_1 = area
             x_spread <- x %>%
              select(features.plot, id_1, scaled_exp) %>% 
               pivot_wider(id_cols = features.plot, names_from = id_1,  values_from = scaled_exp)
         }
         print(x_spread)
         #  > x_spread
           
         #   A tibble: 19 x 7
         #   features.plot    pfc    motor somatosensory parietal temporal      v1
         #   <fct>          <dbl>    <dbl>         <dbl>    <dbl>    <dbl>   <dbl>
         #   1 NR2F1         -1.12  -0.959         -0.432    0.456    1.44    0.612 
         #   2 NFIA          -0.334  1.01           0.814    0.460   -1.71   -0.240 
         #   3 ZBTB18        -0.683 -1.06          -0.346    1.20    -0.413   1.30  
         #   4 TCF4          -1.38  -0.0613        -0.0305   0.0746  -0.334   1.73  
         #   5 SOX11         -0.612 -0.901         -0.856    1.00     1.44   -0.0736
         #   6 NFIX          -1.21  -0.276          0.308    0.136   -0.665   1.71  
         #   7 NFIB          -1.81  -0.320          0.359    0.610    0.111   1.05  
         #   8 HMGB2         -0.308  0.639          0.0911  -1.18    -0.802   1.56  
         #   9 ZNF462        -1.26  -0.738         -0.580    0.399    1.13    1.05  
         #   10 ZIC2          -0.941  0.161         -0.555   -0.581    0.0614  1.85  
         #   11 FOXP1         -1.15  -0.455         -0.431   -0.203    0.505   1.73  
         #   12 FOXO3         -1.29  -0.184          0.898    1.44    -0.567  -0.296 
         #   13 BCL11A        -0.626 -0.947         -0.485   -0.0230   0.230   1.85  
         #   14 SUB1           1.91   0.00897       -0.137   -0.233   -0.556  -0.990 
         #   15 HOPX           1.52   0.650         -0.373   -0.191   -0.160  -1.44  
         #   16 HES1           1.35   0.220          0.236   -0.463   -1.66    0.316 
         #   17 FOS           -0.752  0.118         -0.533   -0.781    0.0706  1.88  
         #   18 DDIT3          0.333 -0.0637        -0.173    1.75    -1.12   -0.729 
         #   19 CEBPD         -0.727 -0.520         -0.345   -0.0106   1.98   -0.382 
           
# 5. CLUSTER GENES ----------------------------------------------------------
  # by similarity
   
   clusterGenes <- function(x, x_spread, method_dist) {    
     
          dist <- x_spread %>% column_to_rownames('features.plot') %>% as.matrix() %>%

                      dist(diag = TRUE, upper = TRUE, method = method_dist)
           
           genes_hclust <- hclust(dist)
           
           # Set order in which genes should be plotted according to genes_hclust$order
           x %<>% mutate(features.plot = features.plot %>% 
                           factor(levels = genes_hclust$labels[genes_hclust$order]),
                          id = as_factor(as.character(id)))
                        
           x %<>% arrange(features.plot, id_1, id_2)
             # mode = paired      
              # %>% mutate(gene_region = as_factor(gene_region))
          return(x)
   }
   
   x_clustered <- clusterGenes(x, x_spread, method_dist = 'euclidean')
  
  # 6. GGPLOT -------------------------------------------------

     message('Making ggplot \n x: \n')
  
  dotplot <- ggplot(x_clustered) + ggtitle(title) +
  
            geom_point(aes(x = features.plot, y = id_1,  # x = gene_region if mode = 'paired'
                 colour = scaled_exp,
                 # colour = avg.exp.scaled,
                 # alpha = scaled_exp,
                 size = pct.exp),
                 shape = 16) + 
          scale_colour_gradientn(colours = viridis(n = 100, option = 'plasma', end = 0.9),
                                 na.value = "white")
    
  
      if(color_scale == 'red_yellow_blue') {
        
        dotplot <- dotplot + 
                      scale_colour_gradientn(colours = rev(heatmap.colors),
                                         # TODO Maybe keep high yellow, low blue gradient.
                                         #values = c(0, 0.3, 1),
                                         # breaks = c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9),
                                         # limits = c(-1.9, 1.9),
                                         oob = scales::squish,
                                         na.value = "white")
      }
  
    #  Refine plot . . . . . . . . . . . . . . . . 
    
    dotplot <- dotplot +
    
      # scale_x_discrete(labels = x$x_label) +
      # scale_alpha(range = c(0.5, 0.9)) +
      scale_size_area(limits = c(2, 90), max_size = max_size) +

      guides(size = guide_legend(direction = "vertical")) +
      
      theme_void() +
      theme(plot.title = element_text(size = 9),
            text = element_text(size = 5, colour = 'grey30'),
            axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1),
            axis.text.y = element_text(size = 5, angle = 0, hjust = 1),
            axis.title = element_blank(),
            legend.position = 'right',
            legend.key.size = unit(0.1, units = 'in'))
        
return(lst(x, x_clustered, dotplot))

} # END makeDotPlot ---------------------------------------------



heatmap.colors <- 
  c("Vermilion" = "#DF4619",
    "Flame"="#E45A17",
    "Marigold"="#F4AD39",
    "Naples Yellow"="#FAD85F",
    "Key Lime"="#F0FD8C",
    "Light Green"="#ADF196",
    "Maximum Blue Green"="#25BFC6",
    "Pacific Blue"="#21AFCB",
    "Blue NCS"="#1A90C1")
