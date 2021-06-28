cheaVenn <- function(region_use_chea,
                     main.pos = c(0.5, 0.9),
                     region_use_markers,
                     df_marker_tfs = tfs$MAST_wilcox_by_group,
                     df_chea,
                     find_top_chea = FALSE,
                     label_print = FALSE,
                     frac_keep = 0.20,
                     table_use = 'Integrated--meanRank',
                     score_threshold) {
  
  message('\nCombination:\n')
  cat(combination <- paste0('ChEA: ', region_use_chea, ' / ', 'RG markers: ' , region_use_markers))
  
  # Define which TFs to use from the ChEA results table. -------------------------------------------
  
  # a) Find top n TFs -------------------------------------------------------------------
  
  df_chea_split <- df_chea %>% dfilter(table == table_use & 
                                         query_genes %>% str_detect(paste0('neuron_.*_', region_use_chea))) %>% 
    split(.$query_genes)
  
  df_chea_top <- df_chea_split %>% map(~ mutate(.x, score = as.numeric(score)) %>%
                                         top_frac(n = frac_keep, wt = -score)) %>%
    rbindlist %>% distinct(tf, .keep_all = TRUE)
  
  # b) Keep TFs with score < threshold ----------------------------------------------------------
  
  # df_chea <- df_chea %>% mutate(score = as.numeric(score)) %>%
  #                         dfilter(table == table_use &
  #                           cell_type == 'neuron' & region == region_use &
  #                             score <= score_threshold)
  
  # # Count how many TFs above score threshold for each stage
  # tfs_pred_per_stage <- df_chea %>% group_by(query_genes) %>% 
  #                         summarise(n_distinct(tf))
  # 
  tfs_chea <- df_chea_top %>% count(tf) %>% .$tf
  
  warning('\nTFs ChEA:\n')
  cat_n(tfs_chea)
  
  # Get RG marker TFs for this area/region ---------------------------------------------------------
  
  df_marker_tfs <- df_marker_tfs %>% ungroup %>% 
    dfilter(cell_type == 'rg' & region == region_use_markers)
  
  tfs_markers <- df_marker_tfs %>% unlist(x = .$genes, use.names = FALSE) %>% sort %>% unique
  
  # Calculate overlap --------------------------------------------------------
  
  compare_list <- lst(tfs_markers, tfs_chea)
  
  # VennDiagram
  
  # venn_result <- VennDiagram::calculate.overlap(x = lst(tfs_rg_markers, tfs_predicted))
  venn_result <- get.venn.partitions(x = compare_list)
  tf_markers_only <- venn_result %>% 
    dfilter(tfs_markers == TRUE & tfs_chea == FALSE) %>% pull(`..values..`)
  
  # venndetail
  
  venn_detail <- VennDetail::venndetail(x = compare_list)
  venn_detail_result <- result(venn_detail)
  
  input <- venn_detail@input %>% enframe
  result <- venn_detail@result %>% group_by(Subset) %>% summarise(genes = list(as.character(Detail)))
  
  venn_detail_df <- rbindlist(lst(input, result), idcol = 'type') %>%
    mutate(n_genes = map_int(value, length)) # value %>% map(length) %>% unlist
  
  print(venn_detail_df %>% dfilter(name == 'Shared'))

  n_shared <- venn_detail_df %>% dfilter(name == 'Shared') %>% pull(n_genes)
  if(nrow( venn_detail_df %>% dfilter(name == 'Shared')) == 0 ){ n_shared <- 0.01 }
  message('\nn_shared: ')
  print(n_shared)
  cat('\nn_shared: ', n_shared)

  venn_detail_df %<>% mutate(pct_shared = round(n_shared / n_genes, 2))
  
  #  Make Venn diagram 
  # a) VennDiagram::venn.diagram
  
  # p <- venn.diagram(x = compare_list, force.unique = TRUE, 
  #                   filename = NULL, main = paste('TFs ', region_use),
  #                   ) %>% grobTree(name = region_use)
  #                           # How to make objects from VennDiagram compatible with 
  #                           # cowplot plot_grid?
  #                           #https://stackoverflow.com/a/51006231/4463919
  # 
  
  # b) Add item labels (adapted from RAM::venn.group)
  p_venn <- vennGroup(title = combination,
                      main.pos = main.pos,
                      compare_list, 
                      label = label_print, 
                      fill = c("orange", "blue"),
                      cat.pos = c(-1, 1),
                      cat.dist = c(0.02, 0.02),
                      cat.cex = 1,
                      lab.cex = 1, width = 20, height = 20) %>%
    
    grobTree(name = combination)
  
  # ggsave(plot = p, file.path(dir_dropbox, paste0('venn_tf_', region_use, '.pdf')), 
  #       width = 5, height = 5)
  
  # -------------------------------------------
  # TF score distribution of selected ChEA TFs 
  
  scoreDistribution <- function(df) {
    
    p <- df %>% 
      ggplot(aes(x = score)) +
      geom_density(aes(fill = query_genes, colour = query_genes), alpha = 0.3) +
      geom_density(alpha = 0.2, colour = 'grey20', linetype = 'dashed', show.legend = FALSE) +
      theme_minimal()
    
    #  ggsave(plot = p, 
    #     filename = file.path(dir_dropbox, paste0('tf_score_distribution_', region_use, '.pdf')), 
    #     width = 7, height = 5)  
  }
  
  p_scores <- df_chea_top %>% scoreDistribution
  
  return(lst(venn_detail_df, compare_list, venn_result, venn_detail, venn_detail_result, 
             plots = lst(p_venn, p_scores)))
  
}

# . . . . . . . . . . .



vennGroup <- function (vectors, 
                       title,
                       main.pos,
                       main.cex = 1,
                       cat.cex = 1.5, cex = 1, 
                       cat.pos = NULL, cat.dist = NULL, 
                       label = TRUE, lab.cex = 1, lab.col = "black", 
                       fill = NULL, 
          file = NULL, ext = NULL, width = 8, height = 8) 
{
  save <- !is.null(file)
  if (save) {
    .get.dev(file, ext, height = height, width = width)
  }
  if (!requireNamespace("VennDiagram")) {
    stop("package 'VennDiagram' is required for this function")
  }
  if (!requireNamespace("RColorBrewer")) {
    stop("package 'RColorBrewer' is required for this function")
  }
  if (!requireNamespace("grid")) {
    stop("package 'grid' is required to use this function")
  }
  len <- length(vectors)
  if (is.null(fill)) {
    if (len == 2) {
      fill = c("lightpink", "lightblue")
    }
    else {
      fill = RColorBrewer::brewer.pal(len, "Pastel1")
    }
  }
  else {
    if (length(fill) == len) {
      fill = fill
    }
    else if (length(fill) > len) {
      warning(paste("more colors being provided than required, will ignore ", 
                    length(fill) - len, " colors", sep = ""))
      fill = fill[1:len]
    }
    else {
      warning("not enough colors being provided, will use default")
      if (len == 2) {
        fill = c("lightpink", "lightblue")
      }
      else {
        fill = RColorBrewer::brewer.pal(len, "Pastel1")
      }
    }
  }
  if (len > 2 && label) {
    warning("currently only support 2 groups to have actual item labels; will only use numbers")
  }
  else if (len > 5 || len < 2) {
    stop("please provide 2 to 5 vectors")
  }
  
  alpha = rep(0.5, len)
  
  if (!is.null(cat.pos) && !is.null(cat.dist)) {
    
    v <- VennDiagram::venn.diagram(x = vectors, 
                                   force.unique = TRUE, 
                                   ext.text = FALSE,
                                   fill = fill, 
                                   alpha = alpha, 
                                   
                                   main = title,
                                   main.pos = main.pos,
                                   main.fontfamily = 'serif',
                                   main.cex = main.cex,
                                   fontfamily = rep('sans', 3),
                                   cat.fontfamily = rep('sans', 2),

                                   cat.dist = cat.dist, cat.pos = cat.pos, 
                                   cat.cex = cat.cex, cex = cex, 
                                   filename = NULL)
  }
  
  else if (!is.null(cat.pos) && is.null(cat.dist)) {
    v <- VennDiagram::venn.diagram(vectors, fill = fill, 
                                   alpha = alpha, cat.pos = cat.pos, cat.fontface = "bold", 
                                   cat.cex = cat.cex, cex = cex, filename = NULL)
  }
  else if (is.null(cat.pos) && !is.null(cat.dist)) {
    v <- VennDiagram::venn.diagram(main = title_use, vectors, fill = fill, 
                                   alpha = alpha, cat.fontface = "bold", cat.dist = cat.dist, 
                                   cat.cex = cat.cex, cex = cex, filename = NULL)
  }
  else {
    v <- VennDiagram::venn.diagram(vectors, fill = fill, 
                                   alpha = alpha, cat.fontface = "bold", cat.cex = cat.cex, 
                                   cex = cex, filename = NULL)
  }
  if (len > 2 && len <= 5) {
    grid::grid.newpage()
    grid::grid.draw(v)
  }
  if (len == 2) {
    if (!label) {
      grid::grid.newpage()
      grid::grid.draw(v)
    }
    else {
      name <- lapply(v, names)
      
      v.labels <- lapply(v, function(i) i$label)
      v.lab <- vector()
      for (i in 1:length(v.labels)) {
        if (length(v.labels[[i]] %in% names(vectors)) != 
            0 && isTRUE(v.labels[[i]] %in% names(vectors))) {
          v.lab <- c(v.lab, v.labels[[i]])
        }
      }
      
      v1 <- vectors[[v.lab[1]]]
      v2 <- vectors[[v.lab[2]]]
      
      v[[5]]$label <- paste(c(v[[5]]$label, 
                              setdiff(v1, v2) %>% head(30)), 
                              collapse = "\n"
                              )
      
      v[[5]]$gp$cex <- lab.cex
      v[[5]]$gp$col <- lab.col
      
      v[[6]]$label <- paste(c(v[[6]]$label, 
                              setdiff(v2, v1)), 
                              collapse = "\n")
      
      v[[6]]$gp$cex <- lab.cex
      v[[6]]$gp$col <- lab.col
      
      # Intersecting genes
      v[[8]]$label <- paste(c(v[[8]]$label, 
                              intersect(v1, v2)), 
                              collapse = "\n")
      v[[8]]$gp$cex <- lab.cex
      v[[8]]$gp$col <- lab.col
      
      grid::grid.newpage()
      grid::grid.draw(v)
    }
  }
  
  if (save) {
    dev.off()
  }
  
  return(v)
  
}

