p_load(rlang)
p_load(cowplot)

# featurePlotNice

featurePlotNice <- function(seurat.object = seurat.object,
                            version = v3,
                            integrated = FALSE,
                            exp.df  = data,
                            gene = X,
                            cap = cap,
                            pt.size = 0.1, pt.size.grey = 0.01,
                            alpha.group = 0.7, alpha.grey = 0.5,
                            begin = 0,
                            end = 1,
                            color.scale = "viridis",
                            legend.pos = "none"
                              ) {
                                
                   gene <- enquo(gene)
                   exp.df <- enquo(exp.df)
                   version <- enquo(version)
                   
                   print(quo_text(version))
                   
                   if(quo_text(version) == "v2") {
                    
                   umap.embeddings <- seurat.object@dr$umap@cell.embeddings %>%
                     data.frame %>% rownames_to_column("cell.name")
                   
                   exp.df <- seurat.object %>% slot(quo_text(exp.df))
                   
                   }
                 
                   if(quo_text(version) == "v3") {
                    
                         umap.embeddings <- seurat.object@reductions$umap@cell.embeddings %>%
                           data.frame %>% rownames_to_column("cell.name")
                   
                         if(integrated == FALSE) {
                           exp.df <- seurat.object@assays$RNA %>% slot(quo_text(exp.df))
                         }
                   
                         if(integrated == TRUE) {
                           exp.df <- seurat.object@assays$integrated %>% slot(quo_text(exp.df))
                                                                              
                         }
                   
                   }
                   
                   
                   
     exp.df <- exp.df[as_name(gene), ] %>% as.data.frame %>%
                  set_colnames("exp") %>%
                    rownames_to_column("cell.name")
  

     head(umap.embeddings) %>% print
    
     meta <- seurat.object@meta.data %>% rownames_to_column ("cell.name")

     head(meta) %>% print
    
     umap.metadata <- left_join(meta, umap.embeddings)
     umap.metadata %>% head %>% print
    
 #   gene.text <- quo_text(gene)
    
 #   df <- exp.df %>% dplyr::select(cell.name, !!gene) %>%
    
 df <- left_join(umap.metadata, exp.df, by = "cell.name") %>% arrange(exp)
 
      vln_plot <- df %>% ggplot() + 
                  geom_violin(aes(x = Area, y = exp, fill = Area))

 print(head(df))

 df <- df %>% mutate(exp.cap = case_when(exp > cap ~ cap,
                                              TRUE ~ exp)
              )
 
 head(df)

 p <- ggplot(df) +

   # Background cells
   geom_point(data = dplyr::filter(df, exp == 0),
              aes(UMAP_1, UMAP_2),
              colour = "grey80", size = pt.size.grey,
              alpha = alpha.grey, shape = 16) +

   # Cells with non-zero expression of the gene.
   geom_point(data = dplyr::filter(df, exp > 0),
              aes(UMAP_1, UMAP_2,
                  colour = exp.cap),
              size = pt.size, 
              alpha = alpha.group,
              shape = 16) +

 # geom_vline(xintercept = quan)

 ggtitle(as_name(gene)) +
   
   scale_color_viridis_c(option = "magma",
                         begin = begin,
                         end = end, 
                         name = "expression",
                         guide = guide_colourbar(ticks = TRUE)
                         breaks = c(0.1, cap),
                         labels =c("min", "max")) +
                         # breaks = seq(0, 6, by = 1)) +fffdffdf
   
   # scale_alpha(range = c(0.3, 0.8))
   # guides(alpha = FALSE) +
   
   theme(legend.position = legend.pos,
         legend.key.width = unit(0.25,"cm"),
         # legend.text = element_blank(),
         legend.title = element_text(size = 6),
         panel.grid = element_blank(),
         axis.title = element_blank(),
         axis.text = element_blank(),
         axis.ticks = element_blank(),
         panel.background = element_blank(),
         axis.line = element_blank(),
         panel.border = element_blank(),
         aspect.ratio = 1
   )

 return(list(vln_plot, p))
}

# Example

# featurePlotNice(seurat.object = wb.small,
#               exp.df = data,
#               gene = SOX2,
#               legend.pos = "right",
#               cap = 3.8,
#               pt.size = 0.25,
#               pt.size.grey = 0.2,
#               alpha.group = 0.7,
#               alpha.grey = 0.3,
#               begin = 0,
#               end = 0.95)


# fx: histogramExpression 
histogramExpression <-  function(gene = SOX2, seurat.object = seurat.object,
                                integrated = FALSE, version = v3, exp.df = data) {

          gene <- enquo(gene)
          version <- enquo(version)
          exp.df <- enquo(exp.df)
          # gene <- quo(SOX2)

        if(quo_text(version) == "v2") {
        
          gene.exp <- seurat.object@data[quo_text(gene), ] %>%
            as.data.frame() %>% setNames(quo_text(gene))

            }

        if(quo_text(version) == "v3") {
        
              if(integrated == FALSE) {
                gene.exp <- seurat.object@assays$RNA %>% slot(quo_text(exp.df)) %>%
                  # Gene names are rows in @data.
                  .[quo_text(gene), ] %>%
                  as.data.frame() %>% setNames(quo_text(gene))
              }

              if(integrated == TRUE) {
                gene.exp <- seurat.object@assays$integrated %>% slot(quo_text(exp.df)) %>%
                  .[quo_text(gene), ] %>%
                  as.data.frame() %>% setNames(quo_text(gene))
              }

        }

          print(gene.exp %>% str)

          max.gene <- gene.exp %>% pull(!!gene) %>% max
          print(max.gene)

          p <-  ggplot(gene.exp %>% dplyr::filter(!!gene > 0.1 )) +
                        geom_histogram(aes(x = !!gene), fill = "grey50", bins = 100, alpha = 0.5) +
                        scale_x_continuous(breaks = seq(0, max.gene, by = 0.5)) +
                        theme_minimal()  +
                        theme(title = element_blank(),
                              axis.title = element_blank(),
                              axis.text = element_text(size = 5, angle = 45)
                              )

          return(p)
        }

       # Example

       # histogramExpression( seurat.object = wb.small,
       #                      gene = SOX2, 
       #                     exp.df = data)


# fx: FeatPlotHist

featPlotHist <- function(seurat.object,
                         exp.df = data,
                         gene = SOX2,
                         alpha.group = 0.7,
                         alpha.grey = 0.3,
                         pt.size.grey = 0.2,
                         pt.size = 0.3,
                         cap = cap,
                         legend.pos = "right",
                         version = v3,
                         integrated = FALSE,
                         begin = 0,
                         end = 0.95) {

  gene <- enquo(gene)
  exp.df <- enquo(exp.df)
  version <- enquo(version)

  p1 <- featurePlotNice(seurat.object = seurat.object,
                        integrated = integrated,
                        exp.df = !!exp.df,
                        gene = !!gene,
                        legend.pos = legend.pos,
                        cap = cap,
                        pt.size = pt.size,
                        pt.size.grey = pt.size.grey,
                        alpha.group = alpha.group,
                        alpha.grey = alpha.grey,
                        begin = begin,
                        end = end,
                        version = !!version)
  

  p2 <- histogramExpression(!!gene, seurat.object = seurat.object,
                            integrated = integrated, version = !!version, exp.df = !!exp.df)

  # print(grid.arrange(p1, p2, ncol = 2, nrow = 1))
  # print(p1 + annotation_custom(ggplotGrob(p2), xmin = 1, xmax = 3,
  #                                       ymin = -0.3, ymax = 0.6)
  # )

  print(p1[[1]])
  

  print(ggdraw() +
          draw_plot(p1[[2]]) +
          draw_plot(p2, x = 0 , y = 0.8, 
                    width = 0.2, height = 0.2)
  )


  
}

## Example

# fPlots <- list(featPlotHist(wb.small, cap = 3.75)
#               featPlotHist(gene = SATB2, wb.small, cap = 3.5, begin = 0, end = 0.7, alpha.grey = 0.3, alpha.group = 0.7)
#               featPlotHist(seurat.object = wb.small, gene = VIM, cap = 5.5)
#               featPlotHist(gene = EOMES, seurat.object = wb.small, cap = 3.25)
#)

# plot_grid(plotlist = fPlots)

# genes <- c(SOX2, 
#           SATB2, BCLB11, NEUN, TBR1,
#           PPP1R17, TBR2)

makePlot <- function(gene, cap = Inf) {
  
                gene = enquo(gene)
                
                featPlotHist(gene = !!gene,
                             cap = cap,
                             seurat.object, 
                             begin = 0, end = 0.95,
                             alpha.grey = 0.3, alpha.group = 0.7)
  
}

# makePlot(gene=SOX2, cap = 3.9)
