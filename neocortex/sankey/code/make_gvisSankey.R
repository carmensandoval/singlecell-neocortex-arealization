# fxn: make_gvisSankey

#' @param gvis_df_x a data frame with columns: 

# item1                  <chr> "rg_pfc_early", "rg_pfc_early", "rg_pfc_early", "rg_pfc_early",...
# item2                  <chr> "neuron_pfc_early", "neuron_msp_early", "neuron_temporal_early"...
# n                      <dbl> 66, 16, 4, 13, 20, 38, 30, 22, 17, 27, 13, 10, 28, 101, 13, 46,...
# celltype.source        <fct> rg, rg, rg, rg, rg, rg, rg, rg, rg, rg, rg, rg, rg, rg, rg, rg,...
# region.source          <fct> pfc, pfc, pfc, pfc, pfc, pfc, pfc, pfc, pfc, pfc, pfc, pfc, msp...
# stage.source           <fct> early, early, early, early, early, early, early, early, early, ...
# color.region.source    <chr> "#e72718", "#e72718", "#e72718", "#e72718", "#e72718", "#e72718...
# color.rg.region.source <chr> "#FF2B1A", "#FF2B1A", "#FF2B1A", "#FF2B1A", "#FF2B1A", "#FF2B1A...
# celltype.target        <fct> neuron, neuron, neuron, neuron, neuron, neuron, neuron, neuron,...
# region.target          <fct> pfc, msp, temporal, v1, pfc, msp, temporal, v1, pfc, msp, tempo...
# stage.target           <fct> early, early, early, early, mid, mid, mid, mid, late, late, lat...
# color.region.target    <chr> "#e72718", "orange", "#ff61e8", "#1b30f3", "#e72718", "orange",...
# color.rg.region.target <chr> "#FF2B1A", "#FFA500", "#FF61E8", "#1C32FF", "#FF2B1A", "#FFA500...

make_gvisSankey <- function(chart_id, gvis_df_x, out_dir) {
  
  gvis_df_x %<>% 
    mutate(item2 = paste0("to_", item2)) %>% 
    # region_source and region_target should be factors with levels
    # set to pfc, msp, temporal, v1.
      arrange(region_source, region_target)
  # TODO fix color issue [NA is not a valid color string] 2020-10-23
  # remove to_ in item2 for making colors

  groups.order <- gvis_df_x %>% 
    split(f = factor(.$item1, levels = unique(.$item1))) %>% 
      map(~ c(unique(.$item1), .$item2)) %>% 
        unlist(use.names = FALSE) %>% unique %>% 
          tibble(group = . ) %>% 
            mutate(group = str_remove(group, "to_"), 
                   celltype_region = str_extract(group, 
                                                '^[:alnum:]+_[:alnum:]+'))
                    
  groups.order <- left_join(groups.order, 
                            groups.metadata %>% 
                              dplyr::select(-starts_with("color")),
                            by =c("celltype_region" = "group")) %>%
                  
                  left_join(colors.area, 
                            by = (c("region" = "category")))
                          
  colors <- paste("['", 
                 paste(groups.order$color.region, 
                 collapse = "\',\'"),"']", 
                 sep = "")
  
  plot <- gvisSankey(chartid = paste("sankey", chart_id, sep = "_"),
                     data =  gvis_df_x %>% dplyr::select(item1, item2, n), 
                  # Doesn't like having all those extra columns
                from = "item1", 
                to = "item2", 
                weight = "n",
                options = list(width = 1000,
                                height = 800,
                                backgroundColor = "black",
                                tooltip = "{ 
                                      textStyle: { color: 'black' }, 
                                      showColorCode: 'True',
                                      isHTML: 'True'
                                      }",
                                
                                sankey = paste( "{ iterations: 0 ,
                                                         node: { colors : ", colors, ", 
                                                                width : 25,
                                                                nodePadding : 35,
                                                                interactivity: 'True'
                                                                              },
                                                      
                                                  label: { color: '#871b47' } , 
                                                  link: { colorMode: 'gradient' }
                          
                                                 }")
                                        )
            )
  print(plot, tag = "html", file = file.path(out_dir, paste0(chart_id, ".html")))
  return(list(groups.order, colors, gvis_df_x, plot))
}