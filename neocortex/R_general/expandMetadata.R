expandMetadata <- function(data_frame) {
            
  data_frame %>% set_names(names(.) %>% tolower) %>%
  # data_frame can be ex. @meta_data or markers df.
    mutate(cell_type = tolower(cell_type) %>% factor(levels = c('rg', 'ipc', 'neuron', 
                                                                'opc', 'cr', 'astrocyte', 'dividing')),
           #structure = tolower(structure),
           individual = as.factor(tolower(individual)), 
           # ^ Use base R `as.factor` so it orders individuals alphabetically/numerically
           # (instead of in the order in which they appear with as_factor)
           stage = fct_collapse(individual,
                                early = c('gw14', 'gw16', 'gw17'),
                                 mid = c('gw18_2', 'gw18', 'gw19_2', 
                                         'gw19', 'gw20', 'gw20_31', 'gw20_34'), #,'gw20_31and34'),
                                 late = c('gw22', 'gw22t', #'gw22both',
                                          'gw25')),
           area = factor(tolower(area), 
                         levels = c('pfc', 'motor', 'somatosensory', 
                                    'parietal', 'temporal', 'v1')),
           region = fct_collapse(area, 
                                 msp = c('motor', 'somatosensory', 'parietal'))
           ) %>%
           unite('stage_region', stage, region, sep = '_', remove = FALSE) %>%
           unite('stage_area', stage, area, sep = '_', remove = FALSE) %>%
           unite('celltype_stage', cell_type, stage, sep = '_', remove = FALSE) %>%
           unite('celltype_region', cell_type, region, sep = '_', remove = FALSE)

}
