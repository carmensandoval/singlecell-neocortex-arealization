# Add gw19_2 neurons to ncx_full
# 2020-12-12

# {r 2020-12-12_1733}

load('../../../data/gw19_2.RData')

gw19_2 <- UpdateSeuratObject(gw19_2)

annot_gw19_2 <- read_tsv(file.path(data.dir, 'tbls/gw19_2_clusteridentity_wholebrain_annotated.txt')) %>%
  mutate(area = cell_name %>% str_replace('gw19_2_([:alnum:]+)_[:alpha:]+', '\\1') %>%
           str_replace('ss|SS', 'somatosensory'),
         across(.cols = c(cell_type, area), tolower),
         structure = str_replace(area, areas.rgx, 'neocortex'),
         cell_type = str_replace(cell_type, 'excitatory neuron', 'neuron'))

annot_gw19_2_ncx <- annot_gw19_2 %>% dplyr::filter(structure == 'neocortex')

gw19_2@meta.data %<>% rownames_to_column('cell_name') %>% 
  set_names(tolower(names(.))) %>% 
  mutate(individual = 'gw19_2') %>% 
  left_join(annot_gw19_2_ncx %>% select(cell_name, cell_type)) %>%
  select(cell_name, structure, area, cell_type, everything()) %>%
  # select(-structure.1, -area.1) %>%
  mutate(across(.cols = c(structure, area, cell_type), tolower)) %>%
  set_rownames(.$cell_name)

gw19_2_ncx_neuron <- subset(gw19_2, subset = (structure == 'neocortex' & cell_type == 'neuron'))

ncx_full <- merge(x = ncx_full, gw19_2_ncx_neuron)
# Some cell names are duplicated across objects provided. Renaming to enforce unique cell names. WTF.

ncx_full@meta.data %<>% select(ngene:region) %>% expandMetadata()
sms('Updated ncx_full')

# 2020-12-12
# save_rds('../../../data/exn_lineage/neocortex_exn_seuratobj.rds')

# Check numbers per cell type in gw19_2 neocortex cells.
# ncx_full@meta.data %>% 
#  dfilter(individual == 'gw19_2') %>% count(cell_type)

# 1 neocortex        cr    63
# 2 neocortex  dividing  2499
# 3 neocortex       ipc  2763
# 4 neocortex    neuron 16380
# 5 neocortex       opc   807
# 6 neocortex        rg  2496