set_here("~/cse-phd/second-trimester/neocortex")
p_load(here)
here()

conflict_prefer("combine", "dplyr")
source(here("tfs_and_gene_function/R/funcs/makeDotPlot.R"))

markers_caldwell <- read_tsv(here("markers/caldwell_neuron_table.tsv"), comment = "//")
markers_allen <- read_tsv(here("markers/brainspan_RNAseq.tsv"), comment = "//")
markers <- c(markers_caldwell[,1], markers_allen[,1]) %>% combine %>% unique

p_load(Seurat)
neocortex <- read_rds(here("data/210322_neocortex_full_v3.2.3.rds"))
cortical_areas <- c("pfc", "motor", "somatosensory", "parietal", "temporal", "v1")

Idents(neocortex) <- neocortex@meta.data$cell_type


dotplot <- makeDotPlot(seurat_object = neocortex, 
                       mode = 'genes_by_area',
                       idents_use = NULL,
                       # group = 'rg_early',
                       # area = 'pfc; v1', 
                       genes_list = markers,
                       group_by = 'cortical_area',
                       split_by = NULL,
                       scale = FALSE,
                       title = 'All neocortex cells',
                       #color_scale = 'red_yellow_blue', 
                       color_scale = 'viridis',
                       scale_by = 'gene',
                       max_size = 4,
                       quantiles = c(0, 0.9))

fs::dir_create(out_dir <- here("dotplots_area_validation/out"))

ggsave(plot = dotplot$dotplot, 
       filename = file.path(out_dir, 
                            paste0("markers_viridis_", Sys.time(), ".pdf")),
       width = length(markers)/5.5, height = 5/3, 
       units = "in", device = "pdf")
