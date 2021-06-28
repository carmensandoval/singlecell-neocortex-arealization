Permissions Size User      Date Modified Name
.rw-r--r--@ 8.1M csandoval 13 May 17:26  PFC_gw16.csv
.rw-r--r--@ 4.3M csandoval 13 May 17:26  SS_gw16.csv
.rw-r--r--@ 2.7M csandoval 13 May 17:26  Temp_gw16.csv
.rw-r--r--@  10M csandoval 13 May 17:26  V1_gw16.csv

cell_by_feature_matrices/gw20:
  Permissions Size User      Date Modified Name
.rw-r--r--@ 1.1M csandoval 13 May 17:29  PFC_gw20.csv
.rwxr--r--@  18M csandoval 13 May 17:29  Somatosensory_gw20.csv
.rwxr--r--@  13M csandoval 13 May 17:30  Temporal_gw20.csv
.rwxr--r--@  10M csandoval 13 May 17:30  V1_gw20.csvsoma



here::i_am("./vitessce.R")

pfc <- read_csv("cell_by_feature_matrices/PFC.csv")  
pfc %>% glimpse

cleanSpatialDF <- function(spatial_df) {
  
spatial_df %<>% select(cell_no = `Cell No.`, 
               nucleus_x = `Nucleus location (X)`, nucleus_y = `Nucleus location (Y)`, 
               nucleus_perimeter = `Nucleus perimeter`,
               everything())
}

plotSpatial <- function(spatial_df, orientation = "vertical") {
                        # max_value_SOX2 = Inf, 
                        # max_value_SATB2 = Inf) {
  
  # spatial_df %<>%
  # mutate_at(.vars = c("nucleus_x", "nucleus_y"), .funs = ~ .x * -1) %>% 
  # filter(!SATB2 == 0)
  
  p1 <- ggplot(spatial_df) + 
    geom_point(aes(x = nucleus_x, y = nucleus_y, colour = SOX2), size = 0.5) +
    scale_color_viridis_c(option = "C",
                          limits = c(0, 10))
                          
  p2 <- ggplot(spatial_df) + 
    geom_point(aes(x = nucleus_x, y = nucleus_y, colour = SATB2), size = 0.5) +
    scale_color_viridis_c(option = "C",
                          limits = c(0, 50))
  
  if(orientation == "vertical") {
    return(p1 / p2) } else {
      return(p1 + p2)
    }
  
  }

processSpatial <- function(spatial_df, orientation = "vertical") {
  
  spatial_df %>%
    cleanSpatialDF() %>% 
    plotSpatial(orientation = orientation)
}

temporal <- read_csv("cell_by_feature_matrices/Temporal.csv")  
temporal <- cleanSpatialDF(temporal)
plotSpatial(temporal)

v1 <- read_csv("cell_by_feature_matrices/V1.csv")
v1 <- cleanSpatialDF(v1)
plotSpatial(v1)

somatosensory <- read_csv("cell_by_feature_matrices/Somatosensory.csv") %>% 
  cleanSpatialDF() %>% 
  somatosensory %>%  select(cell_no, nucleus_x = nucleus_y, nucleus_y = nucleus_x, everything()) %>% 
  mutate_at(.vars = c("nucleus_x", "nucleus_y"), .funs = ~ .x * -1) 

plotSpatial(somatosensory)

somatosensory %>% 
  select(cell_no, nucleus_x = nucleus_y, nucleus_y = nucleus_x, 
         everything()) %>%
write_csv("cell_by_feature_matrices/Somatosensory.csv")

# ---

# GW16
conflict_prefer("set_names", "magrittr")
cell_by_feature_files <- list.files("~/vitessce-human-neocortex-example/cell_by_feature_matrices/gw16",
                                    full.names = TRUE)

cell_by_feature <- map(cell_by_feature_files, read_csv) %>% 
  set_names(value = str_replace(string = cell_by_feature_files, p
                                attern = ".*gw16/([:alnum:]+)_GW16.csv",
                                "\\1"))

plots <- cell_by_feature %>% 
  imap( ~ plotSpatial(.x) %>% 
          ggsave(filename = paste0(.y, ".pdf")))


cell_by_feature[1:3] %<>% map(~ mutate(.x,
tmp_x = nucleus_x,
nucleus_x = nucleus_y,
nucleus_y = tmp_x) %>%
  select( - tmp_x)
)

cell_by_feature$PFC %<>% mutate(nucleus_y = - nucleus_y)

plots <- cell_by_feature %>% 
  imap( ~ plotSpatial(.x) %>% 
          ggsave(filename = paste0(.y, ".pdf")))


# Save corrected dataframes
cell_by_feature %>% 
  imap( ~ write_csv(.x, 
                    file = paste0("~/vitessce-human-neocortex-example/cell_by_feature_matrices/gw16/",
                                                        .y , "_gw16.csv")))

