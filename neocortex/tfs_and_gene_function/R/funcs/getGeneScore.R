# Calculate enrichment ratio, filter genes w. adj p-value < 0.05, sort table.
# Positive values indicate that the feature is more highly expressed in the first group.

# date: 2020-11-27
# last modified: 2020-11-27

getGeneScore <- function(markers_df) {

  # markers %<>% rownames_to_column("gene") %>%
  markers_df %<>% mutate(enrichment.ratio = pct.1 / (pct.2 + .000001),
                        diff.pct = pct.1 - pct.2,
                        gene.score = avg_logFC * enrichment.ratio)

  markers_df %<>% 
    # filter(p_val_adj <= 0.1) %>% 
     # dplyr::filter(pct.1 >= 0.5 | pct.2 >= 0.5) %>%
        dplyr::select(gene, cluster, cell_type, individual, individual,
                      gene.score, p_val_adj, diff.pct, enrichment.ratio,
                      pct.1, pct.2, avg_logFC) %>%
          mutate(across(is.numeric, .f = ~ round(.x, 5)))

          # print(markers_df)

          markers_df %>% group_by(cluster, cell_type, individual) %>%
            arrange(desc(gene.score), .by_group = TRUE)

# write_tsv(markers, path =  paste0("../out/DEgenes_", cluster_a, "_vs_", cluster_b, "scale.data.tsv"))
}
