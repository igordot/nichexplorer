# simplify and clean up the original Seurat object for nichExplorer Shiny app

library(tidyverse)
library(Seurat)

# load the original Seurat object (v2)
so = readRDS("/path/seurat_obj.rds")

# expression matrix
exp_mat = so@data
saveRDS(exp_mat, "data.exp.rds")

# meta data table
meta_tbl = so@meta.data %>%
  rownames_to_column("cell") %>% as_tibble() %>%
  mutate(
    # rename and reorder the clusters to match the figures
    cluster = str_replace(ident11named, ".*-", ""),
    cluster = factor(cluster, levels = c("V1", "V2", "P1", "P2", "P3", "P4", "P5", "O1", "O2", "O3", "C")),
    # longer cluster names
    cluster_long =
      fct_recode(
        cluster,
        "V1 (Ly6a)" = "V1",
        "V2 (Stab2)" = "V2",
        "P1 (Mgp)" = "P1",
        "P2 (Lpl)" = "P2",
        "P3 (Wif1)" = "P3",
        "P4 (Spp1 Ibsp)" = "P4",
        "P5 (Gas6 Hp)" = "P5",
        "O1 (Col16a1 Tnn)" = "O1",
        "O2 (Fbn1 Igf1)" = "O2",
        "O3 (Bglap Car3)" = "O3",
        "C (Cycling)" = "C"
      ),
    # rename and reorder the treatment to match the figures
    treatment = str_replace(ident2, "^FU", "5FU"),
    treatment = factor(treatment, levels = c("CTRL", "5FU"))
  ) %>%
  select(cell, sample_name = orig.ident, cluster, cluster_long, treatment)
tsne_tbl = so@dr$tsne@cell.embeddings %>%
  round(3) %>% as.data.frame() %>% rownames_to_column("cell") %>% as_tibble()
meta_tbl = full_join(meta_tbl, tsne_tbl, by = "cell") %>% arrange(cell)
saveRDS(meta_tbl, "data.meta.rds")


