library(Seurat)
library(tidyverse)

data <- schard::h5ad2seurat(paste(getwd(),'/data/pbmc_ct0.60_noribo.h5ad', sep=""))

ss <- subset(data, subset = cell_type_0.60 %in% c("CD16+ monocytes"))

ss <- preprocess(ss)

loads <- Loadings(ss, reduction="pca")

loads_df <- as.data.frame(loads)

loads_df$gene <- rownames(loads_df)

top_genes_long <- loads_df %>%
  pivot_longer(cols = starts_with("PC_"), names_to = "PC", values_to = "loading") %>%
  filter(PC %in% c("PC_1", "PC_2", "PC_3")) %>%
  group_by(PC) %>%
  slice_max(order_by = abs(loading), n = 10) %>%
  ungroup()

# Plot
ggplot(top_genes_long, aes(x = reorder(gene, loading), y = loading, fill = PC)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~ PC, scales = "free_y") +
  labs(title = "CD16+ Monocytes: Top 10 PCA Loadings for PCs 1–3", x = "Gene", y = "Loading") +
  scale_fill_manual(values = c("PC_1" = "#1f77b4", "PC_2" = "#ff7f0e", "PC_3" = "#9467bd")) + 
  scale_y_continuous(breaks = seq(-0.1, 0.1, by=0.1)) +
  theme_minimal()
