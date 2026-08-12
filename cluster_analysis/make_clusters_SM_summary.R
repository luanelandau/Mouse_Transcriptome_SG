library(dplyr)
library(readr)

clusters <- read_csv("clusters_SM.csv", show_col_types = FALSE)

summary_table <- clusters |>
  group_by(Cluster) |>
  summarise(
    Number_of_genes = n(),
    Genes = paste(Geneid, collapse = ", "),
    .groups = "drop"
  ) |>
  rename(Cluster_name = Cluster)

write_csv(summary_table, "clusters_SM_summary.csv")
