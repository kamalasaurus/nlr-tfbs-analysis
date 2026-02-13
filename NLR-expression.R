
raman_stuff <- read.csv2('/Users/kamal/Documents/NYU/raman-expression-classification-analysis/hits_summarized.tsv', header = TRUE, 
          sep = "\t", 
          dec = ".", 
          stringsAsFactors = FALSE, 
          na.strings = c("", "NA", "NaN")
)

# create a df of 2 columns from the raman_stuff df with geneID and Incomplete, then convert Incomplete to logical where '-' is FALSE and anything else is TRUE
raman_incomplete <- data.frame(gene = raman_stuff$geneID, incomplete = raman_stuff$Incomplete != '-')


library(Seurat)
library(ggplot2)

get_data <- function(seurat_object) {
  return(seurat_object@assays$RNA@data)
}

normalize_data <- function(seurat_object) {
  return(Seurat::NormalizeData(
    object             = seurat_object,
    assay              = "RNA",
    normalization.method = "LogNormalize",
    scale.factor       = 10000
  ))
}

aggregate_data <- function(seurat_object, group_by) {
  return(Seurat::AggregateExpression(
    object   = seurat_object,
    assays   = "RNA",
    slot     = "data",
    group.by = group_by
  ))
}

rosette_21d <- readRDS("/Volumes/PhData2/SalkExpressionData/downloads/GSE226097_rosette_21d_230221.rds")
rosette_30d <- readRDS("/Volumes/PhData2/SalkExpressionData/downloads/GSE226097_rosette_30d_230221.rds")
seedling_12d <- readRDS("/Volumes/PhData2/SalkExpressionData/downloads/GSE226097_seedling_12d_230221.rds")
seedling_6d <- readRDS("/Volumes/PhData2/SalkExpressionData/downloads/GSE226097_seedling_6d_230221.rds")
seedling_3d <- readRDS("/Volumes/PhData2/SalkExpressionData/downloads/GSE226097_seedling_3d_230221.rds")

conserved_gene_ids <- readLines('/Volumes/PhData2/SalkExpressionData/unique_ortho_27_genes.txt')



object_names <- c("rosette_21d", "rosette_30d", "seedling_12d", "seedling_6d", "seedling_3d")

normalized <- c(rosette_21d, rosette_30d, seedling_12d, seedling_6d, seedling_3d) |>
  lapply(normalize_data)

cell_type_aggregates <- normalized |>
  lapply(\(obj) aggregate_data(obj, group_by = "CellType")) |>
  lapply(\(obj) obj$RNA)

seurat_clusters_aggregates <- normalized |>
  lapply(\(obj) aggregate_data(obj, group_by = "seurat_clusters")) |>
  lapply(\(obj) obj$RNA)

names(normalized) <- object_names
names(cell_type_aggregates) <- object_names
names(seurat_clusters_aggregates) <- object_names

mesophyll_markers <- c('AT1G29910', 'AT3G01500', 'AT5G16560', 'AT4G12970', 'AT2G42840')
names(mesophyll_markers) <- c('CAB3', 'CA1', 'KAN1', 'EPFL9', 'PDF1')

sorted_mesophyll_markers <- seurat_clusters_aggregates$rosette_21d[mesophyll_markers,] |>
  as.matrix() |>
  apply(1, (\(x) names(sort(x)))) |>
  t()

seurat_clusters_aggregates$rosette_21d['AT1G29910',] |> sort() |> names()


sorted_all <- lapply(
  seurat_clusters_aggregates,              # loop over each list element
  function(mat) {                          # mat is a gene × cluster matrix
    mat[mesophyll_markers, ] |>            # 1. Subset rows by marker genes
      as.matrix() |>                       # 2. Ensure it’s a plain matrix
      apply(1, function(x)                # 3. For each row...
        names(sort(x, decreasing = T))                    #    …sort values & pull names
      ) |>
      t()                                  # 4. Transpose back to genes × rank
  }
)

sorted_values <- lapply(
  seurat_clusters_aggregates,              # loop over each list element
  function(mat) {                          # mat is a gene × cluster matrix
    mat[mesophyll_markers, ] |>            # 1. Subset rows by marker genes
      as.matrix() |>                       # 2. Ensure it’s a plain matrix
      apply(1, function(x) sort(x, decreasing = T)) |>
      t()                                  # 4. Transpose back to genes × rank
  }
)

# different number of clusters by age
# for the 21day: 3, 0, 2, 10, 1
# for the 30days: 0, 1, 11, 6, 2
# for the 12days: 0, 1
# for the 6days: 3, 1, 0, 2
# for the 3days: 0, 2, 3

# DimPlot(rosette_21d, reduction = "umap", label = T)
# DimPlot(rosette_30d, reduction = "umap", label = T)
# DimPlot(seedling_12d, reduction = "umap", label = T)
# DimPlot(seedling_6d, reduction = "umap", label = T)
# DimPlot(seedling_3d, reduction = "umap", label = T)


# RIN4 AT3G25070

seurat_clusters_aggregates$rosette_21d['AT3G25070',] |> (\(x) sort(x, decreasing = T))()
#2, 1, 0, 3 
# sums to 5954.708
seurat_clusters_aggregates$rosette_30d['AT3G25070',] |> (\(x) sort(x, decreasing = T))()
#1, 4, 0, 5, 3
# 9609.956
seurat_clusters_aggregates$seedling_12d['AT3G25070',] |> (\(x) sort(x, decreasing = T))()
#0, 1, 2, 5, 3, 4 
# 11983.48
seurat_clusters_aggregates$seedling_6d['AT3G25070',] |> (\(x) sort(x, decreasing = T))()
#0, 2, 1, 3, 5
# 5490
seurat_clusters_aggregates$seedling_3d['AT3G25070',] |> (\(x) sort(x, decreasing = T))()
#1, 0, 3, 9, 8, 6, 2
# 1910


missing_ids_21 <- conserved_gene_ids[! conserved_gene_ids %in% rownames(seurat_clusters_aggregates$rosette_21d)]
missing_ids_30 <- conserved_gene_ids[! conserved_gene_ids %in% rownames(seurat_clusters_aggregates$rosette_30d)]
missing_ids_12 <- conserved_gene_ids[! conserved_gene_ids %in% rownames(seurat_clusters_aggregates$seedling_12d)]
missing_ids_6 <- conserved_gene_ids[! conserved_gene_ids %in% rownames(seurat_clusters_aggregates$seedling_6d)]
missing_ids_3 <- conserved_gene_ids[! conserved_gene_ids %in% rownames(seurat_clusters_aggregates$seedling_3d)]


# 1100 conserved ids are not in the rosette_21d data
# 1083 from rosette_30d
# 832 from seedling_12d
# 835 from seedling_6d
# 974 from seedling_3d

valid_ids_21 <- intersect(conserved_gene_ids, rownames(seurat_clusters_aggregates$rosette_21d))
valid_ids_30 <- intersect(conserved_gene_ids, rownames(seurat_clusters_aggregates$rosette_30d))
valid_ids_12 <- intersect(conserved_gene_ids, rownames(seurat_clusters_aggregates$seedling_12d))
valid_ids_6 <- intersect(conserved_gene_ids, rownames(seurat_clusters_aggregates$seedling_6d))
valid_ids_3 <- intersect(conserved_gene_ids, rownames(seurat_clusters_aggregates$seedling_3d))

guarded_log2 <- function(x) {
  return(log2(x + 1))
}

expressed_conserved_21 <- seurat_clusters_aggregates$rosette_21d[valid_ids_21, c('g2', 'g1', 'g0', 'g3')] |>
  apply(1, sum) |>
  guarded_log2() |>
  sort(decreasing = T)
expressed_conserved_30 <- seurat_clusters_aggregates$rosette_30d[valid_ids_30, c('g1', 'g4', 'g0', 'g5', 'g3')] |>
  apply(1, sum) |>
  guarded_log2() |>
  sort(decreasing = T)
expressed_conserved_12 <- seurat_clusters_aggregates$seedling_12d[valid_ids_12, c('g0', 'g1', 'g2', 'g5', 'g3', 'g4')] |>
  apply(1, sum) |>
  guarded_log2() |>
  sort(decreasing = T)
expressed_conserved_6 <- seurat_clusters_aggregates$seedling_6d[valid_ids_6, c('g0', 'g2', 'g1', 'g3', 'g5')] |>
  apply(1, sum) |>
  guarded_log2() |>
  sort(decreasing = T)
expressed_conserved_3 <- seurat_clusters_aggregates$seedling_3d[valid_ids_3, c('g1', 'g0', 'g3', 'g9', 'g8', 'g6', 'g2')] |>
  apply(1, sum) |>
  guarded_log2() |>
  sort(decreasing = T)


all_21 <- seurat_clusters_aggregates$rosette_21d[, c('g2', 'g1', 'g0', 'g3')] |>
  apply(1, sum) |>
  guarded_log2() |>
  sort(decreasing = T)
all_30 <- seurat_clusters_aggregates$rosette_30d[, c('g1', 'g4', 'g0', 'g5', 'g3')] |>
  apply(1, sum) |>
  guarded_log2() |>
  sort(decreasing = T)
all_12 <- seurat_clusters_aggregates$seedling_12d[, c('g0', 'g1', 'g2', 'g5', 'g3', 'g4')] |>
  apply(1, sum) |>
  guarded_log2() |>
  sort(decreasing = T)
all_6 <- seurat_clusters_aggregates$seedling_6d[, c('g0', 'g2', 'g1', 'g3', 'g5')] |>
  apply(1, sum) |>
  guarded_log2() |>
  sort(decreasing = T)
all_3 <- seurat_clusters_aggregates$seedling_3d[, c('g1', 'g0', 'g3', 'g9', 'g8', 'g6', 'g2')] |>
  apply(1, sum) |>
  guarded_log2() |>
  sort(decreasing = T)



nlr_list <- readLines('/Volumes/PhData2/SalkExpressionData/french_nlrs.txt')
laflamme_list <- c('AT3G07040', 'AT4G26090', 'AT5G45250', 'AT1G12220', 'AT5G46470', 'AT3G50950', 'AT5G18360', 'AT1G50180')

nlr_list <- c(nlr_list, laflamme_list)

nlr_list_21 <- nlr_list[nlr_list %in% rownames(seurat_clusters_aggregates$rosette_21d)]
nlr_list_30 <- nlr_list[nlr_list %in% rownames(seurat_clusters_aggregates$rosette_30d)]
nlr_list_12 <- nlr_list[nlr_list %in% rownames(seurat_clusters_aggregates$seedling_12d)]
nlr_list_6 <- nlr_list[nlr_list %in% rownames(seurat_clusters_aggregates$seedling_6d)]
nlr_list_3 <- nlr_list[nlr_list %in% rownames(seurat_clusters_aggregates$seedling_3d)]

nlr_expression_21 <- seurat_clusters_aggregates$rosette_21d[nlr_list_21, c('g2', 'g1', 'g0', 'g3')] |>
  apply(1, sum) |>
  guarded_log2() |>
  sort(decreasing = T)
nlr_expression_30 <- seurat_clusters_aggregates$rosette_30d[nlr_list_30, c('g1', 'g4', 'g0', 'g5', 'g3')] |>
  apply(1, sum) |>
  guarded_log2() |>
  sort(decreasing = T)
nlr_expression_12 <- seurat_clusters_aggregates$seedling_12d[nlr_list_12, c('g0', 'g1', 'g2', 'g5', 'g3', 'g4')] |>
  apply(1, sum) |>
  guarded_log2() |>
  sort(decreasing = T)
nlr_expression_6 <- seurat_clusters_aggregates$seedling_6d[nlr_list_6, c('g0', 'g2', 'g1', 'g3', 'g5')] |>
  apply(1, sum) |>
  guarded_log2() |>
  sort(decreasing = T)
nlr_expression_3 <- seurat_clusters_aggregates$seedling_3d[nlr_list_3, c('g1', 'g0', 'g3', 'g9', 'g8', 'g6', 'g2')] |>
  apply(1, sum) |>
  guarded_log2() |>
  sort(decreasing = T)

# make a ggplot of all_21 and nlr_expression_21 on the same axes -- make it a violin plot with points
all_genes <- sort(unique(c(names(all_3), names(all_6), names(all_12), names(all_21), names(all_30))))
get_expr <- function(vec, gene) {
  if (gene %in% names(vec)) return(vec[gene])
  else return(0)
}
# make all_data a dataframe with the names of the genes as the rownames, the columns are the ages, and 0 if the gene is not expressed
all_data_df <- data.frame(
  all_3  = sapply(all_genes, get_expr, vec = all_3),
  all_6  = sapply(all_genes, get_expr, vec = all_6),
  all_12 = sapply(all_genes, get_expr, vec = all_12),
  all_21 = sapply(all_genes, get_expr, vec = all_21),
  all_30 = sapply(all_genes, get_expr, vec = all_30)
)
new_names <- sapply(rownames(all_data_df), function(x) strsplit(x, "\\.")[[1]][1])
rownames(all_data_df) <- new_names


expressed_conserved_nlrs_21 <- expressed_conserved_21[nlr_list]
expressed_conserved_nlrs_21 <- expressed_conserved_nlrs_21[!is.na(expressed_conserved_nlrs_21)]
expressed_conserved_nlrs_30 <- expressed_conserved_30[nlr_list]
expressed_conserved_nlrs_30 <- expressed_conserved_nlrs_30[!is.na(expressed_conserved_nlrs_30)]
expressed_conserved_nlrs_12 <- expressed_conserved_12[nlr_list]
expressed_conserved_nlrs_12 <- expressed_conserved_nlrs_12[!is.na(expressed_conserved_nlrs_12)]
expressed_conserved_nlrs_6 <- expressed_conserved_6[nlr_list]
expressed_conserved_nlrs_6 <- expressed_conserved_nlrs_6[!is.na(expressed_conserved_nlrs_6)]
expressed_conserved_nlrs_3 <- expressed_conserved_3[nlr_list]
expressed_conserved_nlrs_3 <- expressed_conserved_nlrs_3[!is.na(expressed_conserved_nlrs_3)]

# this one doesn't use the expressed_conserved_nlrs_ variables, but rather the nlr_expression_ variables
# it's a larger set, not sure if it's better or worse

#expressed_conserved_nlrs_21 <- expressed_conserved_21[nlr_list]
#expressed_conserved_nlrs_21 <- expressed_conserved_21[!is.na(expressed_conserved_nlrs_21)]
#expressed_conserved_nlrs_30 <- expressed_conserved_30[nlr_list]
#expressed_conserved_nlrs_30 <- expressed_conserved_30[!is.na(expressed_conserved_nlrs_30)]
#expressed_conserved_nlrs_12 <- expressed_conserved_12[nlr_list]
#expressed_conserved_nlrs_12 <- expressed_conserved_12[!is.na(expressed_conserved_nlrs_12)]
#expressed_conserved_nlrs_6 <- expressed_conserved_6[nlr_list]
#expressed_conserved_nlrs_6 <- expressed_conserved_6[!is.na(expressed_conserved_nlrs_6)]
#expressed_conserved_nlrs_3 <- expressed_conserved_3[nlr_list]
#expressed_conserved_nlrs_3 <- expressed_conserved_3[!is.na(expressed_conserved_nlrs_3)]


# create a dataframe of the nlr expression levels
nlr_expression_df <- data.frame(
  gene = nlr_list,
  seedling_3d = nlr_expression_3[nlr_list],
  seedling_6d = nlr_expression_6[nlr_list],
  seedling_12d = nlr_expression_12[nlr_list],
  rosette_21d = nlr_expression_21[nlr_list],
  rosette_30d = nlr_expression_30[nlr_list]
)

nlr_expression_df <- nlr_expression_df[rowSums(is.na(nlr_expression_df)) != ncol(nlr_expression_df), ]
nlr_expression_df[is.na(nlr_expression_df)] <- 0
nlr_expression_df <- nlr_expression_df[!duplicated(nlr_expression_df$gene), ]
rownames(nlr_expression_df) <- nlr_expression_df$gene
#nlr_expression_df <- nlr_expression_df[, -1]

## ------------------------------------------------------------
## Column-wise means of all_data_df and NLR z-scores
## ------------------------------------------------------------
## We compute per-age (column) means and standard deviations across all genes
## in all_data_df, then scale each NLR expression value relative to those
## background distributions (same age). Result: nlr_expression_zscores_df.

all_col_means <- colMeans(all_data_df, na.rm = TRUE)
all_col_sds   <- sapply(all_data_df, sd, na.rm = TRUE)

# Map NLR dataframe column names to all_data_df columns
nlr_cols  <- c("seedling_3d","seedling_6d","seedling_12d","rosette_21d","rosette_30d")
all_cols  <- c("all_3","all_6","all_12","all_21","all_30")

# Ensure all expected columns exist before proceeding
stop_if_missing <- function(required, present, label){
  missing <- setdiff(required, present)
  if (length(missing) > 0) stop(paste("Missing", label, "columns:", paste(missing, collapse=", ")))
}
stop_if_missing(all_cols,   colnames(all_data_df),       "all_data_df")
stop_if_missing(nlr_cols,   colnames(nlr_expression_df), "nlr_expression_df")

z_mat <- sapply(seq_along(nlr_cols), function(i){
  nlr_col <- nlr_cols[i]
  all_col <- all_cols[i]
  mu  <- all_col_means[all_col]
  sdv <- all_col_sds[all_col]
  if (is.na(sdv) || sdv == 0) {
    rep(NA_real_, nrow(nlr_expression_df))
  } else {
    (nlr_expression_df[[nlr_col]] - mu) / sdv
  }
})
colnames(z_mat) <- nlr_cols
rownames(z_mat) <- rownames(nlr_expression_df)

nlr_expression_zscores_df <- as.data.frame(z_mat)
nlr_expression_zscores_df$gene <- rownames(nlr_expression_zscores_df)
nlr_expression_zscores_df <- nlr_expression_zscores_df[, c("gene", nlr_cols)]

# (Optional) preview summary
nlr_z_summary <- sapply(nlr_expression_zscores_df[ , nlr_cols], summary)

## ------------------------------------------------------------
## Long-format z-scores + row-wise mean sorting + heatmap
## ------------------------------------------------------------
library(dplyr)
library(tidyr)
underscore_label <- function(vals) sub('^[^_]+_', '', vals)

# Compute row-wise mean z-score across ages
nlr_expression_zscores_df$row_mean <- rowMeans(nlr_expression_zscores_df[, nlr_cols], na.rm = TRUE)

nlr_expression_zscores_df_ordered <- nlr_expression_zscores_df |>
  arrange(desc(row_mean))

write.csv(nlr_expression_df, file = "nlr_expression_unordered.csv", row.names = FALSE)

write.csv(nlr_expression_zscores_df_ordered, file = "nlr_expression_zscores_ordered.csv", row.names = FALSE)

# Order genes by descending average z-score
ordered_genes <- nlr_expression_zscores_df |>
  arrange(row_mean) |>
  pull(gene)

# Ensure gene column exists (if object constructed in a fresh session)
if (!'gene' %in% colnames(nlr_expression_zscores_df)) {
  nlr_expression_zscores_df$gene <- rownames(nlr_expression_zscores_df)
}

# Long format for plotting (namespace dplyr/tidyr to avoid select() masking)
nlr_z_long <- nlr_expression_zscores_df |>
  dplyr::select(gene, dplyr::all_of(nlr_cols), row_mean) |>
  tidyr::pivot_longer(cols = dplyr::all_of(nlr_cols), names_to = "age", values_to = "zscore")

# Set factor levels for consistent x (age) order and y (gene) order (high mean on top)
nlr_z_long$age  <- factor(nlr_z_long$age, levels = nlr_cols)
nlr_z_long$gene <- factor(nlr_z_long$gene, levels = ordered_genes)

# -------------------------------------------------------------
# Optional: add a left marker column with blue stars for laflamme_list genes
# without disturbing the main age columns. We prepend a pseudo-column 'marker'.
# -------------------------------------------------------------
if (exists('laflamme_list')) {
  marker_col <- 'laflamme'
  age_levels <- levels(nlr_z_long$age)
  # Only add marker level if it is not already present
  if (!(marker_col %in% age_levels)) {
    new_levels <- c(marker_col, age_levels)
    nlr_z_long$age <- factor(nlr_z_long$age, levels = new_levels)
    # Create empty marker rows (zscore NA) to reserve the column space
    marker_rows <- data.frame(
      gene     = ordered_genes,
      age      = factor(rep(marker_col, length(ordered_genes)), levels = new_levels),
      zscore   = NA_real_,
      row_mean = nlr_expression_zscores_df$row_mean[match(ordered_genes, nlr_expression_zscores_df$gene)]
    )
    nlr_z_long <- rbind(marker_rows, nlr_z_long)
  }
  genes_with_star <- intersect(laflamme_list, ordered_genes)
  if (length(genes_with_star) > 0) {
    star_df <- data.frame(
      gene    = genes_with_star,
      age     = factor(rep(marker_col, length(genes_with_star)), levels = levels(nlr_z_long$age)),
      zscore  = NA_real_,
      row_mean = nlr_expression_zscores_df$row_mean[match(genes_with_star, nlr_expression_zscores_df$gene)]
    )
  } else {
    star_df <- data.frame(
      gene    = character(0),
      age     = factor(character(0), levels = levels(nlr_z_long$age)),
      zscore  = numeric(0),
      row_mean = numeric(0)
    )
  }
} else {
  star_df <- NULL
}

# Heatmap of z-scores: red = most negative, blue = most positive (with optional marker column)
nlr_z_heatmap <- ggplot(nlr_z_long, aes(x = age, y = gene, fill = zscore)) +
  geom_tile(color = NA) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0,
                       name = "Z-score", na.value = "white") +
  labs(title = "NLR Expression Z-scores (Background = All Genes)",
       x = "Age / Condition",
       y = "NLR Gene (sorted by mean z-score)") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 4, angle = 45),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_x_discrete(labels = underscore_label)

if (!is.null(star_df) && nrow(star_df) > 0) {
  nlr_z_heatmap <- nlr_z_heatmap +
    geom_text(data = star_df, aes(x = age, y = gene), label = '*', color = 'blue', size = 2, na.rm = TRUE)
}
  scale_x_discrete(labels = underscore_label)

# Print plot (comment out if sourcing)
print(nlr_z_heatmap)




## plotting ordinary all gene expression and then nlr expression on top
library(tidyr)

long_data <- pivot_longer(all_data_df,
                          cols = everything(),
                          names_to = "age",
                          values_to = "expression")

long_data$age <- factor(long_data$age, levels = c("all_3", "all_6", "all_12", "all_21", "all_30"))

ggplot(long_data, aes(x = age, y = expression)) +
  geom_violin(fill = "lightblue", alpha = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red") +
  labs(title = "Gene Expression by Age",
       x = "Age",
       y = "Expression (log2)") +
  theme_minimal() +
  scale_x_discrete(labels = underscore_label)

nlr_long <- pivot_longer(nlr_expression_df,
                         cols = -gene,
                         names_to = "age",
                         values_to = "expression")

library(dplyr)

nlr_long$age <- recode(nlr_long$age,
                       "seedling_3d" = "all_3",
                       "seedling_6d" = "all_6",
                       "seedling_12d" = "all_12",
                       "rosette_21d" = "all_21",
                       "rosette_30d" = "all_30")

nlr_long$age <- factor(nlr_long$age, levels = c("all_3", "all_6", "all_12", "all_21", "all_30"))

ggplot(long_data, aes(x = age, y = expression)) +
  geom_violin(fill = "lightblue", alpha = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red") +
  geom_point(data = nlr_long, 
             aes(x = age, y = expression), 
             color = "orange", 
             size = 2,
             position = position_jitter(width = 0.1)) +
  labs(title = "Gene Expression by Age",
       x = "Age",
       y = "Expression (log2)") +
  theme_minimal() +
  scale_x_discrete(labels = underscore_label)

guardee_list <- c('AT3G25070', 'AT5G13160', 'AT3G57750', 'AT1G14370')

# create a dataframe of the guardee expression levels
guardee_expression_df <- data.frame(
  gene = guardee_list,
  all_3 = sapply(guardee_list, function(gene) get_expr(all_3, gene)),
  all_6 = sapply(guardee_list, function(gene) get_expr(all_6, gene)),
  all_12 = sapply(guardee_list, function(gene) get_expr(all_12, gene)),
  all_21 = sapply(guardee_list, function(gene) get_expr(all_21, gene)),
  all_30 = sapply(guardee_list, function(gene) get_expr(all_30, gene))
)
new_guardee_names <- sapply(rownames(guardee_expression_df), function(x) strsplit(x, "\\.")[[1]][1])
rownames(guardee_expression_df) <- new_guardee_names

# Transform the guardee expression dataframe into long format.
guardee_long <- pivot_longer(guardee_expression_df, 
                             cols = -gene, 
                             names_to = "age", 
                             values_to = "expression")

# Ensure the factor levels are in the desired order:
guardee_long$age <- factor(guardee_long$age, levels = c("all_3", "all_6", "all_12", "all_21", "all_30"))

# Now add guardee points in violet to the ggplot:
ggplot(long_data, aes(x = age, y = expression)) +
  geom_violin(fill = "lightblue", alpha = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red") +
  # Green points for nlr expression data:
  geom_point(data = nlr_long, 
             aes(x = age, y = expression), 
             color = "green", 
             size = 2,
             position = position_jitter(width = 0.1)) +
  # Violet points for guardee expression data:
  geom_point(data = guardee_long, 
             aes(x = age, y = expression), 
             color = "violet", 
             size = 2,
             position = position_jitter(width = 0.1)) +
  labs(title = "Gene Expression by Age",
       x = "Age",
       y = "Expression (log2)") +
  theme_minimal() +
  scale_x_discrete(labels = underscore_label)

cognate_list <- c('AT3G07040', 'AT1G12220', 'AT3G50950', 'AT3G57710')

cognate_expression_df <- data.frame(
  gene = cognate_list,
  all_3 = sapply(cognate_list, function(gene) get_expr(all_3, gene)),
  all_6 = sapply(cognate_list, function(gene) get_expr(all_6, gene)),
  all_12 = sapply(cognate_list, function(gene) get_expr(all_12, gene)),
  all_21 = sapply(cognate_list, function(gene) get_expr(all_21, gene)),
  all_30 = sapply(cognate_list, function(gene) get_expr(all_30, gene))
)
new_cognate_names <- sapply(rownames(cognate_expression_df), function(x) strsplit(x, "\\.")[[1]][1])
rownames(cognate_expression_df) <- new_cognate_names

cognate_long <- pivot_longer(cognate_expression_df, 
                             cols = -gene, 
                             names_to = "age", 
                             values_to = "expression")

# Ensure the factor levels are in the desired order:
cognate_long$age <- factor(cognate_long$age, levels = c("all_3", "all_6", "all_12", "all_21", "all_30"))

# Now add cognate points in orange to the ggplot:
ggplot(long_data, aes(x = age, y = expression)) +
  geom_violin(fill = "lightblue", alpha = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red") +
  # Green points for nlr expression data:
  geom_point(data = nlr_long, 
             aes(x = age, y = expression), 
             color = "orange", 
             size = 2,
             position = position_jitter(width = 0.1)) +
  # Violet points for guardee expression data:
  geom_point(data = guardee_long, 
             aes(x = age, y = expression), 
             color = "purple", 
             size = 2,
             position = position_jitter(width = 0.1)) +
  # Orange points for cognate expression data:
  geom_point(data = cognate_long, 
             aes(x = age, y = expression), 
             color = "blue", 
             size = 2,
             position = position_jitter(width = 0.1)) +
  labs(title = "Gene Expression by Age",
       x = "Age",
       y = "Expression (log2)") +
  theme_minimal() +
  scale_x_discrete(labels = underscore_label)


laflamme_list <- c('AT3G07040', 'AT4G26090', 'AT5G45250', 'AT1G12220', 'AT5G46470', 'AT3G50950', 'AT5G18360', 'AT1G50180')

laflamme_expression_df <- data.frame(
  gene = laflamme_list,
  all_3 = sapply(laflamme_list, function(gene) get_expr(all_3, gene)),
  all_6 = sapply(laflamme_list, function(gene) get_expr(all_6, gene)),
  all_12 = sapply(laflamme_list, function(gene) get_expr(all_12, gene)),
  all_21 = sapply(laflamme_list, function(gene) get_expr(all_21, gene)),
  all_30 = sapply(laflamme_list, function(gene) get_expr(all_30, gene))
)
new_laflamme_names <- sapply(rownames(laflamme_expression_df), function(x) strsplit(x, "\\.")[[1]][1])
rownames(laflamme_expression_df) <- new_laflamme_names

laflamme_long <- pivot_longer(laflamme_expression_df, 
                              cols = -gene, 
                              names_to = "age", 
                              values_to = "expression")

# Ensure the factor levels are in the desired order:
laflamme_long$age <- factor(laflamme_long$age, levels = c("all_3", "all_6", "all_12", "all_21", "all_30"))

# Now add cognate points in orange to the ggplot:
ggplot(long_data, aes(x = age, y = expression)) +
  geom_violin(fill = "lightblue", alpha = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red") +
  # Green points for nlr expression data:
  geom_point(data = nlr_long, 
             aes(x = age, y = expression), 
             color = "orange", 
             size = 2,
             position = position_jitter(width = 0.1)) +
  # Orange points for cognate expression data:
  geom_point(data = laflamme_long, 
             aes(x = age, y = expression), 
             color = "blue", 
             size = 2,
             position = position_jitter(width = 0.1)) +
  labs(title = "Gene Expression by Age",
       x = "Age",
       y = "Expression (log2)") +
  theme_minimal() +
  scale_x_discrete(labels = underscore_label)

# initial double violin

ggplot(long_data, aes(x = age, y = expression)) +
  geom_violin(fill = "lightblue", alpha = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red") +
  geom_violin(data = nlr_long, 
              aes(x = age, y = expression), 
              fill = "orange", alpha = 0.5, show.legend = FALSE) +
  # Green points for nlr expression data:
  # Orange points for cognate expression data:
  geom_point(data = laflamme_long, 
             aes(x = age, y = expression), 
             color = "blue", 
             size = 2,
             position = position_jitter(width = 0.1)) +
  labs(title = "Gene Expression by Age",
       x = "Age",
       y = "Expression (log2)") +
  theme_minimal() +
  scale_x_discrete(labels = underscore_label)


# create a dataframe and do the same with the expressed_conserved_nlrs



# Calculate scaling factors for violin widths based on histogram frequencies (not max values)
# Get histogram data for background
background_hist <- hist(long_data$expression, plot = FALSE, breaks = 30)
background_max <- max(background_hist$counts)

# Get histogram data for NLR
nlr_hist <- hist(nlr_long$expression, plot = FALSE, breaks = 30)
nlr_max <- max(nlr_hist$counts)


# Use log scaling for width proportionality
log_background_max <- log10(background_max + 1)
log_nlr_max <- log10(nlr_max + 1)
width_scale_factor <- log_nlr_max / log_background_max

# Create the plot with scaled violin widths
library(ggplot2)
library(scales)
library(ggtext)

# for raman_incomplete create dfs with an expression column pulled with get_expr from all_data_df and filter out the incomplete ones

raman_incomplete_expression_df <- data.frame(
  gene = raman_incomplete$gene,
  incomplete = raman_incomplete$incomplete,
  all_3 = sapply(raman_incomplete$gene, function(gene) get_expr(all_3, gene)),
  all_6 = sapply(raman_incomplete$gene, function(gene) get_expr(all_6, gene)),
  all_12 = sapply(raman_incomplete$gene, function(gene) get_expr(all_12, gene)),
  all_21 = sapply(raman_incomplete$gene, function(gene) get_expr(all_21, gene)),
  all_30 = sapply(raman_incomplete$gene, function(gene) get_expr(all_30, gene))
)
new_raman_names <- sapply(rownames(raman_incomplete_expression_df), function(x) strsplit(x, "\\.")[[1]][1])
rownames(raman_incomplete_expression_df) <- new_raman_names

raman_long <- pivot_longer(raman_incomplete_expression_df, 
                              cols = -c(gene, incomplete), 
                              names_to = "age", 
                              values_to = "expression")

# Ensure the factor levels are in the desired order:
raman_long$age <- factor(raman_long$age, levels = c("all_3", "all_6", "all_12", "all_21", "all_30"))



ggplot(long_data, aes(x = age, y = expression)) +
  geom_violin(fill = "lightblue", alpha = 1, width = 1) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red") +
  stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "black") +
  geom_violin(data = nlr_long, 
              aes(x = age, y = expression), 
              fill = "orange", 
              alpha = 0.7, 
              width = width_scale_factor,  # Scale the width proportionally
              show.legend = FALSE) +
  stat_summary(data = nlr_long, fun = mean, geom = "crossbar", width = 0.3, color = "red", linetype = "dotted") +
  stat_summary(data = nlr_long, fun = median, geom = "crossbar", width = 0.3, color = "black", linetype = "dotted") +
  geom_point(data = raman_long, 
             aes(x = age, y = expression, color = incomplete), 
             size = 0.5,
             position = position_jitter(width = 0.1)) +
  scale_x_discrete(labels = c("all_3" = "3 days", 
                              "all_6" = "6 days", 
                              "all_12" = "12 days", 
                              "all_21" = "21 days", 
                              "all_30" = "30 days")) +
  labs(title = "Gene Expression by Age",
       subtitle = sprintf("Histogram Frequencies (log10 compressed)- <br><span style='color:lightblue'>All genes max: %d</span>,<span style='color:orange'>NLR max: %d</span> (width ratio: %.2f)",
                          background_max, nlr_max, width_scale_factor),
       x = "Age",
       y = "Expression (log2)") +
  theme_minimal() +
  theme(
    plot.margin = margin(t = 40, r = 5, b = 5, l = 5, "pt"),  # Add top margin for second axis
    plot.subtitle = ggtext::element_markdown(hjust = 0.5, margin = margin(b = 10))  # Enable HTML rendering
  )

# for raman_long, create a subset for each day for every gene above median expression and then one below median
raman_long_3_above <- raman_long[raman_long$age == 'all_3' & raman_long$expression > median(long_data$expression), ]
raman_long_3_below <- raman_long[raman_long$age == 'all_3' & raman_long$expression <= median(long_data$expression), ]
raman_long_6_above <- raman_long[raman_long$age == 'all_6' & raman_long$expression > median(long_data$expression), ]
raman_long_6_below <- raman_long[raman_long$age == 'all_6' & raman_long$expression <= median(long_data$expression), ]
raman_long_12_above <- raman_long[raman_long$age == 'all_12' & raman_long$expression > median(long_data$expression), ]
raman_long_12_below <- raman_long[raman_long$age == 'all_12' & raman_long$expression <= median(long_data$expression), ]
raman_long_21_above <- raman_long[raman_long$age == 'all_21' & raman_long$expression > median(long_data$expression), ]
raman_long_21_below <- raman_long[raman_long$age == 'all_21' & raman_long$expression <= median(long_data$expression), ]
raman_long_30_above <- raman_long[raman_long$age == 'all_30' & raman_long$expression > median(long_data$expression), ]
raman_long_30_below <- raman_long[raman_long$age == 'all_30' & raman_long$expression <= median(long_data$expression), ]
# count how many incomplete in each, store as data frame with the age, above/below median, and count of incomplete
incomplete_counts_df <- data.frame(
  age = c("3d", "3d", "6d", "6d", "12d", "12d", "21d", "21d", "30d", "30d"),
  expression_level = c("above_median", "below_median", "above_median", "below_median", 
                       "above_median", "below_median", "above_median", "below_median", 
                       "above_median", "below_median"),
  incomplete_count = c(
    sum(raman_long_3_above$incomplete),
    sum(raman_long_3_below$incomplete),
    sum(raman_long_6_above$incomplete),
    sum(raman_long_6_below$incomplete),
    sum(raman_long_12_above$incomplete),
    sum(raman_long_12_below$incomplete),
    sum(raman_long_21_above$incomplete),
    sum(raman_long_21_below$incomplete),
    sum(raman_long_30_above$incomplete),
    sum(raman_long_30_below$incomplete)
  ),
  total_count = c(
    nrow(raman_long_3_above),
    nrow(raman_long_3_below),
    nrow(raman_long_6_above),
    nrow(raman_long_6_below),
    nrow(raman_long_12_above),
    nrow(raman_long_12_below),
    nrow(raman_long_21_above),
    nrow(raman_long_21_below),
    nrow(raman_long_30_above),
    nrow(raman_long_30_below)
  )
)

# Calculate proportion of incomplete genes
incomplete_counts_df$incomplete_proportion <- incomplete_counts_df$incomplete_count / incomplete_counts_df$total_count

print(incomplete_counts_df)

# plot the proportions as a bar plot with age on the x axis, proportion on the y axis, and fill by above/below median
ggplot(incomplete_counts_df, aes(x = age, y = incomplete_proportion, fill = expression_level)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Proportion of Incomplete Genes by Age and Expression Level",
       x = "Age",
       y = "Proportion of Incomplete Genes",
       fill = "Expression Level") +
  theme_minimal()

# generate csvs of nlr_long for each age with the gene names and expression levels where expression is above median, and where expression is below median and > 0, and where expression is 0
write.csv(nlr_long[nlr_long$age == 'all_3' & nlr_long$expression > median(long_data$expression), c('gene', 'expression')], 
          file = 'nlr_expression_3d_above_median.csv', row.names = FALSE)
write.csv(nlr_long[nlr_long$age == 'all_3' & nlr_long$expression <= median(long_data$expression) & nlr_long$expression > 0, c('gene', 'expression')],
          file = 'nlr_expression_3d_below_median.csv', row.names = FALSE)
write.csv(nlr_long[nlr_long$age == 'all_3' & nlr_long$expression == 0, c('gene', 'expression')],
          file = 'nlr_expression_3d_zero.csv', row.names = FALSE)
write.csv(nlr_long[nlr_long$age == 'all_6' & nlr_long$expression > median(long_data$expression), c('gene', 'expression')],
          file = 'nlr_expression_6d_above_median.csv', row.names = FALSE)
write.csv(nlr_long[nlr_long$age == 'all_6' & nlr_long$expression <= median(long_data$expression) & nlr_long$expression > 0, c('gene', 'expression')],
          file = 'nlr_expression_6d_below_median.csv', row.names = FALSE)
write.csv(nlr_long[nlr_long$age == 'all_6' & nlr_long$expression == 0, c('gene', 'expression')],
          file = 'nlr_expression_6d_zero.csv', row.names = FALSE)
write.csv(nlr_long[nlr_long$age == 'all_12' & nlr_long$expression > median(long_data$expression), c('gene', 'expression')],
          file = 'nlr_expression_12d_above_median.csv', row.names = FALSE)
write.csv(nlr_long[nlr_long$age == 'all_12' & nlr_long$expression <= median(long_data$expression) & nlr_long$expression > 0, c('gene', 'expression')],
          file = 'nlr_expression_12d_below_median.csv', row.names = FALSE)
write.csv(nlr_long[nlr_long$age == 'all_12' & nlr_long$expression == 0, c('gene', 'expression')],
          file = 'nlr_expression_12d_zero.csv', row.names = FALSE)
write.csv(nlr_long[nlr_long$age == 'all_21' & nlr_long$expression > median(long_data$expression), c('gene', 'expression')],
          file = 'nlr_expression_21d_above_median.csv', row.names = FALSE)
write.csv(nlr_long[nlr_long$age == 'all_21' & nlr_long$expression <= median(long_data$expression) & nlr_long$expression > 0, c('gene', 'expression')],
          file = 'nlr_expression_21d_below_median.csv', row.names = FALSE)
write.csv(nlr_long[nlr_long$age == 'all_21' & nlr_long$expression == 0, c('gene', 'expression')],
          file = 'nlr_expression_21d_zero.csv', row.names = FALSE)
write.csv(nlr_long[nlr_long$age == 'all_30' & nlr_long$expression > median(long_data$expression), c('gene', 'expression')],
          file = 'nlr_expression_30d_above_median.csv', row.names = FALSE)
write.csv(nlr_long[nlr_long$age == 'all_30' & nlr_long$expression <= median(long_data$expression) & nlr_long$expression > 0, c('gene', 'expression')],
          file = 'nlr_expression_30d_below_median.csv', row.names = FALSE)
write.csv(nlr_long[nlr_long$age == 'all_30' & nlr_long$expression == 0, c('gene', 'expression')],
          file = 'nlr_expression_30d_zero.csv', row.names = FALSE)
# do the same for raman_long
write.csv(raman_long[raman_long$age == 'all_3' & raman_long$expression > median(long_data$expression), c('gene', 'expression')], 
          file = 'raman_expression_3d_above_median.csv', row.names = FALSE)
write.csv(raman_long[raman_long$age == 'all_3' & raman_long$expression <= median(long_data$expression) & raman_long$expression > 0, c('gene', 'expression')],
          file = 'raman_expression_3d_below_median.csv', row.names = FALSE)
write.csv(raman_long[raman_long$age == 'all_3' & raman_long$expression == 0, c('gene', 'expression')],
          file = 'raman_expression_3d_zero.csv', row.names = FALSE)
write.csv(raman_long[raman_long$age == 'all_6' & raman_long$expression > median(long_data$expression), c('gene', 'expression')],
          file = 'raman_expression_6d_above_median.csv', row.names = FALSE)
write.csv(raman_long[raman_long$age == 'all_6' & raman_long$expression <= median(long_data$expression) & raman_long$expression > 0, c('gene', 'expression')],
          file = 'raman_expression_6d_below_median.csv', row.names = FALSE)
write.csv(raman_long[raman_long$age == 'all_6' & raman_long$expression == 0, c('gene', 'expression')],
          file = 'raman_expression_6d_zero.csv', row.names = FALSE)
write.csv(raman_long[raman_long$age == 'all_12' & raman_long$expression > median(long_data$expression), c('gene', 'expression')],
          file = 'raman_expression_12d_above_median.csv', row.names = FALSE)
write.csv(raman_long[raman_long$age == 'all_12' & raman_long$expression <= median(long_data$expression) & raman_long$expression > 0, c('gene', 'expression')],
          file = 'raman_expression_12d_below_median.csv', row.names = FALSE)
write.csv(raman_long[raman_long$age == 'all_12' & raman_long$expression == 0, c('gene', 'expression')],
          file = 'raman_expression_12d_zero.csv', row.names = FALSE)
write.csv(raman_long[raman_long$age == 'all_21' & raman_long$expression > median(long_data$expression), c('gene', 'expression')],
          file = 'raman_expression_21d_above_median.csv', row.names = FALSE)
write.csv(raman_long[raman_long$age == 'all_21' & raman_long$expression <= median(long_data$expression) & raman_long$expression > 0, c('gene', 'expression')],
          file = 'raman_expression_21d_below_median.csv', row.names = FALSE)
write.csv(raman_long[raman_long$age == 'all_21' & raman_long$expression == 0, c('gene', 'expression')],
          file = 'raman_expression_21d_zero.csv', row.names = FALSE)
write.csv(raman_long[raman_long$age == 'all_30' & raman_long$expression > median(long_data$expression), c('gene', 'expression')],
          file = 'raman_expression_30d_above_median.csv', row.names = FALSE)
write.csv(raman_long[raman_long$age == 'all_30' & raman_long$expression <= median(long_data$expression) & raman_long$expression > 0, c('gene', 'expression')],
          file = 'raman_expression_30d_below_median.csv', row.names = FALSE)
write.csv(raman_long[raman_long$age == 'all_30' & raman_long$expression == 0, c('gene', 'expression')],
          file = 'raman_expression_30d_zero.csv', row.names = FALSE)

# Add subcat_Chogag column to raman_long by matching gene names
raman_long$subcat_Choghag <- sapply(raman_long$gene, function(gene) {
  match_idx <- which(raman_stuff$geneID == gene)
  if (length(match_idx) == 0) return(NA_character_)
  else return(as.character(raman_stuff$subcat_Choghag[match_idx[1]]))  # Force to character
}, USE.NAMES = FALSE)

ggplot(long_data, aes(x = age, y = expression)) +
  geom_violin(fill = "lightblue", alpha = 1, width = 1) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red") +
  stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "black") +
  geom_violin(data = nlr_long, 
              aes(x = age, y = expression), 
              fill = "orange", 
              alpha = 0.7, 
              width = width_scale_factor,  # Scale the width proportionally
              show.legend = FALSE) +
  stat_summary(data = nlr_long, fun = mean, geom = "crossbar", width = 0.3, color = "red", linetype = "dotted") +
  stat_summary(data = nlr_long, fun = median, geom = "crossbar", width = 0.3, color = "black", linetype = "dotted") +
  geom_point(data = raman_long, 
             aes(x = age, y = expression, color = subcat_Choghag), 
             size = 0.5,
             position = position_jitter(width = 0.1)) +
  scale_x_discrete(labels = c("all_3" = "3 days", 
                              "all_6" = "6 days", 
                              "all_12" = "12 days", 
                              "all_21" = "21 days", 
                              "all_30" = "30 days")) +
  labs(title = "Gene Expression by Age",
       subtitle = sprintf("Histogram Frequencies (log10 compressed)- <br><span style='color:lightblue'>All genes max: %d</span>,<span style='color:orange'>NLR max: %d</span> (width ratio: %.2f)",
                          background_max, nlr_max, width_scale_factor),
       x = "Age",
       y = "Expression (log2)") +
  theme_minimal() +
  theme(
    plot.margin = margin(t = 40, r = 5, b = 5, l = 5, "pt"),  # Add top margin for second axis
    plot.subtitle = ggtext::element_markdown(hjust = 0.5, margin = margin(b = 10))  # Enable HTML rendering
  )

# remove subcat_Chogag column from raman_long
# raman_long$subcat_Chogag <- NULL

ggplot(long_data, aes(x = age, y = expression)) +
  geom_violin(fill = "lightblue", alpha = 1, width = 1) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red") +
  stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "black") +
  geom_violin(data = raman_long, 
              aes(x = age, y = expression, fill = subcat_Choghag), 
              #fill = "orange", 
              alpha = 0.7, 
              #width = width_scale_factor,  # Scale the width proportionally
              show.legend = TRUE) +
#  stat_summary(data = nlr_long, fun = mean, geom = "crossbar", width = 0.3, color = "red", linetype = "dotted") +
#  stat_summary(data = nlr_long, fun = median, geom = "crossbar", width = 0.3, color = "black", linetype = "dotted") +
#  geom_point(data = raman_long, 
#             aes(x = age, y = expression, color = subcat_Choghag), 
#             size = 0.5,
#             position = position_jitter(width = 0.1)) +
  scale_x_discrete(labels = c("all_3" = "3 days", 
                              "all_6" = "6 days", 
                              "all_12" = "12 days", 
                              "all_21" = "21 days", 
                              "all_30" = "30 days")) +
  labs(title = "Gene Expression by Age",
       subtitle = sprintf("Histogram Frequencies (log10 compressed)- <br><span style='color:lightblue'>All genes max: %d</span>,<span style='color:orange'>NLR max: %d</span> (width ratio: %.2f)",
                          background_max, nlr_max, width_scale_factor),
       x = "Age",
       y = "Expression (log2)") +
  theme_minimal() +
  theme(
    plot.margin = margin(t = 40, r = 5, b = 5, l = 5, "pt"),  # Add top margin for second axis
    plot.subtitle = ggtext::element_markdown(hjust = 0.5, margin = margin(b = 10))  # Enable HTML rendering
  )


# boxplot version
ggplot(long_data, aes(x = age, y = expression)) +
  geom_boxplot(fill = "lightblue", alpha = 1, width = 0.5, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red") +
  stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "black") +
  geom_boxplot(data = raman_long, 
              aes(x = age, y = expression, fill = subcat_Choghag), 
              #fill = "orange", 
              alpha = 0.7, 
              width = 0.5,
              position = position_dodge(width = 0.75),
              outlier.shape = NA,
              show.legend = TRUE) +
  scale_x_discrete(labels = c("all_3" = "3 days", 
                              "all_6" = "6 days", 
                              "all_12" = "12 days", 
                              "all_21" = "21 days", 
                              "all_30" = "30 days")) +
  labs(title = "Gene Expression by Age",
       subtitle = sprintf("Histogram Frequencies (log10 compressed)- <br><span style='color:lightblue'>All genes max: %d</span>,<span style='color:orange'>NLR max: %d</span> (width ratio: %.2f)",
                          background_max, nlr_max, width_scale_factor),
       x = "Age",
       y = "Expression (log2)") +
  theme_minimal() +
  theme(
    plot.margin = margin(t = 40, r = 5, b = 5, l = 5, "pt"),  # Add top margin for second axis
    plot.subtitle = ggtext::element_markdown(hjust = 0.5, margin = margin(b = 10))  # Enable HTML rendering
  )






raman_long_full_nlrs <- raman_long %>% filter(!(subcat_Choghag %in% c('TIR', 'NL', NA, 'LR')))

ggplot(long_data, aes(x = age, y = expression)) +
  geom_boxplot(fill = "lightblue", alpha = 1, width = 0.5, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red") +
  stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "black") +
  geom_boxplot(data = raman_long_full_nlrs, 
               aes(x = age, y = expression, fill = subcat_Choghag), 
               alpha = 0.7, 
               width = 0.5,
               position = position_dodge(width = 0.75),
               outlier.shape = NA,
               show.legend = TRUE) +
  scale_x_discrete(labels = c("all_3" = "3 days", 
                              "all_6" = "6 days", 
                              "all_12" = "12 days", 
                              "all_21" = "21 days", 
                              "all_30" = "30 days")) +
  labs(title = "Gene Expression by Age",
       subtitle = sprintf("Histogram Frequencies (log10 compressed)- <br><span style='color:lightblue'>All genes max: %d</span>,<span style='color:orange'>NLR max: %d</span> (width ratio: %.2f)",
                          background_max, nlr_max, width_scale_factor),
       x = "Age",
       y = "Expression (log2)") +
  theme_minimal() +
  theme(
    plot.margin = margin(t = 40, r = 5, b = 5, l = 5, "pt"),  # Add top margin for second axis
    plot.subtitle = ggtext::element_markdown(hjust = 0.5, margin = margin(b = 10))  # Enable HTML rendering
  )


cnl_ids <- readLines("/Users/kamal/Documents/NYU/raman-expression-classification-analysis/nlf-bg-sample-comparison/nlr_ids/cnl_ids.txt")
rnl_ids <- readLines("/Users/kamal/Documents/NYU/raman-expression-classification-analysis/nlf-bg-sample-comparison/nlr_ids/rnl_ids.txt")
tnl_ids <- readLines("/Users/kamal/Documents/NYU/raman-expression-classification-analysis/nlf-bg-sample-comparison/nlr_ids/tnl_ids.txt")

raman_long_cnl_rnl_tnl <- raman_long_full_nlrs %>% filter(gene %in% c(cnl_ids, rnl_ids, tnl_ids))
# this ends up being the same list as raman_long_full_nlrs with 181 unique genes -- need to filter uniqueness early for
# final plot


# now make csvs by subcat_Choghag and age for raman_long that are above and below median expression of all genes
for (age in levels(raman_long$age)) {
  for (subcat in unique(raman_long$subcat_Choghag)) {
    if (is.na(subcat)) next  # Skip NA subcategories
    above_median_df <- raman_long[raman_long$age == age & raman_long$subcat_Choghag == subcat & raman_long$expression > median(long_data$expression), c('gene', 'expression')]
    below_median_df <- raman_long[raman_long$age == age & raman_long$subcat_Choghag == subcat & raman_long$expression <= median(long_data$expression) & raman_long$expression > 0, c('gene', 'expression')]
    zero_df <- raman_long[raman_long$age == age & raman_long$subcat_Choghag == subcat & raman_long$expression == 0, c('gene', 'expression')]
    
    write.csv(above_median_df, file = sprintf('raman_expression_%s_%s_above_median.csv', gsub("all_", "", age), subcat), row.names = FALSE)
    write.csv(below_median_df, file = sprintf('raman_expression_%s_%s_below_median.csv', gsub("all_", "", age), subcat), row.names = FALSE)
    write.csv(zero_df, file = sprintf('raman_expression_%s_%s_zero.csv', gsub("all_", "", age), subcat), row.names = FALSE)
  }
}

nlr_gene_list_by_type <- raman_stuff[, c('geneID', 'subcat_Choghag')]
write.csv(nlr_gene_list_by_type, file = 'nlr_gene_list_by_type.csv', row.names = FALSE, quote = FALSE)
