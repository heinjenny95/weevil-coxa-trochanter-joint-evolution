# ============================================================
# Cluster analysis on morphospace PCA scores (PC-space)
# - k-means (k selection: silhouette + gap statistic)
# - hierarchical clustering (Ward.D2; k selection: silhouette)
# - model-based clustering (mclust; k selection: BIC)
#
# Inputs: German CSV with ; separator (read.csv2)
# Output: PDFs + CSVs in Results/Clustering/
# ============================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(cluster)      # silhouette
  library(factoextra)   # nice plotting helpers (kmeans/hclust/silhouette/gap)
  library(mclust)       # model-based clustering
})

# ------------------- SETTINGS -------------------
pca_path <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/PCA_scores_with_specimen_id.csv"

out_dir <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Results/Clustering"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

k_pcs <- 6                 # use PC1..PCk
k_range <- 2:10            # candidate k to evaluate
seed <- 1

# If TRUE: standardize PCs before clustering (recommended for kmeans/hclust)
# PCA scores often already comparable, but scaling is a safe default.
do_scale <- TRUE

# ------------------- READ DATA -------------------
df <- read.csv2(pca_path, stringsAsFactors = FALSE, check.names = FALSE)

# detect PC columns
pc_cols <- grep("^PC[0-9]+$", names(df), value = TRUE)
if (length(pc_cols) < 2) stop("Found fewer than 2 PC columns in PCA table.")
pcs_use <- pc_cols[seq_len(min(k_pcs, length(pc_cols)))]

# ensure numeric
for (cc in pcs_use) df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))

# require specimen_id for output labeling
if (!("specimen_id" %in% names(df))) stop("specimen_id column not found in PCA CSV.")

# drop rows with missing in selected PCs
df2 <- df %>%
  filter(if_all(all_of(pcs_use), ~ !is.na(.x))) %>%
  select(specimen_id, all_of(pcs_use))

cat("Using PCs:", paste(pcs_use, collapse = ", "), "\n")
cat("N specimens:", nrow(df2), "\n")

# matrix for clustering
X <- as.matrix(df2[, pcs_use, drop = FALSE])
if (do_scale) X <- scale(X)

# for plotting in 2D
PC1_name <- pcs_use[1]
PC2_name <- pcs_use[2]

# helper: save ggplot to pdf
save_pdf <- function(plot_obj, file, w = 7.5, h = 5.5) {
  ggsave(filename = file, plot = plot_obj, width = w, height = h, units = "in")
}

# ============================================================
# 1) k-means clustering
#   - choose k via average silhouette & gap statistic
# ============================================================

set.seed(seed)

# --- silhouette across k ---
sil_scores <- map_dfr(k_range, function(k) {
  km <- kmeans(X, centers = k, nstart = 50)
  sil <- silhouette(km$cluster, dist(X))
  tibble(k = k, avg_sil = mean(sil[, "sil_width"]))
})

k_best_sil <- sil_scores$k[which.max(sil_scores$avg_sil)]
cat("k-means: best k by silhouette =", k_best_sil, "\n")

p_sil <- ggplot(sil_scores, aes(x = k, y = avg_sil)) +
  geom_line() + geom_point() +
  theme_classic() +
  labs(title = "k-means: average silhouette vs k", x = "k", y = "Average silhouette")

save_pdf(p_sil, file.path(out_dir, "kmeans_silhouette_vs_k.pdf"))

# --- gap statistic (computationally heavier) ---
# factoextra's fviz_nbclust returns plot-ready object
set.seed(seed)
gap_plot <- fviz_nbclust(X, kmeans, method = "gap_stat", k.max = max(k_range), nstart = 50)
save_pdf(gap_plot, file.path(out_dir, "kmeans_gap_statistic.pdf"))

# Fit kmeans with k chosen by silhouette (you can swap to gap-based if you prefer)
set.seed(seed)
km_fit <- kmeans(X, centers = k_best_sil, nstart = 100)
df_km <- df2 %>% mutate(cluster_kmeans = factor(km_fit$cluster))

# Plot: morphospace PC1-PC2 colored by kmeans cluster
p_km <- ggplot(df_km, aes(x = .data[[PC1_name]], y = .data[[PC2_name]], color = cluster_kmeans)) +
  geom_point(size = 2, alpha = 0.9) +
  theme_classic() +
  labs(title = paste0("k-means clustering (k=", k_best_sil, ") on ", paste(pcs_use, collapse = ",")),
       x = PC1_name, y = PC2_name, color = "Cluster")

save_pdf(p_km, file.path(out_dir, paste0("kmeans_clusters_k", k_best_sil, "_PC1_PC2.pdf")))

# Optional silhouette plot for chosen k
sil_km <- silhouette(km_fit$cluster, dist(X))
p_km_silplot <- fviz_silhouette(sil_km) + theme_classic() +
  labs(title = paste0("k-means silhouette plot (k=", k_best_sil, ")"))
save_pdf(p_km_silplot, file.path(out_dir, paste0("kmeans_silhouette_plot_k", k_best_sil, ".pdf")))

# write assignments
write.csv2(df_km, file.path(out_dir, paste0("clusters_kmeans_k", k_best_sil, ".csv")), row.names = FALSE)


# ============================================================
# 2) Hierarchical clustering (Ward.D2)
#   - choose k via average silhouette across k
#   - outputs dendrogram + morphospace scatter
# ============================================================

dmat <- dist(X)
hc <- hclust(dmat, method = "ward.D2")

# silhouette across k for hclust cuts
sil_scores_hc <- map_dfr(k_range, function(k) {
  cl <- cutree(hc, k = k)
  sil <- silhouette(cl, dmat)
  tibble(k = k, avg_sil = mean(sil[, "sil_width"]))
})

k_best_hc <- sil_scores_hc$k[which.max(sil_scores_hc$avg_sil)]
cat("hclust: best k by silhouette =", k_best_hc, "\n")

p_sil_hc <- ggplot(sil_scores_hc, aes(x = k, y = avg_sil)) +
  geom_line() + geom_point() +
  theme_classic() +
  labs(title = "Hierarchical (Ward.D2): average silhouette vs k", x = "k", y = "Average silhouette")

save_pdf(p_sil_hc, file.path(out_dir, "hclust_silhouette_vs_k.pdf"))

# dendrogram (base plotting into PDF)
pdf(file.path(out_dir, "hclust_dendrogram.pdf"), width = 10, height = 6)
plot(hc, labels = FALSE, main = "Hierarchical clustering dendrogram (Ward.D2)", xlab = "", sub = "")
rect.hclust(hc, k = k_best_hc, border = 2:6)
dev.off()

# assignments + scatter
hc_cl <- cutree(hc, k = k_best_hc)
df_hc <- df2 %>% mutate(cluster_hclust = factor(hc_cl))

p_hc <- ggplot(df_hc, aes(x = .data[[PC1_name]], y = .data[[PC2_name]], color = cluster_hclust)) +
  geom_point(size = 2, alpha = 0.9) +
  theme_classic() +
  labs(title = paste0("Hierarchical clustering (Ward.D2, k=", k_best_hc, ") on ", paste(pcs_use, collapse = ",")),
       x = PC1_name, y = PC2_name, color = "Cluster")

save_pdf(p_hc, file.path(out_dir, paste0("hclust_clusters_k", k_best_hc, "_PC1_PC2.pdf")))

# optional silhouette plot for chosen k
sil_hc <- silhouette(hc_cl, dmat)
p_hc_silplot <- fviz_silhouette(sil_hc) + theme_classic() +
  labs(title = paste0("Hierarchical silhouette plot (k=", k_best_hc, ")"))
save_pdf(p_hc_silplot, file.path(out_dir, paste0("hclust_silhouette_plot_k", k_best_hc, ".pdf")))

write.csv2(df_hc, file.path(out_dir, paste0("clusters_hclust_k", k_best_hc, ".csv")), row.names = FALSE)


# ============================================================
# 3) Model-based clustering (mclust)
#   - choose model + G (k) via BIC
#   - outputs BIC plot + morphospace scatter
# ============================================================

# mclust expects unscaled data sometimes; but scaled is okay.
# We'll use X as currently prepared to be consistent.
mc <- Mclust(X, G = k_range)

# BIC plot
pdf(file.path(out_dir, "mclust_BIC.pdf"), width = 9, height = 6)
plot(mc, what = "BIC")
dev.off()

# selected number of clusters
k_best_mc <- mc$G
cat("mclust: selected G by BIC =", k_best_mc, " model =", mc$modelName, "\n")

df_mc <- df2 %>%
  mutate(cluster_mclust = factor(mc$classification))

p_mc <- ggplot(df_mc, aes(x = .data[[PC1_name]], y = .data[[PC2_name]], color = cluster_mclust)) +
  geom_point(size = 2, alpha = 0.9) +
  theme_classic() +
  labs(title = paste0("mclust (G=", k_best_mc, ", model=", mc$modelName, ") on ", paste(pcs_use, collapse = ",")),
       x = PC1_name, y = PC2_name, color = "Cluster")

save_pdf(p_mc, file.path(out_dir, paste0("mclust_clusters_G", k_best_mc, "_PC1_PC2.pdf")))

write.csv2(df_mc, file.path(out_dir, paste0("clusters_mclust_G", k_best_mc, ".csv")), row.names = FALSE)


# ============================================================
# 4) Combine cluster assignments into one file
# ============================================================

df_all <- df2 %>%
  left_join(df_km %>% select(specimen_id, cluster_kmeans), by = "specimen_id") %>%
  left_join(df_hc %>% select(specimen_id, cluster_hclust), by = "specimen_id") %>%
  left_join(df_mc %>% select(specimen_id, cluster_mclust), by = "specimen_id")

write.csv2(df_all, file.path(out_dir, "clusters_all_methods.csv"), row.names = FALSE)

cat("\nDONE. Outputs written to:\n", out_dir, "\n")



#################################################################
####################################################################
##################################################################


# ============================================================
# Combined PC1-PC2 cluster plots (k-means | hclust | mclust)
# Output: single 13 PDF
# ============================================================
# ============================================================
# Combined PC1-PC2 cluster plots (k-means | hclust | mclust)
# Output: single 13 PNG
# ============================================================

library(ggplot2)
library(patchwork)

# ------------------- SETTINGS -------------------
PC1_name <- "PC1"
PC2_name <- "PC2"

out_path <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Results/Clustering/combined_clusters_PC1_PC2.png"

# ------------------- AXIS LIMITS (shared) -------------------
xlims <- range(df_all[[PC1_name]], na.rm = TRUE)
ylims <- range(df_all[[PC2_name]], na.rm = TRUE)

# ------------------- k-means plot -------------------
p_km <- ggplot(
  df_all,
  aes(x = .data[[PC1_name]],
      y = .data[[PC2_name]],
      color = cluster_kmeans)
) +
  geom_point(size = 2, alpha = 0.9) +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  theme_classic() +
  labs(
    title = "k-means",
    x = PC1_name,
    y = PC2_name,
    color = "Cluster"
  )

# ------------------- hierarchical plot -------------------
p_hc <- ggplot(
  df_all,
  aes(x = .data[[PC1_name]],
      y = .data[[PC2_name]],
      color = cluster_hclust)
) +
  geom_point(size = 2, alpha = 0.9) +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  theme_classic() +
  labs(
    title = "Hierarchical (Ward.D2)",
    x = PC1_name,
    y = PC2_name,
    color = "Cluster"
  )

# ------------------- mclust plot -------------------
p_mc <- ggplot(
  df_all,
  aes(x = .data[[PC1_name]],
      y = .data[[PC2_name]],
      color = cluster_mclust)
) +
  geom_point(size = 2, alpha = 0.9) +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  theme_classic() +
  labs(
    title = "Model-based (mclust)",
    x = PC1_name,
    y = PC2_name,
    color = "Cluster"
  )

# ------------------- combine & save -------------------
p_combined <- p_km | p_hc | p_mc

ggsave(
  filename = out_path,
  plot = p_combined,
  width = 4200,      # pixels
  height = 1500,     # pixels
  dpi = 300,
  units = "px"
)

cat("Combined cluster PNG saved to:\n", out_path, "\n")




# ============================================================
# Hopkins statistic for morphospace clusterability
# (robust across hopkins package versions)
# ============================================================

library(tidyverse)
library(hopkins)

# ------------------- SETTINGS -------------------
pca_path <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/PCA_scores_with_specimen_id.csv"
k_pcs <- 6
n_boot <- 100
set.seed(1)

# ------------------- READ DATA -------------------
df <- read.csv2(pca_path, stringsAsFactors = FALSE)

pc_cols <- grep("^PC[0-9]+$", names(df), value = TRUE)
pcs_use <- pc_cols[seq_len(min(k_pcs, length(pc_cols)))]

df2 <- df %>%
  filter(if_all(all_of(pcs_use), ~ !is.na(.x))) %>%
  select(all_of(pcs_use))

X <- as.matrix(df2)
X <- scale(X)

cat("Using PCs:", paste(pcs_use, collapse = ", "), "\n")
cat("N specimens:", nrow(X), "\n")

# ------------------- helper: extract numeric Hopkins value -------------------
get_H <- function(x) {
  if (is.list(x)) {
    if (!is.null(x$H)) return(as.numeric(x$H))
    # fallback: if list but different field name, take first numeric
    nums <- unlist(x)
    return(as.numeric(nums[1]))
  } else {
    return(as.numeric(x))
  }
}

m_val <- floor(0.1 * nrow(X))

# ------------------- SINGLE HOPKINS TEST -------------------
H_single_raw <- hopkins::hopkins(X, m = m_val)
H_single <- get_H(H_single_raw)

cat("Hopkins statistic (single run):", round(H_single, 3), "\n")

# ------------------- BOOTSTRAPPED HOPKINS -------------------
H_vals <- replicate(n_boot, {
  get_H(hopkins::hopkins(X, m = m_val))
})

H_tbl <- tibble(Hopkins = H_vals)

H_summary <- H_tbl %>%
  summarise(
    mean = mean(Hopkins),
    sd   = sd(Hopkins),
    min  = min(Hopkins),
    max  = max(Hopkins)
  )

print(H_summary)

# ------------------- SAVE RESULTS -------------------
out_dir <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Results/Clustering"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv2(H_tbl, file.path(out_dir, "hopkins_bootstrap_values.csv"), row.names = FALSE)
write.csv2(H_summary, file.path(out_dir, "hopkins_summary.csv"), row.names = FALSE)

# ------------------- VISUALIZATION -------------------
png(
  file.path(out_dir, "hopkins_distribution.png"),
  width = 2000,
  height = 1200,
  res = 300
)
hist(
  H_vals,
  breaks = 20,
  col = "grey70",
  border = "white",
  main = "Hopkins statistic (bootstrap)",
  xlab = "Hopkins H"
)
abline(v = 0.5, lty = 2, col = "red")
dev.off()

cat("Hopkins test finished. Results written to:\n", out_dir, "\n")



####################Within Cluster Disparity###################################

# ============================================================
# 5) Within-cluster disparity across clustering methods
#    - computes per-specimen Euclidean distance to cluster centroid
#    - based on the SAME PC-space used for clustering (scaled X)
#    - outputs combined violin plot + CSV
# ============================================================

library(tidyverse)
library(patchwork)

# ------------------- sanity -------------------
if (!exists("X")) stop("Object 'X' not found. Run clustering part first.")
if (!exists("df_km")) stop("Object 'df_km' not found. Run k-means part first.")
if (!exists("df_hc")) stop("Object 'df_hc' not found. Run hclust part first.")
if (!exists("df_mc")) stop("Object 'df_mc' not found. Run mclust part first.")
if (!exists("out_dir")) stop("Object 'out_dir' not found.")

# ------------------- helper function -------------------
calc_dist_to_centroid <- function(X, clusters, specimen_ids, method_name) {
  cl <- as.factor(clusters)
  
  # centroid per cluster
  centroids <- sapply(levels(cl), function(g) {
    colMeans(X[cl == g, , drop = FALSE])
  })
  
  # if only one cluster level survives weirdly, keep matrix shape stable
  centroids <- as.matrix(centroids)
  if (ncol(centroids) == 1) {
    colnames(centroids) <- levels(cl)[1]
  }
  
  # distance of each specimen to own centroid
  dists <- numeric(nrow(X))
  for (i in seq_len(nrow(X))) {
    g <- as.character(cl[i])
    ctr <- centroids[, g]
    dists[i] <- sqrt(sum((X[i, ] - ctr)^2))
  }
  
  tibble(
    specimen_id = specimen_ids,
    method = method_name,
    cluster = cl,
    dist_to_centroid = dists
  )
}

# ------------------- compute distances -------------------
disp_km <- calc_dist_to_centroid(
  X = X,
  clusters = df_km$cluster_kmeans,
  specimen_ids = df_km$specimen_id,
  method_name = "k-means"
)

disp_hc <- calc_dist_to_centroid(
  X = X,
  clusters = df_hc$cluster_hclust,
  specimen_ids = df_hc$specimen_id,
  method_name = "Hierarchical"
)

disp_mc <- calc_dist_to_centroid(
  X = X,
  clusters = df_mc$cluster_mclust,
  specimen_ids = df_mc$specimen_id,
  method_name = "Model-based"
)

disp_all <- bind_rows(disp_km, disp_hc, disp_mc) %>%
  mutate(
    method = factor(method, levels = c("k-means", "Hierarchical", "Model-based"))
  )

# ------------------- summary table -------------------
disp_summary <- disp_all %>%
  group_by(method, cluster) %>%
  summarise(
    n = n(),
    mean_dist = mean(dist_to_centroid),
    median_dist = median(dist_to_centroid),
    sd_dist = sd(dist_to_centroid),
    .groups = "drop"
  )

write.csv2(disp_all, file.path(out_dir, "within_cluster_disparity_all_specimens.csv"), row.names = FALSE)
write.csv2(disp_summary, file.path(out_dir, "within_cluster_disparity_summary.csv"), row.names = FALSE)

# ------------------- plot prep: cluster means for overlay -------------------
mean_points <- disp_all %>%
  group_by(method, cluster) %>%
  summarise(
    mean_dist = mean(dist_to_centroid),
    n = n(),
    .groups = "drop"
  )

# ------------------- combined violin plot (colored) -------------------

p_disp <- ggplot(disp_all, aes(x = method, y = dist_to_centroid, fill = method)) +
  geom_violin(
    color = "black",
    alpha = 0.8,
    width = 0.9,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = "white",
    color = "black"
  ) +
  geom_jitter(
    width = 0.08,
    height = 0,
    alpha = 0.4,
    size = 1.2
  ) +
  geom_point(
    data = mean_points,
    aes(x = method, y = mean_dist, size = n),
    inherit.aes = FALSE,
    shape = 21,
    fill = "black",
    color = "white",
    stroke = 0.3
  ) +
  scale_fill_manual(
    values = c(
      "k-means" = "#4C72B0",        # blau
      "Hierarchical" = "#DD8452",   # orange
      "Model-based" = "#55A868"     # grn
    )
  ) +
  scale_size_continuous(name = "Cluster size (n)") +
  labs(
    title = "Within-cluster disparity across clustering methods",
    x = NULL,
    y = "Distance to cluster centroid"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.position = "none"  # <- wichtig, sonst doppelt
  )

print(p_disp)
