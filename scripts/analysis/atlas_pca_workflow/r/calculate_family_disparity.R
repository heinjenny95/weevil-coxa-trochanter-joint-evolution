# ============================================================
# Disparity (within-family) - BARPLOT with bootstrap CI
# - Uses PCA scores (PC1..PCk)
# - Joins Family from specimen_key
# - Keeps only families with n > 2
# - Disparity = mean squared Euclidean distance to family centroid
# - Outputs:
#     1) disparity_by_family_k{k}.csv
#     2) disparity_by_family_k{k}.pdf (barplot)
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
})

# ------------------- SETTINGS -------------------
pca_path <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/PCA_scores_with_specimen_id.csv"
key_path <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/specimen_key_with_centroid_size.csv"

k_pcs    <- 5
B        <- 5000
ci_level <- 0.95
min_n    <- 3   # n > 2

# Output
out_dir <- dirname(pca_path)

# Optional: family colors (set to NULL to use ggplot default palette)
# If you want exact matching colors, fill this named vector with your real hex codes.
fam_cols <- c(
  "Anthribidae"    = "#F8766D",
  "Attelabidae"    = "#C49A00",
  "Belidae"        = "#7CAE00",
  "Brentidae"      = "#00BA38",
  "Caridae"        = "#00BFC4",
  "Curculionidae"  = "#C77CFF",
  "Nemonychidae"   = "#F564E3"
)

# ------------------- HELPERS -------------------
disparity_msd <- function(X) {
  cen <- colMeans(X)
  d2  <- rowSums((X - matrix(cen, nrow(X), ncol(X), byrow = TRUE))^2)
  mean(d2)
}

boot_ci <- function(X, B = 2000, ci_level = 0.95) {
  n <- nrow(X)
  stats <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    stats[b] <- disparity_msd(X[idx, , drop = FALSE])
  }
  alpha <- (1 - ci_level) / 2
  qs <- quantile(stats, probs = c(alpha, 1 - alpha), names = FALSE)
  setNames(qs, c("ci_low", "ci_high"))
}

# ------------------- READ + JOIN -------------------
pca <- read.csv2(pca_path, stringsAsFactors = FALSE, check.names = FALSE)
key <- read.csv2(key_path, stringsAsFactors = FALSE, check.names = FALSE)

stopifnot("specimen_id" %in% names(pca))
stopifnot("specimen_id" %in% names(key))
stopifnot("Family" %in% names(key))

df <- pca %>%
  left_join(key %>% select(specimen_id, Family), by = "specimen_id")

# detect PC columns
pc_cols <- grep("^PC[0-9]+$", names(df), value = TRUE)
if (length(pc_cols) < 2) stop("Found fewer than 2 PC columns in PCA table.")

# sort PCs numerically (PC1, PC2, PC10 safe)
pc_num <- as.integer(sub("^PC", "", pc_cols))
pc_cols <- pc_cols[order(pc_num)]

pcs_use <- pc_cols[seq_len(min(k_pcs, length(pc_cols)))]

# coerce PCs numeric
for (cc in pcs_use) df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))

# drop incomplete rows
df <- df %>%
  filter(!is.na(Family)) %>%
  filter(if_all(all_of(pcs_use), ~ !is.na(.x))) %>%
  mutate(Family = as.character(Family))

# keep only families with n > 2
fam_counts <- df %>%
  group_by(Family) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))

keep_fams  <- fam_counts %>% filter(n >= min_n) %>% pull(Family)

df2 <- df %>% filter(Family %in% keep_fams)

cat("Families retained (n >=", min_n, "):\n")
print(fam_counts %>% filter(Family %in% keep_fams) %>% arrange(desc(n)))

# ------------------- DISPARITY + BOOTSTRAP (efficient) -------------------
set.seed(1)

disp_out <- df2 %>%
  group_by(Family) %>%
  group_modify(~{
    X <- as.matrix(.x[, pcs_use, drop = FALSE])
    disp <- disparity_msd(X)
    ci <- boot_ci(X, B = B, ci_level = ci_level)
    tibble(
      n = nrow(X),
      disparity = disp,
      ci_low = unname(ci["ci_low"]),
      ci_high = unname(ci["ci_high"])
    )
  }) %>%
  ungroup() %>%
  arrange(desc(disparity)) %>%
  mutate(Family = factor(Family, levels = Family))

# ------------------- WRITE CSV -------------------
out_csv <- file.path(out_dir, paste0("disparity_by_family_k", k_pcs, "_nGT2.csv"))
write.csv2(disp_out, out_csv, row.names = FALSE)
cat("Wrote:", out_csv, "\n")

# ------------------- BARPLOT -------------------
# Warn if color map incomplete (only relevant if fam_cols is not NULL)
if (!is.null(fam_cols)) {
  missing_cols <- setdiff(levels(disp_out$Family), names(fam_cols))
  if (length(missing_cols) > 0) {
    message(" Missing colors for families: ", paste(missing_cols, collapse = ", "))
    message("Add them to fam_cols if you want exact matching colors.")
  }
}

p <- ggplot(disp_out, aes(x = Family, y = disparity, fill = Family)) +
  geom_col(width = 0.75, alpha = 0.95) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, linewidth = 0.5) +
  geom_text(aes(label = paste0("n=", n)), vjust = -0.6, size = 3) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1)
  ) +
  labs(
    title = paste0("Within-family disparity (PC1..PC", k_pcs, "), families with n > 2"),
    x = "Family",
    y = "Disparity (mean squared distance to family centroid)"
  ) +
  expand_limits(y = max(disp_out$ci_high, na.rm = TRUE) * 1.12)

if (!is.null(fam_cols)) {
  p <- p + scale_fill_manual(values = fam_cols, drop = FALSE)
}

# Save as PDF (publication-friendly)
out_pdf <- file.path(out_dir, paste0("disparity_by_family_barplot_k", k_pcs, "_nGT2.pdf"))
pdf(out_pdf, width = 10, height = 6)
print(p)
dev.off()
cat("Wrote:", out_pdf, "\n")

# Optional: also save as PNG
out_png <- file.path(out_dir, paste0("disparity_by_family_barplot_k", k_pcs, "_nGT2.png"))
ggsave(out_png, plot = p, width = 10, height = 6, dpi = 300)
cat("Wrote:", out_png, "\n")

cat("DONE.\n")
