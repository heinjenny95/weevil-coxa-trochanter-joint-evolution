# ============================================================
# Eigenvalue / variance explained plot (5% cutoff)
# Robust: works whether variance is stored as 0-1 or 0-100
# ============================================================

library(tidyverse)

# ------------------- SETTINGS -------------------
var_path <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/PCA_variance.csv"

out_path <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Results/PCA/pca_variance_explained.png"
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

cutoff_pct <- 5  # 5% cutoff (in percent units)

# ------------------- READ DATA -------------------
var_df <- read.csv2(var_path, stringsAsFactors = FALSE)

# expected columns from your export: PC, Varianz_prozent (and optionally Kumulativ_prozent, Eigenvalue)
needed <- c("PC", "Varianz_prozent")
if (!all(needed %in% names(var_df))) {
  stop("pca_variance.csv must contain columns: ", paste(needed, collapse = ", "))
}

# ------------------- CLEAN + AUTO-DETECT SCALE -------------------
var_df <- var_df %>%
  mutate(
    PC_num = readr::parse_number(PC)
  ) %>%
  arrange(PC_num)

# Auto-detect if Varianz_prozent is in [0,1] or [0,100]
# Heuristic: if max <= 1.0001 -> treat as fraction; else percent
max_v <- max(var_df$Varianz_prozent, na.rm = TRUE)
is_fraction <- max_v <= 1.0001

var_df <- var_df %>%
  mutate(
    variance_pct = if (is_fraction) Varianz_prozent * 100 else Varianz_prozent,
    above_cutoff = variance_pct >= cutoff_pct,
    PC = factor(PC, levels = PC)   # keep order after sorting
  )

# ------------------- PLOT -------------------
# Show every 5th label to avoid clutter
break_idx <- seq(1, nrow(var_df), by = 5)
break_labs <- levels(var_df$PC)[break_idx]

p <- ggplot(var_df, aes(x = PC, y = variance_pct, fill = above_cutoff)) +
  geom_col(color = "black", width = 0.8) +
  geom_hline(yintercept = cutoff_pct, linetype = "dashed", color = "red") +
  scale_x_discrete(breaks = break_labs) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Variance explained by principal components",
    x = "Principal component",
    y = "Variance explained (%)"
  )

# ------------------- SAVE -------------------
ggsave(
  filename = out_path,
  plot = p,
  width = 10,
  height = 5,
  dpi = 300
)

cat("Eigenvalue/variance plot saved to:\n", out_path, "\n")

# ------------------- PCs >= cutoff -------------------
pcs_over <- var_df %>%
  filter(variance_pct >= cutoff_pct) %>%
  transmute(PC, variance_pct = round(variance_pct, 3))

cat("\nPCs explaining >= ", cutoff_pct, "% variance:\n", sep = "")
print(pcs_over)
