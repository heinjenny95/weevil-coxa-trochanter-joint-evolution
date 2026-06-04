options(stringsAsFactors = FALSE)

root <- "<MANUSCRIPT_PROJECT_ROOT>"
in_csv <- file.path(root, "05_Results", "05_Coxa_Thickness", "coxa_combined_metrics.csv")
out_dir <- file.path(root, "05_Results", "05_Coxa_Thickness")

out_png <- file.path(out_dir, "coxa_thickness_2x2.png")
out_pdf <- file.path(out_dir, "coxa_thickness_2x2.pdf")
out_multipage_pdf <- file.path(out_dir, "plots_coxa_thickness_nohole.pdf")
stats_csv <- file.path(out_dir, "coxa_thickness_lm_stats.csv")

if (!requireNamespace("tidyverse", quietly = TRUE)) stop("Package tidyverse is required.")
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Package ggrepel is required.")
if (!requireNamespace("patchwork", quietly = TRUE)) stop("Package patchwork is required.")

library(tidyverse)
library(ggrepel)
library(patchwork)

df <- read.csv2(in_csv, stringsAsFactors = FALSE, check.names = FALSE)

needed <- c("specimen", "bbox_diag_um", "median_thickness_um", "first_slice", "last_slice", "mid_slice")
missing <- setdiff(needed, names(df))
if (length(missing) > 0) {
  stop("Missing required column(s): ", paste(missing, collapse = ", "))
}

df <- df %>%
  mutate(
    specimen = as.character(specimen),
    bbox_diag_um = as.numeric(bbox_diag_um),
    median_thickness_um = as.numeric(median_thickness_um),
    first_slice = as.numeric(first_slice),
    last_slice = as.numeric(last_slice),
    mid_slice = as.numeric(mid_slice),
    span_slices = last_slice - first_slice,
    mid_rel = (mid_slice - first_slice) / pmax(span_slices, 1)
  ) %>%
  filter(is.finite(bbox_diag_um), is.finite(median_thickness_um))

label_df <- bind_rows(
  df %>% arrange(desc(median_thickness_um)) %>% slice_head(n = 6),
  df %>% arrange(median_thickness_um) %>% slice_head(n = 6)
) %>%
  distinct(specimen, .keep_all = TRUE)

base_theme <- theme_classic(base_size = 11)

p1 <- ggplot(df, aes(x = bbox_diag_um, y = median_thickness_um)) +
  geom_point(alpha = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "Median cuticle thickness vs Coxa size",
    x = "Coxa size (bbox diagonal, um) [log10]",
    y = "Median cuticle thickness (um) [log10]"
  ) +
  base_theme

p2 <- ggplot(df, aes(x = median_thickness_um)) +
  geom_histogram(bins = 25) +
  labs(
    title = "Distribution of median cuticle thickness",
    x = "Median cuticle thickness (um)",
    y = "Count"
  ) +
  base_theme

p3 <- ggplot(df, aes(x = bbox_diag_um, y = median_thickness_um)) +
  geom_point(alpha = 0.7) +
  geom_point(data = label_df, size = 2) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = specimen),
    max.overlaps = Inf,
    size = 3
  ) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Thickness vs size (labels = extreme thickness values)",
    x = "Coxa size (bbox diagonal, um) [log10]",
    y = "Median cuticle thickness (um) [log10]"
  ) +
  base_theme

p4 <- ggplot(df, aes(x = span_slices, y = mid_rel)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_point(alpha = 0.8) +
  labs(
    title = "Mid-slice QC",
    x = "Coxa extent in Z (last_slice - first_slice)",
    y = "Relative mid position (0..1; ideally ~0.5)"
  ) +
  base_theme

combined_plot <- (p1 + p2) / (p3 + p4)

ggsave(out_png, combined_plot, width = 12, height = 10, dpi = 300)
ggsave(out_pdf, combined_plot, width = 12, height = 10, device = cairo_pdf)

pdf(out_multipage_pdf, width = 8.5, height = 6.5, useDingbats = FALSE)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

model <- lm(log10(median_thickness_um) ~ log10(bbox_diag_um), data = df)
model_sum <- summary(model)

coef_tab <- as.data.frame(model_sum$coefficients)
coef_tab$term <- rownames(coef_tab)
rownames(coef_tab) <- NULL

coef_tab <- coef_tab %>%
  select(term, everything()) %>%
  rename(
    estimate = Estimate,
    std_error = `Std. Error`,
    t_value = `t value`,
    p_value = `Pr(>|t|)`
  ) %>%
  mutate(type = "coefficient")

fstat <- unname(model_sum$fstatistic)
stats_tab <- data.frame(
  type = "model",
  term = "overall_model",
  estimate = NA_real_,
  std_error = NA_real_,
  t_value = NA_real_,
  p_value = pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE),
  r_squared = model_sum$r.squared,
  adj_r_squared = model_sum$adj.r.squared,
  f_statistic = fstat[1],
  df1 = fstat[2],
  df2 = fstat[3],
  residual_se = model_sum$sigma,
  n = nrow(df)
)

coef_tab <- coef_tab %>%
  mutate(
    r_squared = NA_real_,
    adj_r_squared = NA_real_,
    f_statistic = NA_real_,
    df1 = NA_real_,
    df2 = NA_real_,
    residual_se = NA_real_,
    n = NA_real_
  )

write.csv2(bind_rows(coef_tab, stats_tab), stats_csv, row.names = FALSE)

cat("Wrote:", out_png, "\n")
cat("Wrote:", out_pdf, "\n")
cat("Wrote:", out_multipage_pdf, "\n")
cat("Wrote:", stats_csv, "\n")
