# ============================================================
# Extract 1 extreme specimen per PC in BOTH directions:
#   - single minimum AND single maximum per PC (PC1-PC5)
# Output:
#   1) combined CSV: PC1_PC5_extremes_minmax_single.csv
#   2) one CSV per PC: PC1_extremes_minmax_single.csv, ...
# ============================================================

library(dplyr)

# ------------------- SETTINGS -------------------
pca_path <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/PCA_scores_with_specimen_id.csv"
out_dir  <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Results/PCA/PC_extremes_single_minmax"

pcs_to_use <- paste0("PC", 1:5)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------- READ DATA -------------------
df <- read.csv2(pca_path, stringsAsFactors = FALSE, check.names = FALSE)

stopifnot("specimen_id" %in% names(df))
stopifnot(all(pcs_to_use %in% names(df)))

df <- df %>% dplyr::select(specimen_id, dplyr::all_of(pcs_to_use))

# ------------------- FUNCTION: get min + max (single each) -------------------
get_minmax_one_pc <- function(dat, pc) {
  dat_ok <- dat %>% dplyr::filter(!is.na(.data[[pc]]))
  if (nrow(dat_ok) == 0) stop("No non-NA values for ", pc)
  
  min_row <- dat_ok %>%
    dplyr::arrange(.data[[pc]]) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(PC = pc, extreme = "min", PC_score = .data[[pc]]) %>%
    dplyr::select(PC, extreme, specimen_id, PC_score)
  
  max_row <- dat_ok %>%
    dplyr::arrange(dplyr::desc(.data[[pc]])) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(PC = pc, extreme = "max", PC_score = .data[[pc]]) %>%
    dplyr::select(PC, extreme, specimen_id, PC_score)
  
  dplyr::bind_rows(min_row, max_row)
}

# ------------------- RUN -------------------
ext_all <- dplyr::bind_rows(lapply(pcs_to_use, get_minmax_one_pc, dat = df))

# ------------------- WRITE OUTPUTS -------------------
write.csv2(
  ext_all,
  file.path(out_dir, "PC1_PC5_extremes_minmax_single.csv"),
  row.names = FALSE
)

for (pc in pcs_to_use) {
  one <- ext_all %>% dplyr::filter(PC == pc) %>% dplyr::arrange(extreme)
  write.csv2(
    one,
    file.path(out_dir, paste0(pc, "_extremes_minmax_single.csv")),
    row.names = FALSE
  )
}

# ------------------- CONSOLE SUMMARY -------------------
cat("\nDone. Files written to:\n", out_dir, "\n\n")
cat("Extremes per PC:\n")
print(ext_all %>% dplyr::arrange(PC, extreme))
