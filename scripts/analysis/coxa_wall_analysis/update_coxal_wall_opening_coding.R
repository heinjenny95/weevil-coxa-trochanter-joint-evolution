root <- "<MANUSCRIPT_PROJECT_ROOT>"

candidate_files <- c(
  "03_Data_and_Inputs/Analysis_Input/specimen_key.csv",
  "03_Data_and_Inputs/Analysis_Input/specimen_key_gsheets.csv",
  "03_Data_and_Inputs/Analysis_Input/specimen_key_with_centroid_size.csv",
  "03_Data_and_Inputs/Analysis_Input/specimen_key_with_centroid_size_gsheets.csv",
  "05_Results/03_Allometry/Allometry/specimen_key_with_centroid_size.csv",
  "05_Results/03_Allometry/Allometry/specimen_key_with_centroid_size_gsheets.csv",
  "05_Results/03_Allometry/Allometry/PCA_scores_with_specimen_id_with_centroid_size.csv",
  "05_Results/03_Allometry/Allometry/PCA_scores_with_specimen_id_with_centroid_size_gsheets.csv",
  "05_Results/03_Allometry/Allometry/allometry_merged_table.csv",
  "05_Results/03_Allometry/Allometry/allometry_merged_table_gsheets.csv",
  "05_Results/03_Allometry/Allometry_Phylogenetic/pgls_input_specimen_level.csv",
  "05_Results/03_Allometry/Allometry_Phylogenetic/pgls_input_specimen_level_gsheets.csv",
  "05_Results/03_Allometry/Allometry_Phylogenetic/pgls_input_aggregated_to_tree_tip.csv",
  "05_Results/03_Allometry/Allometry_Phylogenetic/pgls_input_aggregated_to_tree_tip_gsheets.csv",
  "05_Results/07_PCM/10_Logs/tip_level_dataset_used_for_PCM.csv",
  "05_Results/07_PCM/10_Logs/tip_level_dataset_used_for_PCM_gsheets.csv"
)

paths <- file.path(root, candidate_files)
paths <- paths[file.exists(paths)]

backup_dir <- file.path(
  root,
  "04_Analysis_Code",
  "coxa_wall_analysis",
  paste0("backup_before_anthribidae_coxal_wall_false_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)

detect_sep <- function(path) {
  first <- readLines(path, n = 1, warn = FALSE, encoding = "UTF-8")
  if (length(first) == 0) return(";")
  semis <- gregexpr(";", first, fixed = TRUE)[[1]]
  commas <- gregexpr(",", first, fixed = TRUE)[[1]]
  n_semis <- sum(semis > 0)
  n_commas <- sum(commas > 0)
  if (n_semis >= n_commas) {
    ";"
  } else {
    ","
  }
}

read_table_chars <- function(path, sep) {
  read.table(
    path,
    sep = sep,
    header = TRUE,
    quote = "\"",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = "character",
    comment.char = "",
    na.strings = character(),
    fileEncoding = "UTF-8-BOM"
  )
}

family_mask <- function(df) {
  fam_cols <- intersect(c("Family", "Family.x", "Family.y"), names(df))
  if (!length(fam_cols)) return(rep(FALSE, nrow(df)))
  Reduce(`|`, lapply(fam_cols, function(nm) trimws(df[[nm]]) == "Anthribidae"))
}

changed <- data.frame(
  file = character(),
  rows_changed = integer(),
  stringsAsFactors = FALSE
)

for (path in paths) {
  sep <- detect_sep(path)
  df <- read_table_chars(path, sep)
  value_cols <- intersect(c("Coxal wall hole", "coxal_wall_opening", "coxal_wall_hole"), names(df))
  if (!length(value_cols)) next

  idx <- family_mask(df)
  if (!any(idx)) next

  old <- df[idx, value_cols[1], drop = TRUE]
  needs_change <- idx
  needs_change[idx] <- trimws(old) != "FALSE"
  if (!any(needs_change)) next

  file.copy(path, file.path(backup_dir, basename(path)), overwrite = TRUE)
  df[needs_change, value_cols[1]] <- "FALSE"
  write.table(
    df,
    file = path,
    sep = sep,
    quote = TRUE,
    row.names = FALSE,
    col.names = TRUE,
    na = "",
    fileEncoding = "UTF-8"
  )

  changed <- rbind(changed, data.frame(
    file = path,
    rows_changed = sum(needs_change),
    stringsAsFactors = FALSE
  ))
}

log_path <- file.path(backup_dir, "changed_files.csv")
write.csv(changed, log_path, row.names = FALSE)

message("Anthribidae Coxal wall hole correction complete.")
message("Backup/log folder: ", backup_dir)
print(changed)
