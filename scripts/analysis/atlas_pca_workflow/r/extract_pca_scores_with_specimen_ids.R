# ============================================================
# Atlas_Momentas.txt + data_set.xml -> PCA_scores_with_specimen_id.csv
# FIXED reshape order: Dim -> Ncp -> Nsubj (Dim is fastest index)
# - Specimen IDs from XML (order-safe)
# - Exports ALL PCs (no limit)
# - German CSV (write.csv2): ; separator + , decimal
# - Exports PCA variance table as CSV
# - Writes PCA PDF incl. sanity plots + scree plot (eigenvalues in %)
# ============================================================

rm(list = ls())

# ------------------- INPUTS AND OUTPUTS -------------------
# Command-line use:
# Rscript extract_pca_scores_with_specimen_ids.R \
#   Atlas_Momentas.txt data_set.xml output_directory
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 3) {
  momenta_file <- args[1]
  xml_file <- args[2]
  out_dir <- args[3]
} else if (length(args) == 0) {
  momenta_file <- "<KERNEL_WIDTH_PROJECT_ROOT>/data/atlasing/Atlas_Momentas.txt"
  xml_file <- "<KERNEL_WIDTH_PROJECT_ROOT>/data/atlasing/data_set.xml"
  out_dir <- "<KERNEL_WIDTH_PROJECT_ROOT>/data/atlasing"
} else {
  stop(
    "Expected either no arguments or: Atlas_Momentas.txt data_set.xml output_directory",
    call. = FALSE
  )
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(out_dir, "PCA_scores_with_specimen_id.csv")
out_var_csv <- file.path(out_dir, "PCA_variance.csv")
out_pdf <- file.path(out_dir, "PCA_simple_plots.pdf")

# ------------------- PACKAGES -------------------
if (!requireNamespace("xml2", quietly = TRUE)) {
  stop("Package 'xml2' is required. Install via: install.packages('xml2')")
}
library(xml2)

# ------------------- HELPERS -------------------
stop_with_hint <- function(msg) stop(paste0("\n ", msg, "\n"), call. = FALSE)

clean_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("[\u00A0]", " ", x)
  x <- gsub("\\s+", " ", x)
  x
}

strip_vtk <- function(x) sub("\\.[Vv][Tt][Kk]$", "", x)

# ------------------- READ specimen IDs from XML (ORDER-SAFE) -------------------
doc <- read_xml(xml_file)
subj_nodes <- xml_find_all(doc, ".//subject")
if (length(subj_nodes) == 0) stop_with_hint("Keine <subject> Nodes in data_set.xml gefunden.")

subj_ids <- clean_id(xml_attr(subj_nodes, "id"))

# fallback if @id missing
if (any(is.na(subj_ids)) || any(subj_ids == "")) {
  fn_nodes <- xml_find_all(doc, ".//subject//filename")
  fn_text  <- clean_id(xml_text(fn_nodes))
  if (length(fn_text) != length(subj_nodes)) {
    stop_with_hint("Konnte weder subject @id noch passende <filename>-Liste extrahieren.")
  }
  subj_ids <- fn_text
}

specimen_ids <- strip_vtk(subj_ids)

# unique check
dups <- specimen_ids[duplicated(specimen_ids)]
if (length(dups) > 0) {
  stop_with_hint(paste0("Doppelte specimen_ids gefunden: ", paste(unique(dups), collapse = ", ")))
}

# ------------------- READ Atlas_Momentas -------------------
lines <- readLines(momenta_file, warn = FALSE)
lines <- lines[nchar(trimws(lines)) > 0]

hdr <- scan(text = lines[1], quiet = TRUE)
if (length(hdr) != 3) stop_with_hint("Header muss genau 3 Werte haben: Nsubj Ncp Dim")

Nsubj <- as.integer(hdr[1])
Ncp   <- as.integer(hdr[2])
Dim   <- as.integer(hdr[3])

message("Parsed header: Nsubj=", Nsubj, " Ncp=", Ncp, " Dim=", Dim)

if (length(specimen_ids) != Nsubj) {
  stop_with_hint(paste0(
    "Mismatch: subjects in XML = ", length(specimen_ids),
    " aber Nsubj in Atlas_Momentas = ", Nsubj
  ))
}

num_data <- scan(text = lines[-1], quiet = TRUE)

expected_len <- Nsubj * Ncp * Dim
if (length(num_data) != expected_len) {
  stop_with_hint(paste0(
    "Datenlnge passt nicht. Erwartet: ", expected_len,
    " Werte, gefunden: ", length(num_data)
  ))
}

# ------------------- CRITICAL FIX: RESHAPE ORDER -------------------
# File order: subject -> CP -> Dim (Dim fastest)
# In R: first dimension fastest => (Dim, Ncp, Nsubj)
momenta_arr <- array(num_data, dim = c(Dim, Ncp, Nsubj))

# Flatten per subject
momenta_mat <- matrix(NA_real_, nrow = Nsubj, ncol = Dim * Ncp)
for (i in seq_len(Nsubj)) {
  momenta_mat[i, ] <- as.vector(momenta_arr[ , , i])
}

# ------------------- PCA (ALL PCs) -------------------
pca <- prcomp(momenta_mat, center = TRUE, scale. = FALSE)

scores <- as.data.frame(pca$x)          # all PCs (max Nsubj-1)
colnames(scores) <- paste0("PC", seq_len(ncol(scores)))

out <- cbind(specimen_id = specimen_ids, scores)

# ------------------- WRITE German CSV (scores) -------------------
write.csv2(out, out_csv, row.names = FALSE, fileEncoding = "UTF-8")
message(" Wrote CSV: ", out_csv)
message(" Rows: ", nrow(out), " | PCs exported: ", ncol(out) - 1)

# ------------------- PCA VARIANCE TABLE (Eigenvalues + %) -------------------
eigenvalues <- pca$sdev^2
var_expl    <- eigenvalues / sum(eigenvalues) * 100
cum_var     <- cumsum(var_expl)

var_tbl <- data.frame(
  PC              = paste0("PC", seq_along(eigenvalues)),
  Eigenvalue      = eigenvalues,
  Varianz_prozent = var_expl,
  Kumulativ_prozent = cum_var,
  stringsAsFactors = FALSE
)

# German CSV for variance table
write.csv2(var_tbl, out_var_csv, row.names = FALSE, fileEncoding = "UTF-8")
message(" Wrote variance CSV: ", out_var_csv)

# ------------------- PCA PLOTS (SANITY CHECK + SCREE) -------------------
pdf(out_pdf, width = 7, height = 7)

# 1) PC1 vs PC2
plot(
  scores$PC1, scores$PC2,
  xlab = sprintf("PC1 (%.1f%%)", var_expl[1]),
  ylab = sprintf("PC2 (%.1f%%)", var_expl[2]),
  pch  = 19, cex = 0.8,
  main = "PCA: PC1 vs PC2"
)
text(scores$PC1, scores$PC2, labels = specimen_ids, pos = 3, cex = 0.5)

# 2) PC2 vs PC3 (only if exists)
if ("PC3" %in% names(scores)) {
  plot(
    scores$PC2, scores$PC3,
    xlab = sprintf("PC2 (%.1f%%)", var_expl[2]),
    ylab = sprintf("PC3 (%.1f%%)", var_expl[3]),
    pch  = 19, cex = 0.8,
    main = "PCA: PC2 vs PC3"
  )
  text(scores$PC2, scores$PC3, labels = specimen_ids, pos = 3, cex = 0.5)
}

# 3) Scree plot: eigenvalues as % variance explained
# (Barplot, weil du explizit "% von den PCs" willst)
k <- length(var_expl)
barplot(
  height = var_expl,
  names.arg = paste0("PC", seq_len(k)),
  las = 2,
  ylab = "Varianz erklrt (%)",
  xlab = "Principal Components",
  main = "Scree Plot: Eigenvalues (Varianz in %)"
)
abline(h = 0)

dev.off()
message(" Wrote PDF: ", out_pdf)

# ------------------- QUICK SELF-CHECK OUTPUT -------------------
first_triplet_file <- num_data[1:3]
first_triplet_arr  <- momenta_arr[, 1, 1]  # Dim=1..3, CP1, Subject1

message("Sanity check (should match closely):")
message("  File first 3 values:  ", paste(signif(first_triplet_file, 6), collapse = "  "))
message("  Array [Dim,CP1,S1]:   ", paste(signif(first_triplet_arr,  6), collapse = "  "))
