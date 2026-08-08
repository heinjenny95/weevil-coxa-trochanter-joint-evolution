#!/usr/bin/env Rscript

# Build the transparent specimen-to-phylogeny mapping used by the comparative
# analyses. The source tree does not contain most sampled genera, so the
# resulting 15 labels are taxonomic proxy tips rather than specimen genera.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: build_proxy_mapping.R <specimen_key.csv> <output.csv>")
}

input_path <- args[[1]]
output_path <- args[[2]]

read_csv_robust <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (!length(lines)) stop("Empty input: ", path)
  lines[[1]] <- sub("^\ufeff", "", lines[[1]])
  if (grepl("^sep\\s*=", lines[[1]], ignore.case = TRUE)) {
    sep <- substr(sub("^sep\\s*=\\s*", "", lines[[1]], ignore.case = TRUE), 1, 1)
    lines <- lines[-1]
  } else {
    header <- lines[[1]]
    sep <- if (lengths(regmatches(header, gregexpr(";", header, fixed = TRUE))) >=
               lengths(regmatches(header, gregexpr(",", header, fixed = TRUE)))) ";" else ","
  }
  textConnectionInput <- textConnection(lines)
  on.exit(close(textConnectionInput), add = TRUE)
  read.table(
    textConnectionInput,
    header = TRUE,
    sep = sep,
    quote = "\"",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

write_csv_clean <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(
    x,
    file = path,
    sep = ";",
    row.names = FALSE,
    col.names = TRUE,
    quote = TRUE,
    na = "",
    fileEncoding = "UTF-8"
  )
}

specimens <- read_csv_robust(input_path)
required <- c("specimen_key", "taxon_binomial", "specimen_id", "tree_tip", "Family")
missing <- setdiff(required, names(specimens))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

specimens$tree_tip[specimens$tree_tip == "Neydus"] <- "Nedyus"
specimens$sampled_genus <- sub(" .*", "", trimws(specimens$taxon_binomial))
specimens$proxy_tip <- trimws(specimens$tree_tip)
specimens$exact_genus_match <- tolower(specimens$sampled_genus) == tolower(specimens$proxy_tip)
tip_counts <- table(specimens$proxy_tip)
specimens$proxy_tip_n_specimens <- as.integer(tip_counts[specimens$proxy_tip])
specimens$assignment_type <- ifelse(
  specimens$exact_genus_match,
  "exact sampled-genus match",
  "broader taxonomic proxy"
)

mapping <- specimens[, c(
  "specimen_key",
  "specimen_id",
  "taxon_binomial",
  "sampled_genus",
  "Family",
  "proxy_tip",
  "proxy_tip_n_specimens",
  "exact_genus_match",
  "assignment_type"
)]
names(mapping)[names(mapping) == "Family"] <- "family"

mapping <- mapping[order(mapping$family, mapping$proxy_tip, mapping$taxon_binomial), ]
write_csv_clean(mapping, output_path)

cat("Specimens:", nrow(mapping), "\n")
cat("Proxy tips:", length(unique(mapping$proxy_tip)), "\n")
cat("Exact sampled-genus matches:", sum(mapping$exact_genus_match), "\n")
cat("Broader taxonomic proxy assignments:", sum(!mapping$exact_genus_match), "\n")
