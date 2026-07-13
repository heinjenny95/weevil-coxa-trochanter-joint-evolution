# ============================================================
# Assign joint_type + screw_state per specimen_id (DE-CSV export)
# - Reads: joint_trait_table.xlsx
# - Normalizes traits (Schraube: TRUE/FALSE/X; others TRUE/FALSE)
# - Assigns:
#     joint_type (4 mechanical categories; ASCII-safe labels)
#     screw_state (clear / ambiguous)
#     joint_type_strict (NA if Schraube == X, for strict analyses)
# - Exports as German CSV (;) with Windows-1252 encoding (Excel-safe)
# - Optionally merges into specimen_key.csv and exports updated key (DE CSV)
#
# Columns expected in Excel:
# Family, Subfamily, Genus, Species, Tree_tip, specimen_id, Group,
# Schraube, Coxal wall hole, Coxal Socket, Windung Coxa
#
# Works with older dplyr (no count(name=...), no count(sort=...))
# ============================================================

rm(list = ls())

# -------------------- SETTINGS --------------------
in_xlsx <- "<BEETLE_JOINTS_ROOT>/joint_trait_table.xlsx"

# OPTIONAL: merge into specimen_key.csv if it exists
specimen_key_path <- "<MANUSCRIPT_PROJECT_ROOT>/analysis_data/Input/specimen_key.csv"
do_merge_specimen_key <- file.exists(specimen_key_path)

out_dir <- file.path(dirname(in_xlsx), "joint_type_export")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# German CSV outputs
out_map_csv <- file.path(out_dir, "specimen_joint_types_DE.csv")
out_key_csv <- file.path(out_dir, "specimen_key_WITH_joint_types_DE.csv")

# Encoding for Excel on Windows
CSV_ENCODING <- "Windows-1252"

# -------------------- PACKAGES --------------------
pkgs <- c("readxl", "dplyr", "tidyr", "stringr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# -------------------- HELPERS --------------------
clean_str <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("[\u00A0]", " ") %>%  # NBSP -> space
    stringr::str_squish() %>%
    stringr::str_trim()
}

norm_tristate <- function(x) {
  y <- clean_str(x)
  y_low <- tolower(y)
  y_low[y_low %in% c("", "na", "n/a", "null", "none")] <- NA_character_
  
  dplyr::case_when(
    is.na(y_low) ~ NA_character_,
    y_low %in% c("true","t","1","yes","y","ja","wahr") ~ "TRUE",
    y_low %in% c("false","f","0","no","n","nein","falsch") ~ "FALSE",
    y_low %in% c("x","unklar","ambiguous","intermediate","maybe") ~ "X",
    TRUE ~ y  # unexpected -> keep for warning
  )
}

norm_bi <- function(x) {
  y <- clean_str(x)
  y_low <- tolower(y)
  y_low[y_low %in% c("", "na", "n/a", "null", "none")] <- NA_character_
  
  dplyr::case_when(
    is.na(y_low) ~ NA_character_,
    y_low %in% c("true","t","1","yes","y","ja","wahr") ~ "TRUE",
    y_low %in% c("false","f","0","no","n","nein","falsch") ~ "FALSE",
    TRUE ~ y
  )
}

# old-dplyr-safe count (no name= support needed)
count_old <- function(data, ..., out_name = "n") {
  data %>%
    dplyr::group_by(...) %>%
    dplyr::summarise(n_tmp = dplyr::n(), .groups = "drop") %>%
    dplyr::rename(!!out_name := n_tmp)
}

# Read German/English CSV robustly (specimen_key can be ; or ,)
read_csv_any <- function(path) {
  # try German ; first
  x <- tryCatch(read.csv2(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(x) && ncol(x) > 1) return(x)
  
  # fallback to comma
  x2 <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(x2)) return(x2)
  
  stop("Could not read CSV: ", path)
}

# -------------------- READ EXCEL --------------------
df_raw <- readxl::read_excel(in_xlsx)
if (nrow(df_raw) == 0) stop("Excel loaded but has 0 rows: ", in_xlsx)

names(df_raw) <- clean_str(names(df_raw))

required <- c(
  "Family","Subfamily","Genus","Species","Tree_tip","specimen_id","Group",
  "Schraube","Coxal wall hole","Coxal Socket","Windung Coxa"
)
missing <- setdiff(required, names(df_raw))
if (length(missing) > 0) stop("Missing required columns in Excel: ", paste(missing, collapse = ", "))

df <- df_raw %>%
  dplyr::transmute(
    specimen_id = clean_str(.data[["specimen_id"]]),
    Family      = clean_str(.data[["Family"]]),
    Subfamily   = clean_str(.data[["Subfamily"]]),
    Genus       = clean_str(.data[["Genus"]]),
    Species     = clean_str(.data[["Species"]]),
    Tree_tip    = clean_str(.data[["Tree_tip"]]),
    Group       = clean_str(.data[["Group"]]),
    
    Schraube           = norm_tristate(.data[["Schraube"]]),
    Coxal_wall_hole          = norm_bi(.data[["Coxal wall hole"]]),
    Coxal_Socket   = norm_bi(.data[["Coxal Socket"]]),
    Windung_Coxa       = norm_bi(.data[["Windung Coxa"]])
  )

# -------------------- QC --------------------
# Duplicate specimen IDs
if (anyDuplicated(df$specimen_id) > 0) {
  dup_ids <- unique(df$specimen_id[duplicated(df$specimen_id)])
  write.csv2(data.frame(specimen_id = dup_ids),
             file.path(out_dir, "WARN_duplicate_specimen_id_DE.csv"),
             row.names = FALSE, fileEncoding = CSV_ENCODING)
  warning("Duplicate specimen_id found. See WARN_duplicate_specimen_id_DE.csv")
}

# Unexpected states
bad_schraube <- df %>% dplyr::filter(!is.na(Schraube) & !Schraube %in% c("TRUE","FALSE","X"))
bad_bin <- df %>% dplyr::filter(
  (!is.na(Coxal_wall_hole) & !Coxal_wall_hole %in% c("TRUE","FALSE")) |
    (!is.na(Coxal_Socket) & !Coxal_Socket %in% c("TRUE","FALSE")) |
    (!is.na(Windung_Coxa) & !Windung_Coxa %in% c("TRUE","FALSE"))
)
if (nrow(bad_schraube) > 0) write.csv2(bad_schraube, file.path(out_dir, "WARN_unexpected_states_schraube_DE.csv"),
                                       row.names = FALSE, fileEncoding = CSV_ENCODING)
if (nrow(bad_bin) > 0) write.csv2(bad_bin, file.path(out_dir, "WARN_unexpected_states_binary_DE.csv"),
                                  row.names = FALSE, fileEncoding = CSV_ENCODING)

# Keep only rows with complete traits + specimen_id
df2 <- df %>%
  dplyr::filter(
    specimen_id != "", !is.na(specimen_id),
    !is.na(Schraube),
    !is.na(Coxal_wall_hole),
    !is.na(Coxal_Socket),
    !is.na(Windung_Coxa)
  )

if (nrow(df2) == 0) stop("No rows left after filtering missing traits/specimen_id.")

# -------------------- CLASSIFICATION --------------------
# screw_state: clear vs ambiguous
df2 <- df2 %>%
  dplyr::mutate(
    screw_state = ifelse(Schraube == "X", "ambiguous", "clear")
  )

# joint_type (ASCII-safe labels; hyphen "-" only)
# Logic:
# - True screw-nut: trochanter screw present (TRUE or X) + coxal windung TRUE
# - Unopposed screw: trochanter screw present (TRUE or X) + coxal windung FALSE
# - Socket-guided rotational: no screw + no coxal windung + Einbuchtung TRUE
# - Unconstrained rotational: no screw + no coxal windung + Einbuchtung FALSE
# - Any Schraube FALSE but Windung_Coxa TRUE flagged as inconsistent

df2 <- df2 %>%
  dplyr::mutate(
    joint_type = dplyr::case_when(
      # clear screw
      Schraube == "TRUE"  & Windung_Coxa == "TRUE"  ~ "True screw-nut joint",
      Schraube == "TRUE"  & Windung_Coxa == "FALSE" ~ "Unopposed screw joint",
      
      # ambiguous screw (X): map to same mechanical types, but flagged
      Schraube == "X"     & Windung_Coxa == "TRUE"  ~ "True screw-nut joint",
      Schraube == "X"     & Windung_Coxa == "FALSE" ~ "Unopposed screw joint",
      
      # rotational types (no screw)
      Schraube == "FALSE" & Windung_Coxa == "FALSE" & Coxal_Socket == "TRUE"  ~ "Socket-guided rotational joint",
      Schraube == "FALSE" & Windung_Coxa == "FALSE" & Coxal_Socket == "FALSE" ~ "Unconstrained rotational joint",
      
      # logically inconsistent state
      Schraube == "FALSE" & Windung_Coxa == "TRUE" ~ "INCONSISTENT_coxal_windung_without_trochanter_screw",
      
      TRUE ~ "UNCLASSIFIED"
    ),
    
    # strict version for stats: drop ambiguous screw cases
    joint_type_strict = dplyr::if_else(screw_state == "ambiguous", NA_character_, joint_type)
  )

# Save weird cases explicitly
weird <- df2 %>% dplyr::filter(grepl("^INCONSISTENT|^UNCLASSIFIED$", joint_type))
if (nrow(weird) > 0) {
  write.csv2(weird, file.path(out_dir, "WARN_weird_or_inconsistent_cases_DE.csv"),
             row.names = FALSE, fileEncoding = CSV_ENCODING)
  warning("Some specimens are UNCLASSIFIED/INCONSISTENT. See WARN_weird_or_inconsistent_cases_DE.csv")
}

# -------------------- EXPORT mapping (DE CSV) --------------------
map <- df2 %>%
  dplyr::select(
    specimen_id, Family, Subfamily, Genus, Species, Tree_tip, Group,
    Schraube, Windung_Coxa, Coxal_Socket, Coxal_wall_hole,
    screw_state, joint_type, joint_type_strict
  ) %>%
  dplyr::arrange(Family, Genus, Species, specimen_id)

# German CSV export for Excel
write.csv2(map, out_map_csv, row.names = FALSE, fileEncoding = CSV_ENCODING)
message("Wrote German CSV mapping: ", out_map_csv)

# summary counts
summary_tab <- count_old(df2, joint_type, screw_state, out_name = "n") %>%
  dplyr::arrange(dplyr::desc(n))
write.csv2(summary_tab, file.path(out_dir, "SUMMARY_joint_type_counts_DE.csv"),
           row.names = FALSE, fileEncoding = CSV_ENCODING)

# -------------------- OPTIONAL: merge into specimen_key.csv --------------------
if (do_merge_specimen_key) {
  
  key_raw <- read_csv_any(specimen_key_path)
  names(key_raw) <- clean_str(names(key_raw))
  
  # find specimen_id column in key
  cand <- c("specimen_id", "specimen id", "id", "Specimen_ID", "SpecimenID")
  idx <- which(tolower(names(key_raw)) %in% tolower(cand))
  if (length(idx) == 0) {
    warning("Could not find a specimen_id column in specimen_key.csv (merge skipped).")
  } else {
    key_id_col <- names(key_raw)[idx[1]]
    key_raw[[key_id_col]] <- clean_str(key_raw[[key_id_col]])
    
    key_out <- key_raw %>%
      dplyr::left_join(
        map %>% dplyr::select(specimen_id, screw_state, joint_type, joint_type_strict),
        by = stats::setNames("specimen_id", key_id_col)
      )
    
    write.csv2(key_out, out_key_csv, row.names = FALSE, fileEncoding = CSV_ENCODING)
    message("Wrote merged specimen_key (DE CSV): ", out_key_csv)
  }
  
} else {
  message("specimen_key.csv not found at: ", specimen_key_path, " (merge skipped).")
}

message("Done. Check folder: ", out_dir)
