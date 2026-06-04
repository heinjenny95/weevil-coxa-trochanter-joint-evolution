# ============================================================
# Joint-combo audit (Curculionoidea screw-joint traits)
# Works with older dplyr versions (no count(name=...), no count(sort=...))
#
# Input columns EXACTLY:
# Family, Subfamily, Genus, Species, Tree_tip, specimen_id, Group,
# Schraube, Coxal wall hole, Coxal Socket, Windung Coxa
#
# Outputs:
# - all 24 theoretical combos, observed vs missing (overall)
# - counts + within-family proportions by Family
# - optional: by Subfamily and Group (if present/non-empty)
# - PDF plots
# ============================================================

rm(list = ls())

# -------------------- SETTINGS --------------------
in_xlsx <- "<BEETLE_JOINTS_ROOT>/joint_trait_table.xlsx"
out_dir <- file.path(dirname(in_xlsx), "joint_combo_audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------- PACKAGES --------------------
pkgs <- c("readxl", "dplyr", "tidyr", "stringr", "ggplot2", "forcats")
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

# Schraube: TRUE/FALSE/X
norm_tristate <- function(x) {
  y <- clean_str(x)
  y_low <- tolower(y)
  y_low[y_low %in% c("", "na", "n/a", "null", "none")] <- NA_character_
  
  out <- dplyr::case_when(
    is.na(y_low) ~ NA_character_,
    y_low %in% c("true","t","1","yes","y","ja","wahr") ~ "TRUE",
    y_low %in% c("false","f","0","no","n","nein","falsch") ~ "FALSE",
    y_low %in% c("x","unklar","ambiguous","intermediate","maybe") ~ "X",
    TRUE ~ y  # unknown -> keep as-is (will be flagged)
  )
  out
}

# Binary traits: TRUE/FALSE
norm_bi <- function(x) {
  y <- clean_str(x)
  y_low <- tolower(y)
  y_low[y_low %in% c("", "na", "n/a", "null", "none")] <- NA_character_
  
  out <- dplyr::case_when(
    is.na(y_low) ~ NA_character_,
    y_low %in% c("true","t","1","yes","y","ja","wahr") ~ "TRUE",
    y_low %in% c("false","f","0","no","n","nein","falsch") ~ "FALSE",
    TRUE ~ y
  )
  out
}

# Older dplyr-safe "count"
count_old <- function(data, ..., out_name = "n") {
  # returns a data.frame with columns (...) and out_name
  data %>%
    dplyr::group_by(...) %>%
    dplyr::summarise(n_tmp = dplyr::n(), .groups = "drop") %>%
    dplyr::rename(!!out_name := n_tmp)
}

# -------------------- READ EXCEL --------------------
df_raw <- readxl::read_excel(in_xlsx)
if (nrow(df_raw) == 0) stop("Excel loaded but has 0 rows: ", in_xlsx)

# Normalize column names minimally (NBSP etc.) but keep exact visible spelling
names(df_raw) <- clean_str(names(df_raw))

required <- c(
  "Family","Subfamily","Genus","Species","Tree_tip","specimen_id","Group",
  "Schraube","Coxal wall hole","Coxal Socket","Windung Coxa"
)

missing <- setdiff(required, names(df_raw))
if (length(missing) > 0) {
  stop("Missing required columns in Excel: ", paste(missing, collapse = ", "))
}

df <- df_raw %>%
  dplyr::transmute(
    Family      = clean_str(.data[["Family"]]),
    Subfamily   = clean_str(.data[["Subfamily"]]),
    Genus       = clean_str(.data[["Genus"]]),
    Species     = clean_str(.data[["Species"]]),
    Tree_tip    = clean_str(.data[["Tree_tip"]]),
    specimen_id = clean_str(.data[["specimen_id"]]),
    Group       = clean_str(.data[["Group"]]),
    
    Schraube            = norm_tristate(.data[["Schraube"]]),
    `Coxal wall hole`         = norm_bi(.data[["Coxal wall hole"]]),
    `Coxal Socket`  = norm_bi(.data[["Coxal Socket"]]),
    `Windung Coxa`      = norm_bi(.data[["Windung Coxa"]])
  )

# -------------------- QUALITY CHECKS --------------------
# Duplicate specimen IDs
if (anyDuplicated(df$specimen_id) > 0) {
  dup_ids <- unique(df$specimen_id[duplicated(df$specimen_id)])
  write.csv(
    data.frame(specimen_id = dup_ids),
    file.path(out_dir, "WARN_duplicate_specimen_id.csv"),
    row.names = FALSE
  )
  warning("Duplicate specimen_id found. See WARN_duplicate_specimen_id.csv")
}

# Unexpected states
bad_schraube <- df %>% dplyr::filter(!is.na(Schraube) & !Schraube %in% c("TRUE","FALSE","X"))
bad_bin <- df %>% dplyr::filter(
  (!is.na(`Coxal wall hole`) & !`Coxal wall hole` %in% c("TRUE","FALSE")) |
    (!is.na(`Coxal Socket`) & !`Coxal Socket` %in% c("TRUE","FALSE")) |
    (!is.na(`Windung Coxa`) & !`Windung Coxa` %in% c("TRUE","FALSE"))
)

if (nrow(bad_schraube) > 0) write.csv(bad_schraube, file.path(out_dir, "WARN_unexpected_states_schraube.csv"), row.names = FALSE)
if (nrow(bad_bin) > 0) write.csv(bad_bin, file.path(out_dir, "WARN_unexpected_states_binary.csv"), row.names = FALSE)

# Keep only complete rows for combo enumeration (traits + Family present)
df_complete <- df %>%
  dplyr::filter(
    !is.na(Family), Family != "",
    !is.na(Schraube),
    !is.na(`Coxal wall hole`),
    !is.na(`Coxal Socket`),
    !is.na(`Windung Coxa`)
  )

write.csv(df_complete, file.path(out_dir, "DATA_used_complete_rows.csv"), row.names = FALSE)
if (nrow(df_complete) == 0) stop("No complete rows left after filtering NAs in traits.")

# -------------------- BUILD COMBOS --------------------
df_complete <- df_complete %>%
  dplyr::mutate(
    Schraube = factor(Schraube, levels = c("FALSE","X","TRUE")),
    `Coxal wall hole` = factor(`Coxal wall hole`, levels = c("FALSE","TRUE")),
    `Coxal Socket` = factor(`Coxal Socket`, levels = c("FALSE","TRUE")),
    `Windung Coxa` = factor(`Windung Coxa`, levels = c("FALSE","TRUE")),
    combo = interaction(
      Schraube, `Coxal wall hole`, `Coxal Socket`, `Windung Coxa`,
      sep = "_", drop = FALSE
    )
  )

# Theoretical space = 24 combos
all_combos <- expand.grid(
  Schraube = factor(c("FALSE","X","TRUE"), levels = c("FALSE","X","TRUE")),
  `Coxal wall hole` = factor(c("FALSE","TRUE"), levels = c("FALSE","TRUE")),
  `Coxal Socket` = factor(c("FALSE","TRUE"), levels = c("FALSE","TRUE")),
  `Windung Coxa` = factor(c("FALSE","TRUE"), levels = c("FALSE","TRUE")),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    combo = interaction(
      Schraube, `Coxal wall hole`, `Coxal Socket`, `Windung Coxa`,
      sep = "_", drop = FALSE
    )
  )

# -------------------- OVERALL OBSERVED VS MISSING --------------------
obs_overall <- count_old(df_complete, combo, out_name = "n_observed") %>%
  dplyr::right_join(all_combos %>% dplyr::select(combo), by = "combo") %>%
  dplyr::mutate(n_observed = tidyr::replace_na(n_observed, 0L)) %>%
  dplyr::left_join(all_combos, by = "combo") %>%
  dplyr::mutate(is_missing = (n_observed == 0L)) %>%
  dplyr::arrange(is_missing, dplyr::desc(n_observed), combo)

write.csv(obs_overall, file.path(out_dir, "combos_overall_observed_vs_missing.csv"), row.names = FALSE)

# Missing-only list
write.csv(
  obs_overall %>%
    dplyr::filter(is_missing) %>%
    dplyr::select(combo, Schraube, `Coxal wall hole`, `Coxal Socket`, `Windung Coxa`),
  file.path(out_dir, "combos_missing_only.csv"),
  row.names = FALSE
)

# -------------------- BY FAMILY (COUNTS + WITHIN-FAMILY PROPS) --------------------
fam_levels <- sort(unique(df_complete$Family))

obs_by_family <- count_old(df_complete, Family, combo, out_name = "n") %>%
  dplyr::right_join(
    tidyr::expand_grid(Family = fam_levels, combo = all_combos$combo),
    by = c("Family","combo")
  ) %>%
  dplyr::mutate(n = tidyr::replace_na(n, 0L)) %>%
  dplyr::left_join(all_combos, by = "combo") %>%
  dplyr::group_by(Family) %>%
  dplyr::mutate(
    n_family_total = sum(n),
    prop_in_family = ifelse(n_family_total > 0, n / n_family_total, NA_real_)
  ) %>%
  dplyr::ungroup()

write.csv(obs_by_family, file.path(out_dir, "combos_by_family_counts_and_props.csv"), row.names = FALSE)

presence_by_family <- obs_by_family %>%
  dplyr::mutate(present = n > 0) %>%
  dplyr::select(Family, combo, present, Schraube, `Coxal wall hole`, `Coxal Socket`, `Windung Coxa`, n, prop_in_family)

write.csv(presence_by_family, file.path(out_dir, "combos_by_family_presence.csv"), row.names = FALSE)

# -------------------- OPTIONAL: BY SUBFAMILY / GROUP --------------------
# Subfamily
sub_levels <- sort(unique(df_complete$Subfamily[df_complete$Subfamily != ""]))
if (length(sub_levels) > 0) {
  obs_by_subfamily <- df_complete %>%
    dplyr::filter(Subfamily != "") %>%
    count_old(Subfamily, combo, out_name = "n") %>%
    dplyr::right_join(
      tidyr::expand_grid(Subfamily = sub_levels, combo = all_combos$combo),
      by = c("Subfamily","combo")
    ) %>%
    dplyr::mutate(n = tidyr::replace_na(n, 0L)) %>%
    dplyr::left_join(all_combos, by = "combo") %>%
    dplyr::group_by(Subfamily) %>%
    dplyr::mutate(
      n_subfamily_total = sum(n),
      prop_in_subfamily = ifelse(n_subfamily_total > 0, n / n_subfamily_total, NA_real_)
    ) %>%
    dplyr::ungroup()
  
  write.csv(obs_by_subfamily, file.path(out_dir, "combos_by_subfamily_counts_and_props.csv"), row.names = FALSE)
}

# Group
grp_levels <- sort(unique(df_complete$Group[df_complete$Group != ""]))
if (length(grp_levels) > 0) {
  obs_by_group <- df_complete %>%
    dplyr::filter(Group != "") %>%
    count_old(Group, combo, out_name = "n") %>%
    dplyr::right_join(
      tidyr::expand_grid(Group = grp_levels, combo = all_combos$combo),
      by = c("Group","combo")
    ) %>%
    dplyr::mutate(n = tidyr::replace_na(n, 0L)) %>%
    dplyr::left_join(all_combos, by = "combo") %>%
    dplyr::group_by(Group) %>%
    dplyr::mutate(
      n_group_total = sum(n),
      prop_in_group = ifelse(n_group_total > 0, n / n_group_total, NA_real_)
    ) %>%
    dplyr::ungroup()
  
  write.csv(obs_by_group, file.path(out_dir, "combos_by_group_counts_and_props.csv"), row.names = FALSE)
}

# -------------------- PLOTS (PDF) --------------------
# 1) Overall combo counts (including missing)
p1 <- obs_overall %>%
  dplyr::mutate(combo = forcats::fct_reorder(combo, n_observed, .desc = TRUE)) %>%
  ggplot(aes(x = combo, y = n_observed)) +
  geom_col() +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(
    title = "Observed joint-trait combinations (overall)",
    subtitle = "Full theoretical space (24 combos). Bars at 0 = not observed.",
    x = "Combination (Schraube_CoxalWallHole_CoxalSocket_WindungCoxa)",
    y = "Count"
  )
ggsave(file.path(out_dir, "PLOT_overall_combo_counts.pdf"), p1, width = 14, height = 6)

# 2) Family x combo heatmap (proportions)
p2_data <- obs_by_family %>% dplyr::filter(n_family_total > 0)

p2 <- p2_data %>%
  dplyr::mutate(
    combo = forcats::fct_reorder(combo, n, .desc = TRUE),
    Family = forcats::fct_reorder(Family, n_family_total, .desc = TRUE)
  ) %>%
  ggplot(aes(x = combo, y = Family, fill = prop_in_family)) +
  geom_tile() +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(
    title = "Joint-trait combination composition by family",
    subtitle = "Fill = within-family proportion (n / family total). Empty = 0.",
    x = "Combination",
    y = "Family",
    fill = "Proportion"
  )

ggsave(
  file.path(out_dir, "PLOT_family_combo_heatmap_props.pdf"),
  p2,
  width = 14,
  height = max(6, 0.25 * length(fam_levels))
)

# 3) Family x (Loch/Einbuchtung/Windung) faceted by Schraube (counts)
p3 <- obs_by_family %>%
  dplyr::filter(n_family_total > 0) %>%
  dplyr::mutate(
    combo_short = interaction(`Coxal wall hole`, `Coxal Socket`, `Windung Coxa`, sep = "_", drop = FALSE),
    combo_short = forcats::fct_reorder(combo_short, n, .fun = sum, .desc = TRUE),
    Family = forcats::fct_reorder(Family, n_family_total, .desc = TRUE)
  ) %>%
  ggplot(aes(x = combo_short, y = Family, fill = n)) +
  geom_tile() +
  facet_wrap(~ Schraube, nrow = 1) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(
    title = "Family distribution split by Schraube state",
    subtitle = "Tiles show counts per family and (Loch_Einbuchtung_Windung) combo.",
    x = "Loch_Einbuchtung_Windung",
    y = "Family",
    fill = "Count"
  )

ggsave(
  file.path(out_dir, "PLOT_family_heatmap_by_schraube_counts.pdf"),
  p3,
  width = 14,
  height = max(6, 0.25 * length(fam_levels))
)

# -------------------- SUMMARY --------------------
n_obs <- sum(obs_overall$n_observed > 0)
n_miss <- sum(obs_overall$n_observed == 0)

summary_lines <- c(
  paste0("Input: ", in_xlsx),
  paste0("Rows (raw): ", nrow(df_raw)),
  paste0("Rows (complete traits + family): ", nrow(df_complete)),
  paste0("Theoretical combos: 24"),
  paste0("Observed combos: ", n_obs),
  paste0("Missing combos: ", n_miss),
  "",
  "Top 10 observed combos (combo : n):",
  paste0(capture.output(print(
    head(obs_overall %>% dplyr::filter(n_observed > 0) %>% dplyr::select(combo, n_observed), 10),
    row.names = FALSE
  )), collapse = "\n"),
  "",
  "Output folder:",
  out_dir
)

writeLines(summary_lines, file.path(out_dir, "SUMMARY.txt"))

message("Done. Outputs in: ", out_dir)
