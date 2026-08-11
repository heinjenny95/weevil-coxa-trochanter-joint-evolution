# Propagate conditional robust-helix fit uncertainty through the principal
# downstream association tests. Bootstrap replicate numbers are aligned across
# specimens, yielding 200 complete measurement-realisation datasets.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ape)
  library(nlme)
  library(vegan)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop(
    "Usage: Rscript propagate_robust_geometry_uncertainty.R ",
    "<metrics.csv> <bootstrap_draws.csv> <PCA.csv> <joint_types.csv> ",
    "<specimen_key.csv> <tree.tre> <ecology.csv> <output_dir>"
  )
}

metrics_path <- args[[1]]
draws_path <- args[[2]]
pca_path <- args[[3]]
joint_path <- args[[4]]
key_path <- args[[5]]
tree_path <- args[[6]]
ecology_path <- args[[7]]
output_dir <- args[[8]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260810)

id_aliases <- c(
  "308_lisshorhoptrus_oryzophilus_aligned" = "308_lisshorhoptrus_oryzophilus_trochanter_aligned",
  "310_ormiscus_saltator_trochanter_mirrored_aligned" = "310_ormiscus_saltator_trochanter_aligned",
  "t_pseudonastus_trochanter_aligned" = "t_pseudonasutus_trochanter_aligned"
)

canonical_id <- function(x) {
  x <- trimws(as.character(x))
  hit <- x %in% names(id_aliases)
  x[hit] <- unname(id_aliases[x[hit]])
  x
}

write_clean <- function(x, filename) {
  write.csv(x, file.path(output_dir, filename), row.names = FALSE, na = "")
}

metrics <- read.csv(metrics_path, stringsAsFactors = FALSE, check.names = FALSE)
draws <- read.csv(draws_path, stringsAsFactors = FALSE, check.names = FALSE)
pca <- read.csv2(pca_path, stringsAsFactors = FALSE, check.names = FALSE)
joint <- read.csv(joint_path, stringsAsFactors = FALSE, check.names = FALSE)
specimen_key <- read.csv2(key_path, stringsAsFactors = FALSE, check.names = FALSE)
ecology <- read.csv(ecology_path, stringsAsFactors = FALSE, check.names = FALSE)
tree <- ape::read.tree(tree_path)
tree$node.label <- NULL

metrics$specimen_id <- canonical_id(metrics$specimen_id)
draws$specimen_id <- canonical_id(draws$specimen_id)
pca$specimen_id <- canonical_id(pca$specimen_id)
joint$specimen_id <- canonical_id(joint$specimen_id)
specimen_key$specimen_id <- canonical_id(specimen_key$specimen_id)

numeric_columns <- c(
  "abs_winding_angle_deg", "fitted_pitch_360",
  "endpoint_axial_span_new_axis", "helix_rms_relative_to_radius"
)
for (column in intersect(numeric_columns, names(metrics))) {
  metrics[[column]] <- as.numeric(metrics[[column]])
}
draw_numeric <- c(
  "bootstrap_replicate", "abs_winding_angle_deg", "fitted_pitch_360",
  "endpoint_axial_span", "fit_radius", "helix_rms"
)
for (column in intersect(draw_numeric, names(draws))) {
  draws[[column]] <- as.numeric(draws[[column]])
}
for (column in c("PC1", "PC2")) pca[[column]] <- as.numeric(pca[[column]])

metrics$analysis_set_primary_adequate <- as.logical(metrics$analysis_set_primary_adequate)
metrics$analysis_set_strict_good <- as.logical(metrics$analysis_set_strict_good)

if (nrow(metrics) != 64 || nrow(draws) != 12800) {
  stop("Unexpected robust geometry dimensions")
}
if (anyDuplicated(metrics$specimen_id)) stop("Duplicate canonical specimen IDs in metrics")
if (anyDuplicated(pca$specimen_id)) stop("Duplicate specimen IDs in PCA table")

analysis_sets <- list(
  all = metrics$specimen_id,
  primary_adequate = metrics$specimen_id[metrics$analysis_set_primary_adequate],
  strict_good = metrics$specimen_id[metrics$analysis_set_strict_good]
)

point_geometry <- metrics %>%
  transmute(
    specimen_id,
    angle = abs_winding_angle_deg,
    pitch = fitted_pitch_360,
    axial_span = endpoint_axial_span_new_axis
  )

draw_geometry <- draws %>%
  transmute(
    specimen_id,
    replicate = as.integer(bootstrap_replicate),
    angle = abs_winding_angle_deg,
    pitch = fitted_pitch_360,
    axial_span = endpoint_axial_span
  )

geometry_for <- function(set_name, replicate = 0L) {
  ids <- analysis_sets[[set_name]]
  if (replicate == 0L) {
    point_geometry %>% filter(specimen_id %in% ids)
  } else {
    draw_geometry %>% filter(specimen_id %in% ids, .data$replicate == .env$replicate)
  }
}

model_statistics <- function(model) {
  sm <- summary(model)
  fstat <- sm$fstatistic
  data.frame(
    n = nobs(model),
    r_squared = sm$r.squared,
    adjusted_r_squared = sm$adj.r.squared,
    f_statistic = unname(fstat[[1]]),
    df1 = unname(fstat[[2]]),
    df2 = unname(fstat[[3]]),
    model_p = pf(fstat[[1]], fstat[[2]], fstat[[3]], lower.tail = FALSE)
  )
}

coefficient_rows <- function(model) {
  table <- as.data.frame(summary(model)$coefficients)
  table$term <- rownames(table)
  rownames(table) <- NULL
  table %>%
    transmute(
      term,
      estimate = Estimate,
      std_error = `Std. Error`,
      statistic = `t value`,
      p_value = `Pr(>|t|)`
    )
}

rubin_combine <- function(estimates, standard_errors) {
  keep <- is.finite(estimates) & is.finite(standard_errors)
  estimates <- estimates[keep]
  standard_errors <- standard_errors[keep]
  m <- length(estimates)
  if (m < 2) {
    return(data.frame(
      m = m, estimate = NA_real_, std_error = NA_real_, df = NA_real_,
      statistic = NA_real_, p_value = NA_real_, ci_low = NA_real_,
      ci_high = NA_real_, within_variance = NA_real_, between_variance = NA_real_
    ))
  }
  qbar <- mean(estimates)
  ubar <- mean(standard_errors^2)
  between <- stats::var(estimates)
  total <- ubar + (1 + 1 / m) * between
  se_total <- sqrt(total)
  df <- if (between <= .Machine$double.eps) Inf else {
    (m - 1) * (1 + ubar / ((1 + 1 / m) * between))^2
  }
  statistic <- qbar / se_total
  p_value <- 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE)
  critical <- stats::qt(0.975, df = df)
  data.frame(
    m = m, estimate = qbar, std_error = se_total, df = df,
    statistic = statistic, p_value = p_value,
    ci_low = qbar - critical * se_total,
    ci_high = qbar + critical * se_total,
    within_variance = ubar, between_variance = between
  )
}

# -------------------------------------------------------------------------
# Specimen-level shape--geometry models
# -------------------------------------------------------------------------

shape_definitions <- data.frame(
  model = c("angle_full", "pitch_full", "axial_span_full", "angle_main_region"),
  response = c("angle", "pitch", "axial_span", "angle"),
  subset = c("full", "full", "full", "PC1_lt_0.1"),
  stringsAsFactors = FALSE
)

shape_model_rows <- list()
shape_coefficient_rows <- list()
shape_index <- 1L
coefficient_index <- 1L

for (set_name in names(analysis_sets)) {
  for (replicate in 0:200) {
    geometry <- geometry_for(set_name, replicate)
    dat_base <- inner_join(
      pca %>% select(specimen_id, PC1, PC2), geometry, by = "specimen_id"
    )
    for (definition_index in seq_len(nrow(shape_definitions))) {
      definition <- shape_definitions[definition_index, ]
      dat <- dat_base
      if (definition$subset == "PC1_lt_0.1") dat <- dat %>% filter(PC1 < 0.1)
      dat <- dat %>%
        filter(if_all(all_of(c(definition$response, "PC1", "PC2")), is.finite))
      model <- lm(reformulate(c("PC1", "PC2"), response = definition$response), data = dat)
      metadata <- data.frame(
        analysis_set = set_name,
        replicate = replicate,
        draw_type = ifelse(replicate == 0, "point_estimate", "measurement_bootstrap"),
        model = definition$model,
        response = definition$response,
        subset = definition$subset,
        stringsAsFactors = FALSE
      )
      shape_model_rows[[shape_index]] <- bind_cols(metadata, model_statistics(model))
      shape_index <- shape_index + 1L
      shape_coefficient_rows[[coefficient_index]] <- bind_cols(metadata, coefficient_rows(model))
      coefficient_index <- coefficient_index + 1L
    }
  }
}

shape_models <- bind_rows(shape_model_rows)
shape_coefficients <- bind_rows(shape_coefficient_rows)
write_clean(shape_models, "shape_model_draws.csv")
write_clean(shape_coefficients, "shape_coefficient_draws.csv")

shape_model_summary <- shape_models %>%
  filter(replicate > 0) %>%
  group_by(analysis_set, model, response, subset) %>%
  summarise(
    bootstrap_replicates = n(),
    r_squared_median = median(r_squared),
    r_squared_ci_low = quantile(r_squared, 0.025),
    r_squared_ci_high = quantile(r_squared, 0.975),
    model_p_median = median(model_p),
    model_p_ci_low = quantile(model_p, 0.025),
    model_p_ci_high = quantile(model_p, 0.975),
    fraction_model_p_below_0.05 = mean(model_p < 0.05),
    .groups = "drop"
  ) %>%
  left_join(
    shape_models %>%
      filter(replicate == 0) %>%
      select(analysis_set, model, point_n = n, point_r_squared = r_squared, point_model_p = model_p),
    by = c("analysis_set", "model")
  )
write_clean(shape_model_summary, "shape_model_uncertainty_summary.csv")

shape_rubin <- shape_coefficients %>%
  filter(replicate > 0, term != "(Intercept)") %>%
  group_by(analysis_set, model, response, subset, term) %>%
  group_modify(~ rubin_combine(.x$estimate, .x$std_error)) %>%
  ungroup() %>%
  left_join(
    shape_coefficients %>%
      filter(replicate == 0, term != "(Intercept)") %>%
      select(
        analysis_set, model, term,
        point_estimate = estimate, point_std_error = std_error,
        point_p_value = p_value
      ),
    by = c("analysis_set", "model", "term")
  )
write_clean(shape_rubin, "shape_coefficients_rubin_summary.csv")

# -------------------------------------------------------------------------
# Joint-type tests
# -------------------------------------------------------------------------

is_screw_joint <- function(x) {
  x %in% c("True screw-nut joint", "Unopposed screw configuration")
}

joint_test_rows <- list()
joint_index <- 1L
for (set_name in names(analysis_sets)) {
  for (replicate in 0:200) {
    dat <- geometry_for(set_name, replicate) %>%
      inner_join(joint %>% select(specimen_id, joint_type), by = "specimen_id") %>%
      filter(is_screw_joint(joint_type)) %>%
      mutate(joint_type = factor(joint_type))
    group_counts <- table(dat$joint_type)
    minimum_group_n <- if (length(group_counts)) min(group_counts) else 0
    for (trait in c("angle", "pitch", "axial_span")) {
      result <- data.frame(
        analysis_set = set_name, replicate = replicate,
        draw_type = ifelse(replicate == 0, "point_estimate", "measurement_bootstrap"),
        analysis = "kruskal_wallis", response = trait,
        n = nrow(dat), n_groups = length(group_counts),
        minimum_group_n = minimum_group_n,
        statistic = NA_real_, p_value = NA_real_, r_squared = NA_real_,
        status = "insufficient_group_size", stringsAsFactors = FALSE
      )
      if (length(group_counts) >= 2 && minimum_group_n >= 2) {
        test <- kruskal.test(dat[[trait]] ~ dat$joint_type)
        result$statistic <- unname(test$statistic)
        result$p_value <- test$p.value
        result$status <- "ok"
      }
      joint_test_rows[[joint_index]] <- result
      joint_index <- joint_index + 1L
    }
    permanova_result <- data.frame(
      analysis_set = set_name, replicate = replicate,
      draw_type = ifelse(replicate == 0, "point_estimate", "measurement_bootstrap"),
      analysis = "PERMANOVA", response = "angle_pitch_axial_span",
      n = nrow(dat), n_groups = length(group_counts),
      minimum_group_n = minimum_group_n,
      statistic = NA_real_, p_value = NA_real_, r_squared = NA_real_,
      status = "insufficient_group_size", stringsAsFactors = FALSE
    )
    if (length(group_counts) >= 2 && minimum_group_n >= 2) {
      matrix_scaled <- scale(dat[, c("angle", "pitch", "axial_span")])
      distance <- dist(matrix_scaled)
      permutations <- ifelse(replicate == 0, 999, 199)
      set.seed(20260810 + replicate)
      fit <- vegan::adonis2(distance ~ joint_type, data = dat, permutations = permutations)
      permanova_result$statistic <- fit$F[[1]]
      permanova_result$p_value <- fit$`Pr(>F)`[[1]]
      permanova_result$r_squared <- fit$R2[[1]]
      permanova_result$status <- "ok"
    }
    joint_test_rows[[joint_index]] <- permanova_result
    joint_index <- joint_index + 1L
  }
}

joint_tests <- bind_rows(joint_test_rows)
write_clean(joint_tests, "joint_type_test_draws.csv")
joint_summary <- joint_tests %>%
  filter(replicate > 0) %>%
  group_by(analysis_set, analysis, response, status) %>%
  summarise(
    bootstrap_replicates = n(),
    statistic_median = ifelse(all(is.na(statistic)), NA_real_, median(statistic, na.rm = TRUE)),
    statistic_ci_low = ifelse(all(is.na(statistic)), NA_real_, quantile(statistic, 0.025, na.rm = TRUE)),
    statistic_ci_high = ifelse(all(is.na(statistic)), NA_real_, quantile(statistic, 0.975, na.rm = TRUE)),
    p_value_median = ifelse(all(is.na(p_value)), NA_real_, median(p_value, na.rm = TRUE)),
    p_value_ci_low = ifelse(all(is.na(p_value)), NA_real_, quantile(p_value, 0.025, na.rm = TRUE)),
    p_value_ci_high = ifelse(all(is.na(p_value)), NA_real_, quantile(p_value, 0.975, na.rm = TRUE)),
    fraction_p_below_0.05 = ifelse(all(is.na(p_value)), NA_real_, mean(p_value < 0.05, na.rm = TRUE)),
    minimum_group_n = min(minimum_group_n),
    .groups = "drop"
  ) %>%
  left_join(
    joint_tests %>%
      filter(replicate == 0) %>%
      select(
        analysis_set, analysis, response,
        point_statistic = statistic, point_p_value = p_value,
        point_r_squared = r_squared, point_status = status
      ),
    by = c("analysis_set", "analysis", "response")
  )
write_clean(joint_summary, "joint_type_uncertainty_summary.csv")

# -------------------------------------------------------------------------
# Tip aggregation and bounded-lambda PGLS
# -------------------------------------------------------------------------

correct_tree_tip <- function(family, tree_tip) {
  output <- as.character(tree_tip)
  output[output == "Neydus"] <- "Nedyus"
  output[family == "Belidae"] <- "Agnesiotis"
  output[family == "Caridae"] <- "Car"
  output
}

build_tip_data <- function(set_name, replicate) {
  geometry_for(set_name, replicate) %>%
    inner_join(specimen_key %>% select(specimen_id, Family, tree_tip), by = "specimen_id") %>%
    inner_join(pca %>% select(specimen_id, PC1, PC2), by = "specimen_id") %>%
    mutate(
      tree_tip = correct_tree_tip(Family, tree_tip),
      tree_label = paste(Family, tree_tip, sep = "___")
    ) %>%
    filter(tree_label %in% tree$tip.label) %>%
    group_by(tree_label, Family, tree_tip) %>%
    summarise(
      PC1 = mean(PC1, na.rm = TRUE),
      PC2 = mean(PC2, na.rm = TRUE),
      angle = mean(angle, na.rm = TRUE),
      pitch = mean(pitch, na.rm = TRUE),
      axial_span = mean(axial_span, na.rm = TRUE),
      n_specimens = n(),
      .groups = "drop"
    )
}

fit_gls_at_lambda <- function(dat, response, predictor, lambda, factor_predictor = FALSE) {
  needed <- c("tree_label", response, predictor)
  model_dat <- dat %>%
    select(all_of(needed)) %>%
    filter(if_all(all_of(needed), ~ !is.na(.x))) %>%
    distinct()
  if (nrow(model_dat) < 5) return(NULL)
  if (factor_predictor) {
    model_dat[[predictor]] <- factor(model_dat[[predictor]])
    counts <- table(model_dat[[predictor]])
    if (length(counts) < 2 || any(counts < 2)) return(NULL)
  }
  phy <- drop.tip(tree, setdiff(tree$tip.label, model_dat$tree_label))
  model_dat <- model_dat[match(phy$tip.label, model_dat$tree_label), , drop = FALSE]
  correlation <- tryCatch(
    ape::corPagel(lambda, phy = phy, fixed = TRUE, form = ~tree_label),
    error = function(e) NULL
  )
  if (is.null(correlation)) return(NULL)
  formula <- reformulate(predictor, response = response)
  fit <- tryCatch(
    suppressWarnings(nlme::gls(
      formula, data = as.data.frame(model_dat), correlation = correlation,
      method = "ML", na.action = na.omit
    )),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  list(fit = fit, n_taxa = nrow(model_dat))
}

optimize_lambda <- function(dat, response, predictor, factor_predictor = FALSE) {
  objective <- function(lambda) {
    fitted <- fit_gls_at_lambda(dat, response, predictor, lambda, factor_predictor)
    if (is.null(fitted)) return(Inf)
    value <- tryCatch(-as.numeric(logLik(fitted$fit)), error = function(e) Inf)
    ifelse(is.finite(value), value, Inf)
  }
  optimized <- tryCatch(optimize(objective, c(0, 1), tol = 1e-7), error = function(e) NULL)
  candidates <- c(0, 1)
  if (!is.null(optimized) && is.finite(optimized$minimum)) {
    candidates <- c(candidates, optimized$minimum)
  }
  values <- vapply(candidates, objective, numeric(1))
  if (all(!is.finite(values))) return(NA_real_)
  candidates[which.min(values)]
}

extract_gls_rows <- function(fitted, response, predictor, lambda) {
  if (is.null(fitted)) return(data.frame())
  table <- as.data.frame(summary(fitted$fit)$tTable)
  table$term <- rownames(table)
  rownames(table) <- NULL
  table %>%
    filter(term != "(Intercept)") %>%
    transmute(
      response, predictor, term,
      estimate = Value,
      std_error = `Std.Error`,
      statistic = `t-value`,
      p_value = `p-value`,
      lambda = lambda,
      n_taxa = fitted$n_taxa
    )
}

pgls_models <- expand.grid(
  response = c("PC1", "PC2"),
  predictor = c("angle", "pitch", "axial_span"),
  stringsAsFactors = FALSE
)
pgls_rows <- list()
pgls_index <- 1L
for (set_name in names(analysis_sets)) {
  point_tip <- build_tip_data(set_name, 0L)
  for (model_index in seq_len(nrow(pgls_models))) {
    response <- pgls_models$response[[model_index]]
    predictor <- pgls_models$predictor[[model_index]]
    lambda <- optimize_lambda(point_tip, response, predictor)
    if (!is.finite(lambda)) next
    for (replicate in 0:200) {
      dat <- if (replicate == 0) point_tip else build_tip_data(set_name, replicate)
      extracted <- extract_gls_rows(
        fit_gls_at_lambda(dat, response, predictor, lambda),
        response, predictor, lambda
      )
      if (nrow(extracted)) {
        pgls_rows[[pgls_index]] <- extracted %>%
          mutate(
            analysis_set = set_name,
            replicate = replicate,
            draw_type = ifelse(replicate == 0, "point_estimate", "measurement_bootstrap"),
            lambda_strategy = "point_estimate_lambda_fixed_across_measurement_draws"
          )
        pgls_index <- pgls_index + 1L
      }
    }
  }
}
pgls_draws <- bind_rows(pgls_rows)
write_clean(pgls_draws, "pgls_shape_geometry_draws.csv")
pgls_rubin <- pgls_draws %>%
  filter(replicate > 0) %>%
  group_by(analysis_set, response, predictor, term, lambda, lambda_strategy) %>%
  group_modify(~ rubin_combine(.x$estimate, .x$std_error)) %>%
  ungroup() %>%
  left_join(
    pgls_draws %>%
      filter(replicate == 0) %>%
      select(
        analysis_set, response, predictor, term,
        point_estimate = estimate, point_std_error = std_error,
        point_p_value = p_value, point_n_taxa = n_taxa
      ),
    by = c("analysis_set", "response", "predictor", "term")
  )
write_clean(pgls_rubin, "pgls_shape_geometry_rubin_summary.csv")

# -------------------------------------------------------------------------
# Ecology factor PGLS for the two geometry responses used historically
# -------------------------------------------------------------------------

ecology2 <- ecology %>%
  mutate(
    host_lineage_broad = case_when(
      host_lineage_simple == "angiosperm" ~ "angiosperm",
      host_lineage_simple == "gymnosperm" ~ "gymnosperm",
      TRUE ~ NA_character_
    ),
    woody_association_broad = case_when(
      woody_association_simple %in% c("woody", "woody_or_dead_plant_associated") ~ "woody",
      woody_association_simple %in% c("mostly_nonwoody", "nonwoody_or_shrub_associated") ~ "nonwoody",
      TRUE ~ NA_character_
    ),
    larval_lifestyle_broad = case_when(
      larval_lifestyle_simple == "internal" ~ "internal",
      larval_lifestyle_simple %in% c("mixed", "external_root_feeding") ~ "other",
      TRUE ~ NA_character_
    ),
    fungal_association_broad = case_when(
      fungal_symbiosis_simple %in% c("yes", "yes_or_common") ~ "yes",
      fungal_symbiosis_simple == "no" ~ "no",
      TRUE ~ NA_character_
    )
  )

ecology_predictors <- c(
  "host_lineage_broad", "woody_association_broad",
  "larval_lifestyle_broad", "fungal_association_broad"
)
ecology_rows <- list()
ecology_index <- 1L
for (set_name in c("primary_adequate", "strict_good")) {
  point_tip <- build_tip_data(set_name, 0L) %>% left_join(ecology2, by = "tree_label")
  for (response in c("angle", "axial_span")) {
    for (predictor in ecology_predictors) {
      lambda <- optimize_lambda(point_tip, response, predictor, factor_predictor = TRUE)
      if (!is.finite(lambda)) next
      for (replicate in 0:200) {
        dat <- if (replicate == 0) point_tip else {
          build_tip_data(set_name, replicate) %>% left_join(ecology2, by = "tree_label")
        }
        extracted <- extract_gls_rows(
          fit_gls_at_lambda(dat, response, predictor, lambda, factor_predictor = TRUE),
          response, predictor, lambda
        )
        if (nrow(extracted)) {
          ecology_rows[[ecology_index]] <- extracted %>%
            mutate(
              analysis_set = set_name,
              replicate = replicate,
              draw_type = ifelse(replicate == 0, "point_estimate", "measurement_bootstrap"),
              lambda_strategy = "point_estimate_lambda_fixed_across_measurement_draws"
            )
          ecology_index <- ecology_index + 1L
        }
      }
    }
  }
}
ecology_draws <- bind_rows(ecology_rows)
write_clean(ecology_draws, "ecology_pgls_draws.csv")
if (nrow(ecology_draws)) {
  ecology_rubin <- ecology_draws %>%
    filter(replicate > 0) %>%
    group_by(analysis_set, response, predictor, term, lambda, lambda_strategy) %>%
    group_modify(~ rubin_combine(.x$estimate, .x$std_error)) %>%
    ungroup() %>%
    left_join(
      ecology_draws %>%
        filter(replicate == 0) %>%
        select(
          analysis_set, response, predictor, term,
          point_estimate = estimate, point_std_error = std_error,
          point_p_value = p_value, point_n_taxa = n_taxa
        ),
      by = c("analysis_set", "response", "predictor", "term")
    )
  write_clean(ecology_rubin, "ecology_pgls_rubin_summary.csv")
}

manifest <- c(
  "Robust helix uncertainty propagation",
  paste("Created:", format(Sys.time(), tz = "UTC")),
  paste("Metrics:", normalizePath(metrics_path)),
  paste("Bootstrap draws:", normalizePath(draws_path)),
  paste("Bootstrap replicates:", length(unique(draws$bootstrap_replicate))),
  paste("Analysis set sizes:", paste(names(analysis_sets), lengths(analysis_sets), collapse = "; ")),
  "Primary set: helix RMS / fitted radius <= 0.10.",
  "Strict set: no provisional geometry-quality warnings.",
  "PGLS lambda is optimized on the point-estimate geometry and then held fixed across measurement-bootstrap draws.",
  "Uncertainty remains conditional on the traced semilandmarks and excludes manual placement repeatability."
)
writeLines(manifest, file.path(output_dir, "uncertainty_propagation_manifest.txt"))

cat("Completed robust-geometry uncertainty propagation.\n")
cat("Output:", output_dir, "\n")
