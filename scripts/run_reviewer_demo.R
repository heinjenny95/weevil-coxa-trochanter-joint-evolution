# Rerun representative study analyses from shipped data and compare their
# numerical tables with v1.3.2. No original result files are overwritten.
script_arg <- grep("^--file=", commandArgs(), value=TRUE)
repo <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash="/")))
args <- commandArgs(trailingOnly=TRUE)
out <- if (length(args)) args[[1]] else file.path(repo, "reviewer-demo-output")
lib <- if (length(args) >= 2) args[[2]] else file.path(repo, ".reviewer-library")
if (!dir.exists(lib)) stop("Install dependencies first with scripts/install_R_dependencies.R.")
lib <- normalizePath(lib, winslash="/")
.libPaths(c(lib, .Library))
Sys.setenv(R_LIBS_USER=lib, R_LIBS="", R_LIBS_SITE="", RGL_USE_NULL="true")
pins <- c(RRPP="2.1.2", geomorph="4.0.10")
for (p in names(pins)) {
  if (as.character(packageVersion(p)) != pins[[p]]) stop("Use ", p, " ", pins[[p]], " for this reproduction check.")
}
dir.create(out, recursive=TRUE, showWarnings=FALSE)
out <- normalizePath(out, winslash="/")
if (out == repo || startsWith(repo, paste0(out, "/")) || startsWith(out, paste0(repo, "/data/"))) {
  stop("Choose an output directory separate from the released data.")
}
dir.create(file.path(out,"logs"),showWarnings=FALSE)
rscript <- file.path(R.home("bin"), "Rscript")
if (.Platform$OS.type == "windows") rscript <- paste0(rscript, ".exe")
data <- file.path(repo,"data")
src <- file.path(data,"supplementary_source_data")
pc <- file.path(src,"S01_PCA_and_Morphospace/PCA_scores_with_specimen_id.csv")
key <- file.path(data,"metadata/specimen_key.csv")
merged <- file.path(data,"screw_geometry/figure_source_data/allometry_merged_table.csv")
tree <- file.path(data,"phylogeny/P01_Trees/01_primary_tree_grafen.tre")
geometry <- file.path(data,"screw_geometry/robust_geometry_primary_adequate.csv")
timings <- list()
run <- function(name, script, arguments) {
  t0 <- proc.time()[["elapsed"]]
  status <- system2(rscript, shQuote(c("--vanilla", file.path(repo,"scripts/analysis",script),arguments)),
                    stdout=file.path(out,"logs",paste0(name,".log")), stderr=file.path(out,"logs",paste0(name,".log")))
  seconds <- proc.time()[["elapsed"]]-t0
  timings[[name]] <<- data.frame(workflow=name,seconds=seconds,exit_status=status)
  cat(sprintf("%s: %.2f s, exit %s\n", name, seconds, status))
  if (status != 0) stop("Workflow failed; inspect ", file.path(out,"logs",paste0(name,".log")))
}
t0 <- proc.time()[["elapsed"]]
run("proxy_mapping","revision_audit/build_proxy_mapping.R",c(key,file.path(out,"proxy_mapping.csv")))
run("shape_geometry","screw_geometry_analysis/screw_geometry_extraction_workflow.R",c(pc,geometry,file.path(out,"shape_geometry")))
run("full_shape_allometry","atlas_pca_workflow/r/test_full_shape_allometry.R",c(merged,file.path(out,"full_shape_allometry")))
run("phylogenetic_allometry","allometry_analysis/run_unified_phylogenetic_allometry.R",c(merged,key,tree,file.path(out,"phylogenetic_allometry")))
run("phylogenetic_multivariate_allometry","phylogenetic_comparative_analysis/run_phylogenetic_multivariate_allometry.R",c(merged,key,tree,file.path(out,"phylogenetic_multivariate_allometry")))
run("family_disparity","atlas_pca_workflow/r/calculate_family_disparity.R",c(pc,key,file.path(out,"family_disparity")))
run("coxa_thickness","coxa_wall_analysis/analyse_3d_coxa_wall_thickness.R",c(
  "--metrics",file.path(src,"S04_Coxa_Thickness/coxa_3d_wall_thickness_metrics.csv"),
  "--key",key,"--out-dir",file.path(out,"coxa_thickness"),"--figure-dir",file.path(out,"coxa_thickness/figures"),
  "--boundary-qc",file.path(data,"reviewer_support/coxa_mask_boundary_qc.csv")))

read_table <- function(path) {
  first <- readLines(path,n=1,warn=FALSE,encoding="UTF-8")
  read.table(path,header=TRUE,sep=if(grepl(";",first,fixed=TRUE)) ";" else ",",
             quote="\"",comment.char="",check.names=FALSE,colClasses="character",na.strings=c("NA","NaN",""))
}
checks <- list()
check <- function(actual,expected,filter_all_masks=FALSE) {
  x <- read_table(file.path(out,actual)); y <- read_table(file.path(repo,expected))
  if (filter_all_masks) x <- x[x$subset == "all_masks",,drop=FALSE]
  columns <- intersect(names(x),names(y))
  if (!nrow(y) || nrow(x)!=nrow(y) || !length(columns)) stop("Table dimensions differ: ",actual)
  bad <- character()
  for (column in columns) {
    xx <- x[[column]]; yy <- y[[column]]
    nx <- suppressWarnings(as.numeric(gsub(",",".",xx,fixed=TRUE)))
    ny <- suppressWarnings(as.numeric(gsub(",",".",yy,fixed=TRUE)))
    if (all(is.na(xx) | !is.na(nx)) && all(is.na(yy) | !is.na(ny))) {
      ok <- (is.na(nx)&is.na(ny)) | (!is.na(nx)&!is.na(ny)&abs(nx-ny)<=1e-10+1e-7*abs(ny))
    } else ok <- (is.na(xx)&is.na(yy)) | (!is.na(xx)&!is.na(yy)&xx==yy)
    if (!all(ok)) bad <- c(bad,column)
  }
  checks[[length(checks)+1]] <<- data.frame(output=actual,reference=expected,rows=nrow(x),columns_checked=length(columns),passed=!length(bad))
  if (length(bad)) stop("Reference mismatch in ",actual,": ",paste(bad,collapse=", "))
}
check("proxy_mapping.csv","data/metadata/taxonomic_proxy_mapping.csv")
check("shape_geometry/Table_regression_summary.csv","data/screw_geometry/regression_summary.csv")
check("shape_geometry/Table_plotting_dataset.csv","data/screw_geometry/shape_geometry_analysis_dataset.csv")
for (name in c("allometry_rrpp_multivariate_results.csv","allometry_procD_lm_results.csv",
               "allometry_rrpp_PC1_to_PC5_sensitivity_results.csv","allometry_univariate_all_PC_results.csv",
               "allometry_univariate_PC1_to_PC5_results.csv","allometry_PC_allometric_contributions.csv",
               "allometry_full_shape_regression_scores.csv","allometry_analysis_scope.csv")) {
  check(paste0("full_shape_allometry/",name),paste0("data/supplementary_source_data/S03_Allometry/",name))
}
for (name in c("pgls_allometry_all_traits.csv","pgls_allometry_all_traits_anova.csv",
               "pgls_allometry_tip_data_shape.csv","pgls_allometry_tip_data_geometry.csv",
               "pgls_results_main_traits.csv","pgls_geometry_main_traits_main_dataset.csv")) {
  check(paste0("phylogenetic_allometry/",name),paste0("data/supplementary_source_data/S03_Allometry/",name))
}
check("phylogenetic_multivariate_allometry/phylogenetic_multivariate_allometry_results.csv",
      "data/supplementary_source_data/S03_Allometry/phylogenetic_multivariate_allometry_results.csv")
for (name in c("family_disparity_summary_pc1_pc5.csv","family_disparity_rarefied_n4_pc1_pc5.csv",
               "family_disparity_global_test_pc1_pc5.csv","family_disparity_pairwise_tests_pc1_pc5.csv")) {
  check(paste0("family_disparity/",name),paste0("data/supplementary_source_data/S01_PCA_and_Morphospace/",name))
}
for (name in c("coxa_3d_wall_thickness_model_stats.csv","coxa_3d_wall_thickness_coefficients.csv",
               "coxa_3d_wall_thickness_sensitivity_tests.csv")) {
  check(paste0("coxa_thickness/",name),paste0("data/supplementary_source_data/S04_Coxa_Thickness/",name),TRUE)
}
report <- do.call(rbind,checks)
write.csv(report,file.path(out,"verification.csv"),row.names=FALSE)
write.csv(do.call(rbind,timings),file.path(out,"timings.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
cat(sprintf("PASS: %d reference tables verified in %.2f seconds.\n",nrow(report),proc.time()[["elapsed"]]-t0))
cat("Rerun outputs, logs and verification.csv are in: ",out,"\n",sep="")
