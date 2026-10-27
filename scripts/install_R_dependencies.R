# Install study-analysis dependencies into a separate library, preserving the
# user's existing packages. Run from any working directory with Rscript.
script_arg <- grep("^--file=", commandArgs(), value=TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash="/")
repo_root <- dirname(dirname(script_path))
args <- commandArgs(trailingOnly=TRUE)
lib <- if (length(args)) args[[1]] else file.path(repo_root, ".reviewer-library")
dir.create(lib, recursive=TRUE, showWarnings=FALSE)
lib <- normalizePath(lib, winslash="/")
.libPaths(c(lib, .Library))
Sys.setenv(RGL_USE_NULL="true")
options(repos=c(CRAN="https://cloud.r-project.org"), timeout=300)
t0 <- proc.time()[["elapsed"]]

declared <- trimws(readLines(file.path(repo_root, "r_package_inventory.txt")))
base_packages <- rownames(installed.packages(lib.loc=.Library))
pins <- c(RRPP="2.1.2", geomorph="4.0.10")
# rgl and jpeg are dependencies of the pinned historical geomorph release;
# they are absent from the dependency graph of the current CRAN version.
packages <- setdiff(unique(c(declared[nzchar(declared)], "rgl", "jpeg")), c(base_packages, names(pins)))
missing <- setdiff(packages, rownames(installed.packages(lib.loc=.libPaths())))
if (length(missing)) install.packages(missing, lib=lib, dependencies=NA)
for (package in names(pins)) {
  url <- sprintf("https://cran.r-project.org/src/contrib/Archive/%s/%s_%s.tar.gz",
                 package, package, pins[[package]])
  install.packages(url, lib=lib, repos=NULL, type="source")
}
needed <- setdiff(unique(c(declared, "rgl", "jpeg")), "")
status <- vapply(needed, requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))
if (!all(status)) stop("Dependencies unavailable: ", paste(names(status)[!status], collapse=", "))
for (package in names(pins)) {
  if (as.character(packageVersion(package)) != pins[[package]]) {
    stop("Required historical version was not installed: ", package)
  }
}
elapsed <- proc.time()[["elapsed"]] - t0
ip <- installed.packages(lib.loc=.libPaths())
write.csv(as.data.frame(ip[,c("Package", "Version", "Built"),drop=FALSE]),
          file.path(lib, "installed_versions.csv"), row.names=FALSE)
writeLines(capture.output(sessionInfo()), file.path(lib, "sessionInfo.txt"))
cat(sprintf("\nInstallation verified in %.1f seconds.\nLibrary: %s\n", elapsed, lib))
cat("For analysis runs, set R_LIBS_USER to this directory before starting Rscript.\n")
