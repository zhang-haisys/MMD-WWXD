# Copy to config.local.R and edit for your machine. config.local.R is ignored by Git.
project_root <- normalizePath(".", mustWork = FALSE)
data_dir <- file.path(project_root, "data")
results_dir <- file.path(project_root, "results")
helper_dir <- file.path(project_root, "R", "vendor")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
