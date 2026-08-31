#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[startsWith(args, file_arg)][1])
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

r_files <- c(
  list.files(file.path(repo_root, "R"), pattern = "[.]R$", recursive = TRUE, full.names = TRUE),
  list.files(repo_root, pattern = "[.]R$", recursive = FALSE, full.names = TRUE),
  list.files(file.path(repo_root, "validation"), pattern = "[.]R$", recursive = TRUE, full.names = TRUE),
  list.files(file.path(repo_root, "tests"), pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
)
r_files <- sort(unique(normalizePath(r_files, mustWork = FALSE)))
r_files <- r_files[file.exists(r_files)]

if (!length(r_files)) {
  stop("No R scripts were found.")
}

errors <- character()
for (path in r_files) {
  result <- tryCatch(
    {
      parse(file = path, keep.source = TRUE)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  if (!is.null(result)) {
    errors <- c(errors, paste0(path, ": ", result))
  }
}

if (length(errors)) {
  cat("R SYNTAX CHECK: FAIL\n")
  cat(paste0("- ", errors), sep = "\n")
  cat("\n")
  quit(status = 1L)
}

cat("R SYNTAX CHECK: PASS\n")
cat("Scripts parsed:", length(r_files), "\n")
