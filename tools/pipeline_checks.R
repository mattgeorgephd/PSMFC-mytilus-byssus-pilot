## Shared pipeline checks and provenance, sourced by the thread-strength and gene-mechanics
## scripts. Base R only, so it adds no package dependency.

## A check that warns, naming what failed, when `ok` is not TRUE. The script continues.
## The scripts knit with `warning = FALSE`, which drops warnings from both the report and the
## render log, so a failure is also printed to the report (stdout) and to the log (stderr),
## and counted for provenance_lines().
.checks <- new.env()
.checks$n_run  <- 0L
.checks$failed <- character(0)
warn_unless <- function(ok, ...) {
  .checks$n_run <- .checks$n_run + 1L
  if (!isTRUE(ok)) {
    msg <- paste0(...)
    .checks$failed <- c(.checks$failed, msg)
    cat("CHECK FAILED: ", msg, "\n", sep = "")
    cat("CHECK FAILED: ", msg, "\n", sep = "", file = stderr())
    warning(msg, call. = FALSE)
  }
  invisible(isTRUE(ok))
}

## Repository root: the first parent that holds both analysis folders the pipelines join.
psmfc_repo_root <- function(start = getwd()) {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in 1:8) {
    if (all(dir.exists(file.path(d, c("05_differential-expression", "02_thread-strength"))))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Repository root not found: need a parent folder containing 05_differential-expression/ ",
       "and 02_thread-strength/.", call. = FALSE)
}

## Path relative to the repository root, with forward slashes, for provenance files.
rel_path <- function(path, repo_root) {
  p <- normalizePath(path, winslash = "/", mustWork = FALSE)
  r <- paste0(normalizePath(repo_root, winslash = "/", mustWork = FALSE), "/")
  if (startsWith(p, r)) substring(p, nchar(r) + 1) else p
}

.git_out <- function(repo_root, args) {
  if (!nzchar(Sys.which("git"))) return(NA_character_)
  out <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(repo_root), args),
                                           stdout = TRUE, stderr = FALSE)),
                  error = function(e) character(0))
  if (length(out) == 0 || !is.null(attr(out, "status"))) NA_character_ else out
}

## Code commit that ran, and whether any tracked code or input differs from it. Output folders
## (03_analyses) and knitted reports are excluded, since a run rewrites them.
git_state <- function(repo_root) {
  head <- .git_out(repo_root, c("rev-parse", "--short=12", "HEAD"))
  if (is.na(head[1])) return(c(commit = "unavailable (git not found)", dirty = "unknown"))
  changed <- .git_out(repo_root, c("status", "--porcelain", "--untracked-files=no"))
  changed <- if (all(is.na(changed))) character(0) else changed
  changed <- changed[!grepl("03_analyses/|\\.html$", changed)]
  c(commit = head[1],
    dirty  = if (length(changed) == 0) "no" else
      paste0("YES, ", length(changed), " tracked file(s) differ from the commit: ",
             paste(trimws(utils::head(changed, 10)), collapse = "; "),
             if (length(changed) > 10) " ..." else ""))
}

## Provenance block: code commit and dirty state, one MD5 line per input file (repo-relative),
## then the R session. `header` is a character vector of setting lines written first.
provenance_lines <- function(header, inputs, repo_root, packages = character(0)) {
  gs <- git_state(repo_root)
  inputs <- unique(inputs[!is.na(inputs)])
  md5 <- vapply(inputs, function(f)
    if (file.exists(f) && !dir.exists(f)) unname(tools::md5sum(f)) else "MISSING", character(1))
  pkg <- vapply(packages, function(p)
    if (requireNamespace(p, quietly = TRUE)) paste0(p, " ", as.character(utils::packageVersion(p)))
    else paste0(p, " (not installed)"), character(1))
  c(header,
    if (.checks$n_run > 0)
      c(paste0("checks               : ", .checks$n_run, " run, ", length(.checks$failed), " failed"),
        if (length(.checks$failed)) paste0("  FAILED: ", .checks$failed)),
    paste0("code commit          : ", gs[["commit"]]),
    paste0("uncommitted changes  : ", gs[["dirty"]]),
    paste0("R                    : ", R.version.string),
    if (length(pkg)) paste0("packages             : ", paste(pkg, collapse = ", ")),
    "inputs (md5, repo-relative path):",
    paste0("  ", md5, "  ", vapply(inputs, rel_path, character(1), repo_root = repo_root)))
}
