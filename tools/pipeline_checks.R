## Shared run checks and provenance for the analysis scripts.
## Sourced by 02_thread-strength scripts 3 and 4 and 09_gene-mechanics-correlation scripts 20
## and 21. Base R only (plus `tools`), so sourcing it adds no package dependency.
##
##   psmfc_repo_root(start)   the repository root, found by walking up from `start`
##   warn_unless(ok, ...)     a check that reports a failure but never stops the run
##   provenance_lines(settings, inputs, repo_root, packages)
##                            the lines of a RUN_provenance*.txt file: the caller's settings,
##                            the code commit and any uncommitted tracked changes, R and
##                            package versions, the checks run in this session, and an MD5
##                            of every input, with repository-relative paths

## Folders whose presence marks the repository root. Keep in step with the find_repo_root()
## helpers at the top of the 09_gene-mechanics-correlation scripts, which cannot source this
## file until they have found the root.
PSMFC_ROOT_MARKERS <- c("02_thread-strength", "06_differential-expression")

psmfc_repo_root <- function(start = getwd()) {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  repeat {
    if (all(dir.exists(file.path(d, PSMFC_ROOT_MARKERS)))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Repo root not found above ", start, ": need a folder containing ",
       paste0(PSMFC_ROOT_MARKERS, "/", collapse = " and "), ".", call. = FALSE)
}

## ---- checks -------------------------------------------------------------------------
## One tally per R process. Each script is knit in its own process and sources this file
## once, so sourcing resets the tally.
.psmfc_checks <- new.env(parent = emptyenv())
.psmfc_checks$run    <- 0L
.psmfc_checks$failed <- character(0)

## A failed check prints `CHECK FAILED: <message>` and the run carries on. It is written with
## cat() so it lands in the knitted report whatever the chunk's `message` / `warning`
## options are, and to stderr() so it also lands in the render log the batch drivers keep.
warn_unless <- function(ok, ...) {
  .psmfc_checks$run <- .psmfc_checks$run + 1L
  if (isTRUE(ok)) return(invisible(TRUE))
  msg  <- paste0(...)
  line <- paste0("CHECK FAILED: ", msg, "\n")
  .psmfc_checks$failed <- c(.psmfc_checks$failed, msg)
  cat(line)
  cat(line, file = stderr())
  invisible(FALSE)
}

## ---- provenance ---------------------------------------------------------------------
provenance_lines <- function(settings, inputs, repo_root, packages = character(0)) {
  root <- normalizePath(repo_root, winslash = "/", mustWork = FALSE)
  field <- function(label, value) sprintf("%-21s: %s", label, value)

  ## Repository-relative path. A missing file is resolved through its parent folder when
  ## that exists, so `..` segments do not survive into the record.
  rel <- function(p) {
    p <- if (file.exists(p)) normalizePath(p, winslash = "/")
         else if (dir.exists(dirname(p))) file.path(normalizePath(dirname(p), winslash = "/"), basename(p))
         else p
    if (startsWith(p, paste0(root, "/"))) substring(p, nchar(root) + 2L) else p
  }

  ## git state. Returns NULL when git is not installed or this is not a checkout (for
  ## example a downloaded zip), so provenance still writes.
  git <- function(...) {
    out <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(root), ...),
                                             stdout = TRUE, stderr = FALSE)),
                    error = function(e) NULL)
    if (is.null(out) || !is.null(attr(out, "status"))) NULL else out
  }
  commit <- git("rev-parse", "--short=12", "HEAD")
  git_lines <- if (length(commit) == 1L) {
    st <- git("status", "--porcelain", "--untracked-files=no")
    st <- sub("^\\s*(\\S+)\\s+", "\\1 ", st)      # "XY path" -> "X path"
    c(field("code commit", commit),
      field("uncommitted changes",
            if (length(st) == 0L) "NO"
            else paste0("YES, ", length(st), " tracked file(s) differ from the commit: ",
                        paste(utils::head(st, 10L), collapse = "; "),
                        if (length(st) > 10L) " ..." else "")))
  } else {
    field("code commit", "unavailable (not a git checkout, or git is not on the PATH)")
  }

  pkg_version <- function(p) tryCatch(as.character(utils::packageVersion(p)),
                                      error = function(e) "not installed")
  pkg_lines <- if (length(packages))
    field("packages", paste(sprintf("%s %s", packages, vapply(packages, pkg_version, "")),
                            collapse = ", "))

  check_lines <- if (.psmfc_checks$run > 0L)
    field("checks", paste0(.psmfc_checks$run, " run, ", length(.psmfc_checks$failed), " failed",
                           if (length(.psmfc_checks$failed))
                             paste0(": ", paste(.psmfc_checks$failed, collapse = " | "))))

  inputs <- unique(as.character(inputs))
  is_file <- file_test("-f", inputs)
  md5 <- rep("MISSING", length(inputs))
  md5[is_file] <- unname(tools::md5sum(inputs[is_file]))
  input_lines <- sprintf("  %s  %s", md5, vapply(inputs, rel, "", USE.NAMES = FALSE))

  c(settings, git_lines,
    field("R", R.version.string),
    pkg_lines, check_lines,
    "inputs (md5, repo-relative path):", input_lines)
}
