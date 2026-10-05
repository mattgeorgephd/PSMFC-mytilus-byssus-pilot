## Batch-runner helpers shared by every 00_run_*.Rmd (and the repository-level runner).
## Base R plus rmarkdown.
##
## Each script is rendered in its own R process, exactly as knitting it by hand would be:
## nothing leaks between scripts, and a nested rmarkdown::render() cannot collide with the
## runner's own chunk labels. Every step leaves an HTML report and a log next to each other.

## Render one Rmd in a fresh R process. Returns TRUE when it succeeded.
render_step <- function(rmd, out_html, html_dir, params = NULL) {
  dir.create(html_dir, showWarnings = FALSE, recursive = TRUE)
  log_file <- file.path(html_dir, sub("[.]html$", ".log", out_html))
  job <- tempfile(fileext = ".R")
  prm <- if (length(params)) paste0("params = ", paste(deparse(params), collapse = ""), ", ") else ""
  writeLines(sprintf("rmarkdown::render(%s, %soutput_file = %s, output_dir = %s, quiet = TRUE)",
                     shQuote(normalizePath(rmd, winslash = "/")), prm, shQuote(out_html),
                     shQuote(normalizePath(html_dir, winslash = "/"))), job)
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  status <- system2(rscript, shQuote(job), stdout = log_file, stderr = log_file)
  unlink(job)
  status == 0
}

## The numbered steps of a folder, in run order: every 01_*.Rmd ... NN_*.Rmd in code_dir
## (00_* is the runner itself). `steps` is "all" or a vector of step numbers.
list_steps <- function(code_dir, steps = "all") {
  s <- sort(list.files(code_dir, pattern = "^[0-9]{2}_.*[.]Rmd$"), method = "radix")
  s <- s[!startsWith(s, "00_")]
  if (!identical(steps, "all")) {
    want <- sprintf("%02d", as.integer(steps))
    missing <- setdiff(want, substr(s, 1, 2))
    if (length(missing)) stop("no script numbered ", paste(missing, collapse = ", "), " in ", code_dir, call. = FALSE)
    s <- s[substr(s, 1, 2) %in% want]
  }
  s
}

## Run steps in order. `step_params` maps a step number ("03") to its knit parameters;
## `suffix` is appended to report names (e.g. "_F"). Returns the run log; the caller writes it
## and stops if any step failed, so a partial batch never looks complete.
run_steps <- function(code_dir, html_dir, steps = "all", step_params = list(),
                      stop_on_fail = TRUE, suffix = "") {
  scripts <- list_steps(code_dir, steps)
  log <- list()
  for (s in scripts) {
    stem <- tools::file_path_sans_ext(s)
    out  <- paste0(stem, suffix, ".html")
    message(sprintf("==== %s%s ====", s, if (nzchar(suffix)) paste0("  (", sub("^_", "", suffix), ")") else ""))
    t0 <- Sys.time()
    ok <- render_step(file.path(code_dir, s), out, html_dir, step_params[[substr(s, 1, 2)]])
    if (!ok) message("FAILED; see ", file.path(html_dir, sub("[.]html$", ".log", out)))
    log[[length(log) + 1]] <- data.frame(
      step = substr(s, 1, 2), script = s, run = sub("^_", "", suffix), ok = ok,
      minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
      html = if (ok) out else NA_character_, log = sub("[.]html$", ".log", out))
    if (!ok && isTRUE(stop_on_fail)) { message("Stopping (stop_on_fail = TRUE)."); break }
  }
  do.call(rbind, log)
}

## The error of a failed step, from its log: the lines from the first error line among the
## last 80, at most `n` of them (the last `n` lines if no line starts with "Error").
log_tail <- function(log_file, n = 40) {
  if (!file.exists(log_file)) return("(no log)")
  x <- readLines(log_file, warn = FALSE)
  x <- x[nzchar(trimws(x))]
  if (!length(x)) return("(empty log)")
  last <- utils::tail(seq_along(x), 80)
  from <- last[startsWith(x[last], "Error")][1]
  if (is.na(from)) return(utils::tail(x, n))
  x[from:min(length(x), from + n - 1)]
}

## Write the run log and fail the runner if any step failed. The error quotes each failed
## step's error from its log (`log` column, in html_dir), so the stage log and the knit
## console show the cause, not only where to look for it. Knitting captures message() output
## into the report, which is not written when the knit fails, so the error message carries it.
finish_run <- function(run_log, html_dir, what) {
  print(run_log)
  utils::write.csv(run_log, file.path(html_dir, "run_log.csv"), row.names = FALSE)
  if (!all(run_log$ok)) {
    failed <- run_log[!run_log$ok, , drop = FALSE]
    quoted <- unlist(lapply(failed$log, function(f)
      c("", paste0("---- ", f, " ----"), log_tail(file.path(html_dir, f)))))
    op <- options(warning.length = 8170)   # R's longest error message; the default cuts at 1000
    on.exit(options(op))
    stop(what, " step(s) failed: ", paste(failed$script, collapse = ", "),
         ". Logs: ", html_dir, "\n", paste(quoted, collapse = "\n"), call. = FALSE)
  }
  invisible(run_log)
}

## The packages the pipeline's scripts load (library(), require(), requireNamespace() or
## pkg::), from a scan of 02 to 09 and tools/ on 2026-10-05. Add a package here when a script
## starts using one; check_stack() stops when any of them is missing.
PIPELINE_PACKAGES <- c(
  "AnnotationDbi", "BiocManager", "Biostrings", "DESeq2", "DT", "EnvStats", "GO.db", "GOSemSim",
  "GSEABase", "GenomicAlignments", "GenomicRanges", "RColorBrewer", "apeglm", "ashr", "broom",
  "clusterProfiler", "colorspace", "data.table", "dplyr", "emmeans", "enrichplot", "forcats",
  "ggplot2", "ggvenn", "goseq", "gridExtra", "here", "httr", "jsonlite", "kableExtra", "knitr",
  "openxlsx", "org.Hs.eg.db", "patchwork", "pheatmap", "purrr", "readr", "readxl", "rentrez",
  "rmarkdown", "rrvgo", "rtracklayer", "scales", "stringr", "tibble", "tidyr", "tidyverse", "tm",
  "topGO", "tximport")

## Before a run: does this R match the one the analysis is recorded with (renv.lock)?
## Stops, listing every mismatch, when R's minor version or Bioconductor's differs from
## renv.lock, when a package the pipeline loads (`key`, PIPELINE_PACKAGES) is not installed in
## any library on R's path, or, with `need_go`, when GO.db does not hold the GO release of
## record (GO_RELEASE_OF_RECORD, from pipeline_checks.R); 07 and 08 stop on that anyway, but
## only after the stages before them have run. Other renv.lock packages that are missing, and
## packages at another version, are reported, not stopped on: the lock records the whole
## library the pipeline ran with, not only what it loads.
check_stack <- function(lockfile, need_go = TRUE, key = PIPELINE_PACKAGES) {
  txt <- paste(readLines(lockfile, warn = FALSE), collapse = "\n")
  section_version <- function(name) {
    rx <- sprintf('(?s).*?"%s"\\s*:\\s*\\{\\s*"Version"\\s*:\\s*"([^"]+)".*', name)
    if (grepl(rx, txt, perl = TRUE)) sub(rx, "\\1", txt, perl = TRUE)
    else stop(lockfile, " records no ", name, " version", call. = FALSE)
  }
  rec_r <- section_version("R")
  rec_bioc <- section_version("Bioconductor")
  pv <- regmatches(txt, gregexpr('"Package"\\s*:\\s*"[^"]+",\\s*"Version"\\s*:\\s*"[^"]+"', txt))[[1]]
  rec <- setNames(sub('.*"Version"\\s*:\\s*"([^"]+)"$', "\\1", pv), sub('^"Package"\\s*:\\s*"([^"]+)".*', "\\1", pv))
  minor <- function(v) paste(unlist(strsplit(as.character(v), "[.-]"))[1:2], collapse = ".")
  have <- function(p) nzchar(system.file(package = p))

  problems <- character(0)
  r_now <- paste(R.version$major, R.version$minor, sep = ".")
  if (minor(r_now) != minor(rec_r)) problems <- c(problems, sprintf("R %s; recorded R %s", r_now, rec_r))
  bioc_now <- if (have("BiocVersion")) minor(utils::packageVersion("BiocVersion")) else "not known (BiocVersion not installed)"
  if (have("BiocVersion") && bioc_now != rec_bioc)
    problems <- c(problems, sprintf("Bioconductor %s; recorded Bioconductor %s", bioc_now, rec_bioc))
  go_now <- "GO release not checked (no GO stage in this run)"
  if (need_go) {
    go_release <- if (have("GO.db")) tryCatch({
      info <- GO.db::GO_dbInfo()
      info$value[info$name == "GOSOURCEDATE"]
    }, error = function(e) NA_character_) else NA_character_
    go_now <- if (!have("GO.db")) "GO.db not installed" else if (is.na(go_release)) "GO.db does not load" else
      sprintf("GO.db %s, GO release %s", utils::packageVersion("GO.db"), go_release)
    if (!identical(go_release, GO_RELEASE_OF_RECORD))
      problems <- c(problems, sprintf("%s; 07 and 08 need GO release %s (GO.db %s)", go_now, GO_RELEASE_OF_RECORD, rec[["GO.db"]]))
  }
  missing <- names(rec)[!vapply(names(rec), have, logical(1))]
  present <- setdiff(names(rec), missing)
  differ  <- present[vapply(present, function(p) utils::packageVersion(p) != package_version(rec[[p]]), logical(1))]
  key_missing <- key[!vapply(key, have, logical(1))]
  if (length(key_missing))
    problems <- c(problems, sprintf(
      "%d of the %d packages the pipeline loads are not installed in any library on R's path: %s. R's library path (.libPaths()) is %s; if they are in another library (for example the one renv::restore() filled), add it with R_LIBS=<library> in .Renviron and restart R",
      length(key_missing), length(key), paste(key_missing, collapse = ", "),
      paste(.libPaths(), collapse = "; ")))

  message("R ", r_now, " (renv.lock: ", rec_r, "); Bioconductor ", bioc_now, " (renv.lock: ", rec_bioc, "); ", go_now)
  message(length(present), " of ", length(rec), " renv.lock packages installed, ", length(differ),
          " of them at another version",
          if (length(missing)) paste0("; missing: ", paste(utils::head(missing, 15), collapse = ", "),
                                      if (length(missing) > 15) ", ..." else "") else "")
  key_differ <- intersect(intersect(key, names(rec)), differ)
  if (length(key_differ))
    message("Pipeline packages not at the recorded version: ", paste(sprintf("%s %s (renv.lock %s)", key_differ,
      vapply(key_differ, function(p) as.character(utils::packageVersion(p)), ""), rec[key_differ]), collapse = ", "))
  if (length(problems)) {
    op <- options(warning.length = 8170)
    on.exit(options(op))
    stop("This R does not match the one the analysis is recorded with (", basename(lockfile), "):\n",
         paste0("  - ", problems, collapse = "\n"),
         "\nInstall the recorded stack (AGENTS.md, How to run): with R ", rec_r,
         ", install.packages(\"renv\") and renv::restore(lockfile = \"renv.lock\", library = \"<library>\", ",
         "prompt = FALSE), make <library> R's library (R_LIBS in .Renviron) and restart R. To run anyway, ",
         "knit with check_versions: false; the outputs are then not of record, and 07 and 08 still stop ",
         "unless GO.db holds GO release ", GO_RELEASE_OF_RECORD, ".", call. = FALSE)
  }
  invisible(TRUE)
}
