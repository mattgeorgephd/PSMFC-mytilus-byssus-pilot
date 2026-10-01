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
      html = if (ok) out else NA_character_)
    if (!ok && isTRUE(stop_on_fail)) { message("Stopping (stop_on_fail = TRUE)."); break }
  }
  do.call(rbind, log)
}

## Write the run log and fail the runner if any step failed.
finish_run <- function(run_log, html_dir, what) {
  print(run_log)
  utils::write.csv(run_log, file.path(html_dir, "run_log.csv"), row.names = FALSE)
  if (!all(run_log$ok))
    stop(what, " step(s) failed: ", paste(run_log$script[!run_log$ok], collapse = ", "),
         ". See the .log files in ", html_dir, ".", call. = FALSE)
  invisible(run_log)
}
