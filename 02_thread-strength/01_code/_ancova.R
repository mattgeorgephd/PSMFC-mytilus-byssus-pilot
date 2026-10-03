## Per-animal ANCOVA: the thread-strength analysis of record. Sourced by scripts 04 and 05.
##
## Arms were assigned at random, so each metric is modelled as the animal's day-3 level
## adjusted for its own baseline level:
##
##     y_day3 ~ arm + y_baseline        (lm, one row per animal)
##
## on the log scale (every metric is positive and right-skewed), so arm effects read as
## ratios. THREAD_METRICS below defines each metric from a thread-level column of the thread
## summary and how it is summarised per animal and timepoint:
##
##   adhesion_kpa  mean of the log thread adhesions (log of the geometric mean)
##   mean_force    mean of the log thread peak forces (`max_force` of each trace): the
##                 animal's typical thread
##   max_force     the largest thread peak force: the animal's strongest thread
##   pad_area      mean of the log plaque areas
##
## Extension at break is not analysed: the thread was cut near the plaque-distal junction, so
## the length of distal thread under test differed between pulls and extension cannot be
## measured. Only animals with threads at both timepoints enter; the lab-reference animals
## have no day-3 pull and never do. The day-3 treatment control is the reference arm.
##
## Reported per metric:
##   tests           F for arm and for baseline, each adjusted for the other (drop1)
##   coefficients    lm coefficients with 95% CI
##   adjusted_means  arm means at the mean baseline (back-transformed for log metrics)
##   vs_control      OA, OW, DO each vs control: ratio (log) or difference (raw), with 95% CI
##                   and p both unadjusted and Dunnett-adjusted over the three comparisons
##                   (emmeans `adjust = "mvt"`, the exact multivariate-t method; seeded)
##   data            the per-animal rows the model was fitted to

ANCOVA_ARMS  <- c("control", "OA", "OW", "DO")
DUNNETT_SEED <- 20260930

## metric -> thread-level column, per-animal summary, model scale, figure label
THREAD_METRICS <- data.frame(
  metric = c("adhesion_kpa", "mean_force", "max_force", "pad_area"),
  column = c("adhesion_kpa", "max_force",  "max_force", "pad_area"),
  agg    = c("mean",         "mean",       "max",       "mean"),
  scale  = c("log",          "log",        "log",       "log"),
  label  = c("Adhesion", "Mean peak force", "Maximum peak force", "Plaque area"),
  stringsAsFactors = FALSE)
metric_spec <- function(metric) {
  i <- match(metric, THREAD_METRICS$metric)
  if (is.na(i)) stop("Unknown metric '", metric, "'; see THREAD_METRICS in _ancova.R.")
  as.list(THREAD_METRICS[i, ])
}

ancova_data <- function(threads, metric, scale = metric_spec(metric)$scale,
                        column = metric_spec(metric)$column, agg = metric_spec(metric)$agg) {
  stopifnot(scale %in% c("log", "raw"), agg %in% c("mean", "max"), column %in% names(threads))
  d <- threads %>%
    filter(phase %in% c("pre", "post"), as.character(mussel_trt) %in% ANCOVA_ARMS,
           !is.na(.data[[column]]))
  if (scale == "log" && any(d[[column]] <= 0))
    stop(column, " has non-positive values; it cannot be modelled on the log scale.")
  d %>%
    mutate(v = if (scale == "log") log(.data[[column]]) else .data[[column]],
           mussel = as.character(mussel), arm = as.character(mussel_trt)) %>%
    group_by(mussel, arm, phase) %>%
    summarise(value = if (agg == "max") max(v) else mean(v), n_threads = n(), .groups = "drop") %>%
    pivot_wider(names_from = phase, values_from = c(value, n_threads)) %>%
    filter(!is.na(value_pre), !is.na(value_post)) %>%
    transmute(mussel, arm = factor(arm, levels = ANCOVA_ARMS),
              y = value_post, baseline = value_pre,
              n_threads_pre, n_threads_post) %>%
    arrange(arm, mussel)
}

fit_ancova <- function(threads, metric, scale = metric_spec(metric)$scale) {
  is_log <- scale == "log"   # a plain flag: inside mutate()/tibble() `scale` would mean the column
  d <- ancova_data(threads, metric, scale)
  missing_arms <- setdiff(ANCOVA_ARMS, as.character(unique(d$arm)))
  if (length(missing_arms) > 0)
    stop("No paired animals in arm(s) ", paste(missing_arms, collapse = ", "), " for ", metric,
         "; the ANCOVA needs every arm.")
  fit <- lm(y ~ arm + baseline, data = d)

  tests <- as.data.frame(drop1(fit, test = "F")) %>%
    tibble::rownames_to_column("term") %>%
    filter(term != "<none>") %>%
    transmute(metric = metric, scale = scale, term, df = Df, sum_sq = `Sum of Sq`,
              F = `F value`, p = `Pr(>F)`, df_residual = fit$df.residual, n_animals = nrow(d))

  coefficients <- broom::tidy(fit, conf.int = TRUE) %>%
    mutate(metric = metric, scale = scale, .before = 1)

  # Arm means at the mean baseline. The response was logged before fitting, so declare the
  # transformation; without it `type = "response"` would silently stay on the log scale.
  emm <- emmeans(fit, ~ arm)
  if (is_log) emm <- update(emm, tran = "log")

  n_arm <- d %>% count(arm, name = "n_animals")
  baseline_at <- if (is_log) exp(mean(d$baseline)) else mean(d$baseline)
  adjusted_means <- as.data.frame(summary(emm, type = "response")) %>%
    rename(adjusted_mean = any_of(c("response", "emmean"))) %>%
    mutate(arm = factor(as.character(arm), levels = ANCOVA_ARMS)) %>%
    left_join(n_arm, by = "arm") %>%
    mutate(metric = metric, scale = scale, baseline_at = baseline_at, .before = 1)

  con <- contrast(emm, method = "trt.vs.ctrl", ref = 1)
  est_col <- if (is_log) "ratio" else "estimate"
  un <- as.data.frame(summary(con, infer = c(TRUE, TRUE), adjust = "none", type = "response"))
  set.seed(DUNNETT_SEED)
  du <- as.data.frame(summary(con, infer = c(TRUE, TRUE), adjust = "mvt", type = "response"))
  stopifnot(identical(as.character(un$contrast), as.character(du$contrast)))

  vs_control <- tibble(
    metric = metric, scale = scale,
    contrast = as.character(un$contrast),
    arm = sub(" .*$", "", as.character(un$contrast)),
    estimate_type = if (is_log) "ratio arm / control" else "difference arm - control",
    estimate = un[[est_col]],
    pct_change = if (is_log) 100 * (un[[est_col]] - 1) else NA_real_,
    SE = un$SE, df = un$df, t = un$t.ratio,
    ci_low = un$lower.CL, ci_high = un$upper.CL, p_unadjusted = un$p.value,
    ci_low_dunnett = du$lower.CL, ci_high_dunnett = du$upper.CL, p_dunnett = du$p.value)

  list(metric = metric, scale = scale, data = d, fit = fit, tests = tests,
       coefficients = coefficients, adjusted_means = adjusted_means, vs_control = vs_control)
}

write_ancova_tables <- function(results, out_dir) {
  pick <- function(part) bind_rows(lapply(results, `[[`, part))
  write.csv(pick("tests"),          file.path(out_dir, "STATS_ancova_tests.csv"),          row.names = FALSE)
  write.csv(pick("coefficients"),   file.path(out_dir, "STATS_ancova_coefficients.csv"),   row.names = FALSE)
  write.csv(pick("adjusted_means"), file.path(out_dir, "STATS_ancova_adjusted_means.csv"), row.names = FALSE)
  write.csv(pick("vs_control"),     file.path(out_dir, "STATS_ancova_vs_control.csv"),     row.names = FALSE)
  write.csv(bind_rows(lapply(results, function(r)
              r$data %>% mutate(metric = r$metric, scale = r$scale, .before = 1))),
            file.path(out_dir, "DATA_ancova_animals.csv"), row.names = FALSE)
}

save_ancova_diagnostics <- function(results, file) {
  png(file, width = 1600, height = 700 * length(results), res = 150)
  on.exit(dev.off())
  par(mfrow = c(length(results), 2))
  for (r in results) {
    lab <- paste0(r$metric, if (r$scale == "log") " (log)" else "")
    qqnorm(resid(r$fit), main = paste("Normal Q-Q,", lab)); qqline(resid(r$fit))
    plot(fitted(r$fit), resid(r$fit), xlab = "Fitted", ylab = "Residuals",
         main = paste("Residuals vs fitted,", lab)); abline(h = 0, lty = 2)
  }
}
