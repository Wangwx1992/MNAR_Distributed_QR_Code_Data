# NVSS 2024 natality analyses accompanying the paper.
# From the package root:
#   Rscript application/run_application.R --mode fixed
#   Rscript application/plot_application.R
#   Rscript application/run_application.R --mode repeated --workers 4
#   Rscript application/run_application.R --mode natural
# Dependencies and solver settings are defined in simulation/quantile.R
# and simulation/calibration.R. Use a one-thread BLAS for the timed experiment.
# Multicore repeated masking uses fork and requires macOS or Linux.
#
# The analysis input fields come from the NCHS 2024 Natality Public Use File:
# https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/DVS/natality/Nat2024us.zip
# Data are in native RDS format and already include the paper's fixed cohort.

get_script_directory <- function() {
  invoked <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(invoked)) return(dirname(normalizePath(sub("^--file=", "", invoked[[1L]]))))
  stop("Run with Rscript or set the application directory before sourcing.", call. = FALSE)
}

parse_app_args <- function(args, defaults) {
  if (any(args == "--help")) {
    cat(paste0("Usage: Rscript application/run_application.R",
               " --mode fixed|repeated|natural [--out DIRECTORY]\n",
               "Repeated masking: [--reps 100] [--workers 1] [--bootstrap 5000]\n",
               "Natural missingness: [--partitions 20] [--primary-only]\n",
               "Smoke test all modes with --limit ROWS.\n"))
    return(NULL)
  }
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (!startsWith(args[[i]], "--") || !key %in% names(defaults))
      stop("Unknown option: ", args[[i]], call. = FALSE)
    if (is.logical(defaults[[key]])) {
      defaults[[key]] <- TRUE
      i <- i + 1L
    } else {
      if (i == length(args) || startsWith(args[[i + 1L]], "--"))
        stop("Missing value for ", args[[i]], call. = FALSE)
      value <- args[[i + 1L]]
      if (is.integer(defaults[[key]])) {
        value <- suppressWarnings(as.integer(value))
        if (is.na(value) || value < 1L)
          stop("Expected a positive integer for --", key, call. = FALSE)
      }
      defaults[[key]] <- value
      i <- i + 2L
    }
  }
  defaults
}

coefficient_names <- c("intercept", "age_decade", "age_decade_squared", "height_5in",
                       paste0("education_", c(1, 2, 4:8)))

read_fixed_cohort <- function(root, limit = .Machine$integer.max) {
  base <- readRDS(file.path(root, "data", "fixed_cohort.rds"))
  needed <- c("x", "v", "y", "z")
  stopifnot(all(needed %in% names(base)), length(base$y) == 500000L,
            ncol(base$x) == 11L, ncol(base$v) == 5L)
  if (limit >= length(base$y)) return(base)
  rows <- seq_len(limit)
  base$x <- base$x[rows, , drop = FALSE]
  base$v <- base$v[rows, , drop = FALSE]
  base$y <- base$y[rows]
  base$z <- base$z[rows]
  if (!is.null(base$fixed_mask)) base$fixed_mask <- as.logical(base$fixed_mask[rows])
  base
}

load_natural_cohort <- function(root, limit) {
  raw <- readRDS(file.path(root, "data", "natality_fields.rds"))
  take <- raw[, "RESTATUS"] %in% 1:3 & raw[, "DPLURAL"] == 1L &
    raw[, "DWgt_R"] >= 100L & raw[, "DWgt_R"] <= 400L &
    raw[, "M_Ht_In"] >= 30L & raw[, "M_Ht_In"] <= 78L &
    raw[, "MEDUC"] %in% 1:8 & raw[, "F_PWGT"] == 1L & raw[, "F_DWGT"] == 1L
  eligible <- which(take)
  a <- raw[head(eligible, limit), , drop = FALSE]
  rm(raw, take)
  response <- a[, "PWgt_R"] != 999L
  age <- (a[, "MAGER"] - 30) / 10
  height <- (a[, "M_Ht_In"] - 65) / 5
  education <- a[, "MEDUC"]
  z <- tanh((a[, "DWgt_R"] - 190) / 40)
  y <- rep(NA_real_, length(response))
  y[response] <- (a[response, "PWgt_R"] - 160) / 40
  indicators <- vapply(c(1L, 2L, 4:8), function(j) as.double(education == j),
                       numeric(length(response)))
  x <- cbind(1, age, age^2, height, indicators)
  colnames(x) <- coefficient_names
  v <- cbind(1, age, age^2, height, (education - 4) / 2)
  stratum <- function(kind, labels) {
    categories <- sort(unique(labels))
    n <- as.numeric(table(factor(labels, levels = categories)))
    missing <- as.numeric(table(factor(labels[!response], levels = categories)))
    data.frame(stratum = kind, category = categories, n = n, missing_y = missing,
               missing_fraction = missing / n)
  }
  list(x = x, v = v, y = y, z = z, r = response,
       descriptives = rbind(stratum("education", education),
                            stratum("age_band", pmin(a[, "MAGER"] %/% 5L * 5L, 45L))),
       size = data.frame(n = length(response), observed_y = sum(response),
                         missing_y = sum(!response), age_imputed = sum(a[, "MAGE_IMPFLG"] == 1L),
                         race_imputed = sum(a[, "MRACEIMP"] == 1L)))
}

selection_inputs <- function(base, r, mode = "MNAR") {
  if (mode == "MAR") return(list(t = base$v, h = base$v))
  transformed <- rep(NA_real_, length(r))
  transformed[r] <- tanh(base$y[r])
  list(t = cbind(base$v, transformed), h = cbind(base$v, base$z, base$z^2))
}

mask_probabilities <- function(base, missing, gamma) {
  linear <- .3 * base$v[, 2L] - .3 * base$v[, 4L] + .15 * base$v[, 5L] +
    gamma * tanh(base$y)
  intercept <- uniroot(function(a) mean(plogis(a + linear)) - (1 - missing),
                       c(-30, 30), tol = 1e-13)$root
  plogis(intercept + linear)
}

weight_stats <- function(weights) {
  c(min_respondent_probability = 1 / max(weights),
    max_respondent_weight = max(weights),
    effective_sample_size = sum(weights)^2 / sum(weights^2))
}

pounds <- function(beta) {
  value <- 40 * beta
  value[[1L]] <- value[[1L]] + 160
  value
}

fit_app_selection <- function(base, r, partitions = 1L, mode = "MNAR", delta = 0,
                              warm = NULL) {
  design <- selection_inputs(base, r, mode)
  fit_mnar(design$t, design$h, r, partitions = partitions,
           bounds = c(10, rep(5, ncol(design$t) - 1L)),
           tolerance = 1 / length(r), starts = if (is.null(warm)) NULL else list(warm),
           offset = delta * base$z, initialization = "multistart",
           outcome_dependent = mode == "MNAR")
}

require_success <- function(fit, description) {
  if (!isTRUE(fit$converged) || any(!is.finite(fit$beta)))
    stop("Quantile optimization did not converge: ", description, call. = FALSE)
  fit
}

run_fixed_mask <- function(base, out, partitions = 20L) {
  n <- length(base$y)
  if (is.null(base$fixed_mask))
    stop("The original fixed response mask is missing from the cohort input.", call. = FALSE)
  r <- base$fixed_mask
  if (n == 500000L && sum(r) != 299913L)
    stop("Unexpected respondent count in the archived fixed mask.", call. = FALSE)
  if (sum(r) < 200L || all(r)) stop("Fixed mask needs at least 200 respondents and some missingness.")
  x <- base$x[r, , drop = FALSE]
  y <- base$y[r]
  before <- proc.time()[["elapsed"]]
  selection <- fit_app_selection(base, r, partitions)
  calibration_seconds <- proc.time()[["elapsed"]] - before
  if (!isTRUE(selection$converged)) stop("Fixed-mask calibration did not converge.", call. = FALSE)
  weights <- selection$weights
  metrics <- c(n = n, n_observed = sum(r), missing_fraction = mean(!r),
               outcome_dependence = tail(selection$theta, 1L), calibration_seconds = calibration_seconds,
               overid_statistic = selection$overid_statistic,
               overid_df = selection$overid_df, overid_pvalue = selection$overid_pvalue,
               whitened_condition = selection$whitened_condition,
               weight_stats(weights))
  write.csv(as.data.frame(as.list(metrics)), file.path(out, "calibration.csv"), row.names = FALSE)

  settings <- rbind(data.frame(method = "IP", M = 0L),
                    expand.grid(method = c("ADMM", "Renewable"), M = c(20L, 50L, 100L),
                                stringsAsFactors = FALSE))
  taus <- seq(.1, .9, by = .1)
  rows <- vector("list", nrow(settings) * length(taus))
  coefficients <- rows
  count <- 0L
  for (tau in taus) {
    for (j in seq_len(nrow(settings))) {
      method <- settings$method[[j]]
      blocks <- settings$M[[j]]
      before <- proc.time()[["elapsed"]]
      fit <- switch(method,
        IP = fit_qr_central(x, y, weights, tau, n_total = n),
        ADMM = fit_qr_admm(x, y, weights, tau, partitions = blocks, n_total = n,
                          tolerance = 1e-4, stopping = "residual", maxit = 20000L),
        Renewable = fit_qr_renewable(x, y, weights, tau, partitions = blocks, n_total = n,
                                    tolerance = 1e-4, maxit = 5000L))
      elapsed <- proc.time()[["elapsed"]] - before
      if (method == "IP") reference <- fit
      count <- count + 1L
      plotted_qr <- if (method == "ADMM") fit$parallel_seconds else elapsed
      if (length(plotted_qr) != 1L || !is.finite(plotted_qr) || plotted_qr < 0)
        stop("Invalid quantile elapsed time for ", method, call. = FALSE)
      rows[[count]] <- data.frame(method = method, M = blocks, tau = tau,
        height_coefficient_lb = 40 * fit$beta[[4L]],
        max_coefficient_difference_from_IP_lb = 40 * max(abs(fit$beta - reference$beta)),
        weighted_loss = fit$objective, loss_excess_from_IP = fit$objective - reference$objective,
        beta_norm = sqrt(sum(fit$beta^2)), converged = isTRUE(fit$converged),
        calibration_seconds = calibration_seconds, qr_wall_seconds = elapsed,
        qr_figure_seconds = plotted_qr, total_wall_seconds = calibration_seconds + elapsed,
        total_figure_seconds = calibration_seconds + plotted_qr)
      coefficients[[count]] <- data.frame(method = method, M = blocks, tau = tau,
                                          coefficient = coefficient_names, estimate_lb = pounds(fit$beta))
      message(sprintf("%s M=%s tau=%.1f converged=%s", method, blocks, tau, isTRUE(fit$converged)))
    }
  }
  summary <- do.call(rbind, rows)
  write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
  write.csv(do.call(rbind, coefficients), file.path(out, "coefficients.csv"), row.names = FALSE)
  if (any(!summary$converged)) warning("The fixed-mask summary includes fits that did not meet their stopping rules.")
  invisible(summary)
}

smoothed_empirical_loss <- function(x, y, beta, tau, a = .0125) {
  residual <- y - drop(x %*% beta)
  mean(residual * (tau - pnorm(-residual / a)) + a * dnorm(residual / a))
}

summarize_recovery <- function(metrics, bootstrap_reps, reps, out) {
  set.seed(2026119102L)
  resamples <- matrix(sample.int(reps, bootstrap_reps * reps, replace = TRUE), nrow = bootstrap_reps)
  bootstrap_mean <- function(x) rowMeans(matrix(x[resamples], nrow = bootstrap_reps))
  interval <- function(x) unname(quantile(x, c(.025, .975)))
  groups <- split(metrics, list(metrics$missing_fraction, metrics$gamma, metrics$tau), drop = TRUE)
  summary <- list()
  paired <- list()
  for (group in groups) {
    index <- group[1L, c("missing_fraction", "gamma", "tau")]
    key <- function(method) {
      part <- group[group$method == method, , drop = FALSE]
      part <- part[order(part$rep), , drop = FALSE]
      stopifnot(identical(as.integer(part$rep), seq_len(reps)))
      part
    }
    target <- key("MNAR")
    target_boot <- sqrt(bootstrap_mean(target$prediction_mse_pounds2))
    for (method in c("CC", "MAR", "MNAR", "Oracle")) {
      part <- key(method)
      boot <- sqrt(bootstrap_mean(part$prediction_mse_pounds2))
      pred_int <- interval(boot)
      risk_int <- interval(bootstrap_mean(part$excess_risk_pounds))
      rms <- sqrt(mean(part$prediction_mse_pounds2))
      summary[[length(summary) + 1L]] <- cbind(index, data.frame(method = method, repetitions = reps,
        prediction_rmse_pounds = rms, prediction_rmse_mc_low = pred_int[[1L]],
        prediction_rmse_mc_high = pred_int[[2L]], mean_excess_risk_pounds = mean(part$excess_risk_pounds),
        mean_excess_risk_mc_low = risk_int[[1L]], mean_excess_risk_mc_high = risk_int[[2L]],
        coefficient_rmse_pounds = sqrt(mean(part$coefficient_mse_pounds2))))
      if (method != "MNAR") {
        ratio <- sqrt(mean(target$prediction_mse_pounds2)) / rms
        ratio_int <- interval(target_boot / boot)
        mse_diff <- target$prediction_mse_pounds2 - part$prediction_mse_pounds2
        risk_diff <- target$excess_risk_pounds - part$excess_risk_pounds
        mse_int <- interval(bootstrap_mean(mse_diff))
        diff_int <- interval(bootstrap_mean(risk_diff))
        paired[[length(paired) + 1L]] <- cbind(index, data.frame(comparator = method, repetitions = reps,
          prediction_rmse_ratio = ratio, prediction_rmse_ratio_mc_low = ratio_int[[1L]],
          prediction_rmse_ratio_mc_high = ratio_int[[2L]], prediction_rmse_reduction_percent = 100 * (1 - ratio),
          mean_paired_prediction_mse_difference = mean(mse_diff),
          paired_prediction_mse_difference_mc_low = mse_int[[1L]],
          paired_prediction_mse_difference_mc_high = mse_int[[2L]],
          mean_paired_excess_risk_difference = mean(risk_diff),
          paired_excess_risk_difference_mc_low = diff_int[[1L]],
          paired_excess_risk_difference_mc_high = diff_int[[2L]]))
      }
    }
  }
  write.csv(do.call(rbind, summary), file.path(out, "summary_metrics.csv"), row.names = FALSE)
  write.csv(do.call(rbind, paired), file.path(out, "paired_comparisons.csv"), row.names = FALSE)
}

respondent_start <- function(x, y, tau, ols) {
  start <- ols
  start[[1L]] <- start[[1L]] + unname(quantile(y - drop(x %*% ols), tau))
  start
}

run_repeated_masks <- function(base, out, reps, workers, bootstrap_reps) {
  n <- length(base$y)
  tau_values <- c(.1, .5, .9)
  cells <- expand.grid(missing_fraction = c(.2, .4), gamma = c(-1.5, 0, 1.5))
  probability <- lapply(seq_len(nrow(cells)), function(j)
    mask_probabilities(base, cells$missing_fraction[[j]], cells$gamma[[j]]))
  x <- base$x
  y <- base$y
  gram <- crossprod(x) / n
  ols <- qr.solve(x, y)
  reference <- lapply(tau_values, function(tau) {
    fit <- require_success(fit_qr_smooth(x, y, rep(1, n), tau, n_total = n,
      start = respondent_start(x, y, tau, ols), tolerance = n^(-1.5)), "full-record reference")
    list(beta = fit$beta, loss = smoothed_empirical_loss(x, y, fit$beta, tau))
  })
  mar_design <- cbind(x, base$z, base$z^2)
  checkpoints <- file.path(out, "replications")
  dir.create(checkpoints, showWarnings = FALSE)
  run_rep <- function(rep) {
    path <- file.path(checkpoints, sprintf("rep_%03d.rds", rep))
    if (file.exists(path)) {
      checkpoint <- readRDS(path)
      if (isTRUE(checkpoint$n == n) && checkpoint$rep == rep) return(checkpoint)
      stop("Saved replication is incompatible with the current cohort: ", path)
    }
    set.seed(2026110101L + rep)
    uniforms <- runif(n)
    rows <- coeff <- masks <- list()
    for (j in seq_len(nrow(cells))) {
      pi <- probability[[j]]
      # The independent R stream pairs all six settings within each repetition.
      # Figure 3 uses the original fixed mask preserved with the cohort.
      r <- uniforms < pi
      xr <- x[r, , drop = FALSE]
      yr <- y[r]
      if (nrow(xr) <= ncol(xr)) stop("Too few respondents in replication ", rep)
      mar <- fit_mar(mar_design, r)
      mnar <- fit_app_selection(base, r)
      if (!isTRUE(mar$converged) || !isTRUE(mnar$converged))
        stop("Calibration did not converge in replication ", rep, ", cell ", j, call. = FALSE)
      weights <- list(CC = rep(1, sum(r)), MAR = mar$weights, MNAR = mnar$weights,
                      Oracle = 1 / pi[r])
      olsr <- qr.solve(xr, yr)
      key <- cbind(cells[j, , drop = FALSE], data.frame(rep = rep, n = n, n_observed = sum(r)))
      rownames(key) <- NULL
      masks[[j]] <- cbind(key, data.frame(missing_fraction_realized = mean(!r),
        gamma_estimate = tail(mnar$theta, 1L), gradient_norm = sqrt(sum(mnar$gradient^2)),
        overid_statistic = mnar$overid_statistic, overid_pvalue = mnar$overid_pvalue,
        whitened_condition = mnar$whitened_condition))
      for (method in names(weights)) {
        w <- weights[[method]]
        for (k in seq_along(tau_values)) {
          tau <- tau_values[[k]]
          fit <- require_success(fit_qr_smooth(xr, yr, w, tau, n_total = n,
            start = respondent_start(xr, yr, tau, olsr), tolerance = n^(-1.5)), paste(rep, j, method, tau))
          diff <- fit$beta - reference[[k]]$beta
          excess <- 40 * (smoothed_empirical_loss(x, y, fit$beta, tau) - reference[[k]]$loss)
          if (excess < -80 * n^(-1.5)) stop("Reference loss failed accuracy check in replication ", rep)
          rows[[length(rows) + 1L]] <- cbind(key, data.frame(method = method, tau = tau,
            prediction_mse_pounds2 = 40^2 * max(0, drop(crossprod(diff, gram %*% diff))),
            coefficient_mse_pounds2 = mean((40 * diff)^2), excess_risk_pounds = excess,
            as.list(weight_stats(w))))
          coeff[[length(coeff) + 1L]] <- cbind(key, data.frame(method = method, tau = tau,
            coefficient = coefficient_names, estimate_lb = pounds(fit$beta), error_lb = 40 * diff))
        }
      }
    }
    value <- list(n = n, rep = rep, metrics = do.call(rbind, rows),
                  coefficients = do.call(rbind, coeff), masks = do.call(rbind, masks))
    temp <- tempfile(tmpdir = checkpoints)
    saveRDS(value, temp)
    if (!file.rename(temp, path)) stop("Unable to write the replication file: ", path)
    message("Completed repeated-mask replication ", rep)
    value
  }
  if (workers > 1L && .Platform$OS.type == "windows")
    stop("Use --workers 1 on Windows.", call. = FALSE)
  result <- if (workers == 1L) lapply(seq_len(reps), run_rep) else
    parallel::mclapply(seq_len(reps), run_rep, mc.cores = workers, mc.set.seed = FALSE,
                       mc.preschedule = FALSE)
  failed <- vapply(result, inherits, logical(1), "try-error")
  if (any(failed)) stop("At least one repeated-mask fit failed: ", as.character(result[[which(failed)[1L]]]))
  metrics <- do.call(rbind, lapply(result, `[[`, "metrics"))
  write.csv(metrics, file.path(out, "mask_level_metrics.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(result, `[[`, "coefficients")),
            file.path(out, "mask_level_coefficients.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(result, `[[`, "masks")),
            file.path(out, "mask_level_selection.csv"), row.names = FALSE)
  summarize_recovery(metrics, bootstrap_reps, reps, out)
  invisible(metrics)
}

run_natural <- function(base, out, partitions, primary_only) {
  n <- length(base$y)
  r <- base$r
  xr <- base$x[r, , drop = FALSE]
  yr <- base$y[r]
  taus <- c(.1, .5, .9)
  write.csv(base$size, file.path(out, "cohort.csv"), row.names = FALSE)
  write.csv(base$descriptives, file.path(out, "cohort_descriptives.csv"), row.names = FALSE)
  rows <- selections <- list()
  fit_set <- function(method, resolution, delta, selection = NULL, start = NULL) {
    weights <- if (is.null(selection)) rep(1, sum(r)) else selection$weights
    fitted <- vector("list", length(taus))
    for (k in seq_along(taus)) {
      tau <- taus[[k]]
      initial <- if (is.null(start)) c(unname(quantile(yr, tau)), 0, 0, .5, rep(0, 7)) else start[[k]]
      a <- resolution / 40
      fit <- require_success(fit_qr_smooth(xr, yr, weights, tau, n_total = n, a = a,
        start = initial, partitions = partitions, tolerance = n^(-1.5)), paste(method, resolution, delta, tau))
      uncertainty <- qr_inference(base$x, base$y, r, fit$beta, tau, fit = selection,
                                  bandwidth = a, smooth = TRUE, robust = TRUE, partitions = partitions)
      est <- pounds(fit$beta)
      se <- 40 * drop(uncertainty$se)
      rows[[length(rows) + 1L]] <<- data.frame(method = method, tau = tau, a_lb = resolution,
        delta_z = delta, coefficient = coefficient_names, estimate_lb = est, se_lb = se,
        ci_low_lb = est - 1.96 * se, ci_high_lb = est + 1.96 * se,
        weighted_loss = fit$objective, converged = fit$converged)
      fitted[[k]] <- fit$beta
      message(sprintf("%s tau=%.1f a=%.2f delta=%.2f", method, tau, resolution, delta))
    }
    fitted
  }
  fit_set("CC", .5, 0)
  mar <- fit_app_selection(base, r, partitions, mode = "MAR")
  if (!isTRUE(mar$converged)) stop("MAR working calibration failed.")
  selections[[length(selections) + 1L]] <- app_diagnostics(mar, "MAR", 0)
  fit_set("MAR", .5, 0, mar)
  rm(mar)
  mnar <- fit_app_selection(base, r, partitions)
  if (!isTRUE(mnar$converged)) stop("MNAR working calibration failed.")
  selections[[length(selections) + 1L]] <- app_diagnostics(mnar, "MNAR", 0)
  fit_primary <- fit_set("MNAR", .5, 0, mnar)
  if (!primary_only) {
    for (resolution in c(.25, 1)) fit_set("MNAR", resolution, 0, mnar, fit_primary)
    warm <- mnar$theta
    rm(mnar)
    for (delta in c(-.5, -.25, .25, .5)) {
      fit_delta <- fit_app_selection(base, r, partitions, delta = delta, warm = warm)
      if (!isTRUE(fit_delta$converged)) stop("Sensitivity calibration failed at delta = ", delta)
      selections[[length(selections) + 1L]] <- app_diagnostics(fit_delta, "MNAR", delta)
      fit_set("MNAR", .5, delta, fit_delta, fit_primary)
    }
  }
  write.csv(do.call(rbind, rows), file.path(out, "coefficients.csv"), row.names = FALSE)
  write.csv(do.call(rbind, selections), file.path(out, "selection_diagnostics.csv"), row.names = FALSE)
  invisible(rows)
}

app_diagnostics <- function(fit, method, delta) {
  stat <- fit$overid_statistic
  degrees <- fit$overid_df
  data.frame(method = method, delta_z = delta, parameter = seq_along(fit$theta), theta = fit$theta,
             gradient_norm = sqrt(sum(fit$gradient^2)), at_boundary = fit$at_boundary,
             whitened_condition = fit$whitened_condition,
             overid_statistic = stat, overid_df = degrees,
             overid_pvalue = if (degrees > 0) pchisq(stat, degrees, lower.tail = FALSE) else NA_real_,
             converged = fit$converged)
}

application_main <- function(root = get_script_directory(), args = commandArgs(TRUE)) {
  opts <- parse_app_args(args, list(mode = "fixed", out = "", reps = 100L, workers = 1L,
    bootstrap = 5000L, partitions = 20L, limit = .Machine$integer.max, `primary-only` = FALSE))
  if (is.null(opts)) return(invisible(NULL))
  mode <- match.arg(opts$mode, c("fixed", "repeated", "natural"))
  shared <- file.path(root, "..", "simulation")
  source(file.path(shared, "calibration.R"))
  source(file.path(shared, "quantile.R"))
  out <- if (nzchar(opts$out)) opts$out else file.path(root, "results", switch(mode,
    fixed = "fixed_mask", repeated = "repeated_masks", natural = "natural_missingness"))
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  switch(mode,
    fixed = run_fixed_mask(read_fixed_cohort(root, opts$limit), out, opts$partitions),
    repeated = run_repeated_masks(read_fixed_cohort(root, opts$limit), out, opts$reps, opts$workers,
                                  opts$bootstrap),
    natural = run_natural(load_natural_cohort(root, opts$limit), out, opts$partitions,
                          opts[["primary-only"]]))
}

if (sys.nframe() == 0L) application_main()
