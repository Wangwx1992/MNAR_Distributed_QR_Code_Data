# From the extracted package directory:
#   Rscript simulation/run_simulation.R primary [--reps=500] [--workers=1]
#   Rscript simulation/run_simulation.R computing [--reps=100] [--workers=1]
#   Rscript simulation/run_simulation.R all --smoke
#   Rscript simulation/make_displays.R all [--smoke]
# Uses base R. For comparable compute timings, start R with one BLAS thread
# (e.g. OPENBLAS_NUM_THREADS=1).
# Timing values are specific to the current hardware, R and BLAS.
args <- commandArgs(trailingOnly = TRUE)
if (any(args %in% c("-h", "--help"))) {
  cat("Usage: Rscript simulation/run_simulation.R [all|primary|computing]\n",
      "           [--smoke] [--reps=INTEGER] [--workers=INTEGER]\n",
      "Defaults: 500 primary / 100 computing replicates, one worker.\n",
      "--smoke uses smaller designs and separate output directories.\n", sep = "")
  quit(save = "no")
}
pos <- args[!startsWith(args, "--")]
mode <- if (length(pos)) pos[[1L]] else "all"
if (length(pos) > 1L || !mode %in% c("all", "primary", "computing") ||
    any(!grepl("^(primary|computing|all|--smoke|--reps=[0-9]+|--workers=[0-9]+)$", args))) {
  stop("Invalid arguments. See --help.")
}
integer_arg <- function(name, default) {
  values <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(values)) return(default)
  if (length(values) != 1L) stop("Specify --", name, " once.")
  value <- suppressWarnings(as.integer(sub("^[^=]+=", "", values)))
  if (!is.finite(value) || value < 1L) stop("--", name, " requires a positive integer.")
  value
}
workers <- integer_arg("workers", 1L)
if (workers > 1L && .Platform$OS.type == "windows") {
  stop("For a Windows run, specify --workers=1.")
}
smoke <- "--smoke" %in% args
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
here <- dirname(normalizePath(script, mustWork = TRUE))
source(file.path(here, "data_generation.R"))
source(file.path(here, "calibration.R"))
source(file.path(here, "quantile.R"))

elapsed <- function() unname(proc.time()[["elapsed"]])
save_result <- function(value, path) {
  tmp <- tempfile(pattern = ".case-", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(value, tmp)
  if (!file.rename(tmp, path)) stop("Cannot save ", path)
}

primary_case <- function(design, replication, reference) {
  data <- generate_primary(design, replication, reference)
  n <- design$N
  r <- data$r
  x <- data$x
  y <- data$y
  t <- cbind(x[r, 1:3, drop = FALSE], primary_transform(y[r]))
  cal <- fit_mnar(t, data$h, r, partitions = 20L, bounds = 4, tolerance = 1 / n)
  mar <- fit_mar(cbind(x[, 1:3, drop = FALSE], data$z), r)
  methods <- list(MNAR = list(r = r, fit = cal, weights = cal$weights),
                  CC = list(r = r, fit = NULL, weights = rep(1, sum(r))),
                  MAR = list(r = r, fit = mar, weights = mar$weights),
                  Oracle = list(r = r, fit = NULL, weights = 1 / data$pi[r]),
                  Full = list(r = rep(TRUE, n), fit = NULL, weights = rep(1, n)))
  fits <- list()
  for (tau in c(.1, .5, .9)) {
    truth <- primary_truth(tau, design$error)
    for (method in names(methods)) {
      spec <- methods[[method]]
      idx <- spec$r
      fit <- fit_qr_central(x[idx, , drop = FALSE], y[idx], spec$weights, tau, n_total = n)
      if (is.null(spec$fit)) {
        variance <- qr_inference(x, y, idx, fit$beta, tau, weights = spec$weights)
      } else {
        variance <- qr_inference(x, y, idx, fit$beta, tau, fit = spec$fit)
      }
      error <- fit$beta - truth
      fits[[length(fits) + 1L]] <- data.frame(
        scenario = design$scenario, N = n, error_distribution = design$error,
        missing_fraction = design$missing, realized_missing_fraction = mean(!r),
        rep = replication, tau = tau, method = method, coefficient = seq_along(truth) - 1L,
        beta = fit$beta, truth = truth, error = error, se = variance$se,
        covered95 = abs(error) <= qnorm(.975) * variance$se,
        interval_length = 2 * qnorm(.975) * variance$se, fixed_se = variance$fixed_se,
        objective = fit$objective, fit_converged = isTRUE(fit$converged),
        calibration_converged = if (is.null(spec$fit)) TRUE else isTRUE(spec$fit$converged),
        calibration_gradient = if (is.null(spec$fit)) 0 else spec$fit$gradient
      )
    }
  }
  do.call(rbind, fits)
}

computing_case <- function(design, replication, reference) {
  data <- generate_computing(design, replication, reference)
  n <- design$N
  start <- elapsed()
  cal <- fit_mnar(data$t, data$h, data$r, partitions = 20L,
                  bounds = c(4, rep(3, 6L)), tolerance = 1 / n)
  t_w <- elapsed() - start
  x <- data$x
  y <- data$y
  w <- cal$weights
  partitions <- if (smoke) c(5L, 10L) else c(5L, 10L, 20L, 50L, 100L)
  fits <- list()
  for (tau in c(.3, .5, .7)) {
    truth <- computing_truth(data$beta, tau, design$hetero)
    start <- elapsed()
    ip <- fit_qr_central(x, y, w, tau, n_total = n, tolerance = 1e-8)
    ip_seconds <- elapsed() - start
    record <- function(fit, method, M, serial, displayed) {
      data.frame(scenario = design$scenario, N = n, p = design$p,
                 heteroscedastic = design$hetero, rep = replication, tau = tau,
                 method = method, M = M, AE = sum(abs(fit$beta - truth)),
                 T_W = t_w, T_QR = displayed, T_total = t_w + displayed,
                 T_QR_serial = serial, T_total_serial = t_w + serial,
                 fit_converged = isTRUE(fit$converged),
                 calibration_converged = isTRUE(cal$converged),
                 calibration_gradient = cal$gradient,
                 coefficient_norm = sqrt(sum(fit$beta ^ 2)),
                 objective_excess = fit$objective - ip$objective,
                 max_coefficient_difference = max(abs(fit$beta - ip$beta)),
                 l1_coefficient_difference = sum(abs(fit$beta - ip$beta)))
    }
    fits[[length(fits) + 1L]] <- record(ip, "IP", 0L, ip_seconds, ip_seconds)
    order <- if (replication %% 2L) c("ADMM", "Renewable") else c("Renewable", "ADMM")
    for (M in partitions) {
      for (method in order) {
        start <- elapsed()
        if (method == "ADMM") {
          fit <- fit_qr_admm(x, y, w, tau, partitions = M, n_total = n,
                             tolerance = .05, stopping = "coefficient", maxit = 5000L)
          displayed <- fit$parallel_seconds
        } else {
          fit <- fit_qr_renewable(x, y, w, tau, partitions = M, n_total = n,
                                  tolerance = 1e-4, maxit = 5000L)
          displayed <- NA_real_
        }
        serial <- elapsed() - start
        if (method != "ADMM") displayed <- serial
        fits[[length(fits) + 1L]] <- record(fit, method, M, serial, displayed)
      }
    }
  }
  do.call(rbind, fits)
}

run_study <- function(study, designs, repetitions, get_reference, fit_case) {
  stage <- if (smoke) "smoke" else "formal"
  out <- file.path(here, "results", stage, study)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  references <- lapply(seq_len(nrow(designs)), function(i) {
    get_reference(designs[i, , drop = FALSE])
  })
  tasks <- expand.grid(design = seq_len(nrow(designs)), rep = seq_len(repetitions))
  run_task <- function(i) {
    row <- tasks[i, ]
    design <- designs[row$design, , drop = FALSE]
    path <- file.path(out, sprintf("%s_rep%04d.rds", design$scenario, row$rep))
    if (file.exists(path)) return(TRUE)
    tryCatch({
      value <- fit_case(design, row$rep, references[[row$design]])
      save_result(value, path)
      message(study, ": ", basename(path))
      TRUE
    }, error = function(e) {
      message(study, ": ", basename(path), ": ", conditionMessage(e))
      FALSE
    })
  }
  seq <- seq_len(nrow(tasks))
  completed <- if (workers == 1L) lapply(seq, run_task) else
    parallel::mclapply(seq, run_task, mc.cores = workers, mc.preschedule = FALSE,
                       mc.set.seed = FALSE)
  failures <- sum(!vapply(completed, isTRUE, logical(1L)))
  if (failures) stop(failures, " replication(s) failed in ", study,
                     ". Inspect the error messages and rerun to fill missing results.")
}

if (mode %in% c("all", "primary")) {
  reps <- integer_arg("reps", if (smoke) 2L else 500L)
  if (reps > 500L) stop("The primary design has at most 500 replicates.")
  run_study("primary", primary_designs(smoke), reps,
             function(d) primary_reference(d, smoke), primary_case)
}
if (mode %in% c("all", "computing")) {
  reps <- integer_arg("reps", if (smoke) 2L else 100L)
  if (reps > 100L) stop("The computing design has at most 100 replicates.")
  run_study("computing", computing_designs(smoke), reps,
             function(d) computing_reference(d$p, d$hetero, smoke), computing_case)
}
