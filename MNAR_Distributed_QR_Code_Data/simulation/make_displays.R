# Generate the simulation tables and boxplots from run_simulation.R outputs.
# Rscript simulation/make_displays.R [all|primary|computing] [--smoke] [--reps=INTEGER]
args <- commandArgs(trailingOnly = TRUE)
if (any(args %in% c("-h", "--help"))) {
  cat("Usage: Rscript simulation/make_displays.R [all|primary|computing]\n",
      "           [--smoke] [--reps=INTEGER]\n",
      "Use the same study, --smoke and --reps options as run_simulation.R.\n", sep = "")
  quit(save = "no")
}
pos <- args[!startsWith(args, "--")]
mode <- if (length(pos)) pos[[1L]] else "all"
if (length(pos) > 1L || !mode %in% c("all", "primary", "computing") ||
    any(!grepl("^(primary|computing|all|--smoke|--reps=[0-9]+)$", args))) {
  stop("Invalid arguments. See --help.")
}
smoke <- "--smoke" %in% args
n_arg <- grep("^--reps=", args, value = TRUE)
if (length(n_arg) > 1L) stop("Specify --reps once.")
num_reps <- function(default, maximum) {
  n <- if (length(n_arg)) suppressWarnings(as.integer(sub("^--reps=", "", n_arg))) else default
  if (!is.finite(n) || n < 1L || n > maximum) stop("Expected 1 <= repetitions <= ", maximum)
  n
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
here <- dirname(normalizePath(script, mustWork = TRUE))
source(file.path(here, "data_generation.R"))
stage <- if (smoke) "smoke" else "formal"
out <- file.path(here, "results", stage, "displays")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

load_study <- function(study, designs, reps) {
  files <- unlist(lapply(designs$scenario, function(scenario) {
    file.path(here, "results", stage, study, sprintf("%s_rep%04d.rds", scenario, seq_len(reps)))
  }), use.names = FALSE)
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("Missing ", length(missing), " replication files. First: ", missing[1L])
  do.call(rbind, lapply(files, readRDS))
}

primary_methods <- c("MNAR", "CC", "MAR", "Oracle", "Full")
primary_errors <- c("normal", "t3", "lognormal")
primary_taus <- c(.1, .5, .9)

make_primary_table <- function(point, n, reps) {
  data <- point[point$N == n & point$missing_fraction == .4, ]
  caption <- sprintf(r"(\caption{Coefficient root mean squared error and intercept coverage (percent, in parentheses) for $N=%d$ and $40\%%$ missing responses. Each entry uses %d replications.})", n, reps)
  lines <- c(r"(\begin{table}[t])", r"(\centering\small)", caption,
             r"(\label{tab:sim-primary})", r"(\begin{tabular}{llccccc}\toprule)",
             r"(Error & $\tau$ & MNAR & Complete case & MAR IPW & Oracle IPW & Full data\\\midrule)")
  labels <- c("Normal", r"($t_3/\sqrt{3}$)", "Lognormal")
  for (i in seq_along(primary_errors)) {
    for (j in seq_along(primary_taus)) {
      tau <- primary_taus[j]
      cells <- c(if (j == 1L) labels[i] else "", formatC(tau, digits = 1L, format = "f"))
      for (method in primary_methods) {
        rows <- data[data$error_distribution == primary_errors[i] & data$tau == tau &
                     data$method == method, ]
        stopifnot(nrow(rows) == reps, length(unique(rows$rep)) == reps,
                  all(is.finite(rows$coefficient_mse)), !anyNA(rows$intercept_covered))
        cells <- c(cells, sprintf("%.4f (%.1f)", sqrt(mean(rows$coefficient_mse)),
                                  100 * mean(rows$intercept_covered)))
      }
      lines <- c(lines, paste(paste(cells, collapse = " & "), r"(\\)"))
    }
  }
  writeLines(c(lines, r"(\bottomrule\end{tabular})", r"(\end{table})"),
             file.path(out, "simulation_table_primary.tex"))
}

draw_box <- function(values, at, color, width = .105) {
  quartile <- as.numeric(quantile(values, c(.25, .5, .75), type = 7))
  iqr <- quartile[3L] - quartile[1L]
  inside <- values >= quartile[1L] - 1.5 * iqr & values <= quartile[3L] + 1.5 * iqr
  low <- min(values[inside])
  high <- max(values[inside])
  segments(at, low, at, quartile[1L], col = color, lwd = .65)
  segments(at, quartile[3L], at, high, col = color, lwd = .65)
  segments(at - width / 4, low, at + width / 4, low, col = color, lwd = .65)
  segments(at - width / 4, high, at + width / 4, high, col = color, lwd = .65)
  rect(at - width / 2, quartile[1L], at + width / 2, quartile[3L], col = color, border = color)
  segments(at - width / 2, quartile[2L], at + width / 2, quartile[2L], col = "white", lwd = .8)
  if (any(!inside)) points(rep(at, sum(!inside)), values[!inside], col = color, pch = 1L, cex = .35)
}

make_primary_figure <- function(point, n, reps) {
  data <- point[point$N == n & point$missing_fraction == .4, ]
  stopifnot(nrow(data) == 3L * 3L * 5L * reps)
  data$coefficient_rms_error <- sqrt(data$coefficient_mse)
  write.csv(data, file.path(out, "boxplot_data.csv"), row.names = FALSE)
  colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#666666")
  legend_labels <- c("Calibrated MNAR", "Complete case", "MAR IPW", "Known propensity", "Full data")
  titles <- expression("(a) Normal", "(b) " * t[3] / sqrt(3), "(c) Right skewed")
  ymax <- max(1e-8, max(data$coefficient_rms_error) * 1.06)
  tick <- pretty(c(0, ymax), n = 4)
  tick <- tick[tick >= 0 & tick <= ymax]
  output <- file.path(out, "fig_simulation_boxplots.pdf")
  if (capabilities("cairo")) {
    cairo_pdf(output, width = 7.15, height = 3.55, family = "sans")
  } else {
    pdf(output, width = 7.15, height = 3.55, family = "Helvetica", useDingbats = FALSE)
  }
  on.exit(dev.off(), add = TRUE)
  layout(rbind(seq_len(3L), rep(4L, 3L)), heights = c(5, 1))
  for (i in seq_along(primary_errors)) {
    par(mar = c(3.5, if (i == 1L) 4.7 else 1.3, 2.2, .8),
        mgp = c(1.9, .45, 0), tcl = -.22, cex = .83, xaxs = "i", yaxs = "i")
    plot.new()
    plot.window(xlim = c(.55, 3.45), ylim = c(0, ymax))
    abline(h = tick, col = "#e2e2e2", lwd = .5)
    for (j in seq_along(primary_taus)) {
      for (k in seq_along(primary_methods)) {
        group <- data[data$error_distribution == primary_errors[i] &
                      data$tau == primary_taus[j] & data$method == primary_methods[k], ]
        stopifnot(nrow(group) == reps, length(unique(group$rep)) == reps)
        draw_box(group$coefficient_rms_error, j + (k - 3) * .14, colors[k])
      }
    }
    axis(1L, at = 1:3, labels = c("0.1", "0.5", "0.9"))
    axis(2L, at = tick, labels = if (i == 1L) formatC(tick, digits = 2L, format = "f") else FALSE,
         las = 1L)
    title(main = titles[i], xlab = expression("Quantile level " * tau), font.main = 1L)
    if (i == 1L) mtext("Coefficient RMS error per repetition", side = 2L, line = 3.3, cex = .83)
  }
  par(mar = rep(0, 4), cex = .77)
  plot.new()
  legend("center", legend = legend_labels, fill = colors, border = colors,
         bty = "n", ncol = 3L, x.intersp = .6, y.intersp = 1.15)
}

primary_displays <- function() {
  reps <- num_reps(if (smoke) 2L else 500L, 500L)
  designs <- primary_designs(smoke)
  d <- load_study("primary", designs, reps)
  stopifnot(nrow(d) == nrow(designs) * reps * 3L * 5L * 5L,
            all(is.finite(d$error)), !anyNA(d$covered95))
  key <- c("scenario", "N", "error_distribution", "missing_fraction", "rep", "tau", "method")
  point <- d[d$coefficient == 0L, c(key, "covered95"), drop = FALSE]
  names(point)[names(point) == "covered95"] <- "intercept_covered"
  mse <- aggregate(list(coefficient_mse = d$error ^ 2), by = d[key], FUN = mean)
  point <- merge(point, mse, by = key, sort = FALSE)
  sizes <- sort(unique(designs$N))
  make_primary_table(point, sizes[1L], reps)
  make_primary_figure(point, sizes[2L], reps)
  coef_key <- c("scenario", "N", "error_distribution", "missing_fraction", "tau", "method", "coefficient")
  summary <- aggregate(list(coefficient_mse = d$error ^ 2, coverage95 = d$covered95,
                            interval_length = d$interval_length,
                            fit_failures = !d$fit_converged,
                            calibration_failures = !d$calibration_converged),
                       by = d[coef_key], FUN = mean)
  summary$coefficient_rmse <- sqrt(summary$coefficient_mse)
  summary$fit_failures <- round(summary$fit_failures * reps)
  summary$calibration_failures <- round(summary$calibration_failures * reps)
  write.csv(summary, file.path(out, "primary_coefficient_summary.csv"), row.names = FALSE)
  if (any(!d$fit_converged | !d$calibration_converged)) {
    warning("Some primary fits did not satisfy convergence checks. All finite results are shown; see primary_coefficient_summary.csv.")
  }
}

mean_sd <- function(values) {
  if (length(values) > 1L) c(mean(values), sd(values)) else c(mean(values), NA_real_)
}

computing_displays <- function() {
  reps <- num_reps(if (smoke) 2L else 100L, 100L)
  designs <- computing_designs(smoke)
  m_set <- if (smoke) c(5L, 10L) else c(5L, 10L, 20L, 50L, 100L)
  d <- load_study("computing", designs, reps)
  stopifnot(nrow(d) == nrow(designs) * reps * 3L * (1L + 2L * length(m_set)))
  keys <- c("N", "p", "heteroscedastic", "tau", "method", "M")
  metrics <- c("AE", "T_W", "T_QR", "T_total", "T_QR_serial", "T_total_serial",
               "objective_excess", "max_coefficient_difference", "l1_coefficient_difference")
  stopifnot(all(vapply(d[metrics], function(x) all(is.finite(x)), logical(1L))))
  groups <- split(d, interaction(d[keys], drop = TRUE))
  s <- do.call(rbind, lapply(groups, function(group) {
    stopifnot(nrow(group) == reps, length(unique(group$rep)) == reps)
    record <- group[1L, keys]
    record$replications <- nrow(group)
    record$stopping_failures <- sum(!group$fit_converged)
    record$calibration_failures <- sum(!group$calibration_converged)
    for (metric in metrics) {
      stat <- mean_sd(group[[metric]])
      record[[paste0(metric, "_mean")]] <- stat[1L]
      record[[paste0(metric, "_sd")]] <- stat[2L]
    }
    record
  }))
  write.csv(s, file.path(out, "computing_summary.csv"), row.names = FALSE)
  names_method <- c(ADMM = "MNAR-IPW-ADMM", Renewable = "MNAR-IPW-renewable", IP = "Central IP")
  for (hetero in c(FALSE, TRUE)) {
    err <- if (hetero) "heteroscedastic" else "homoscedastic"
    panels <- designs[designs$hetero == hetero, c("N", "p"), drop = FALSE]
    for (tau in c(.3, .5, .7)) {
      code <- sprintf("%02d", as.integer(round(tau * 10)))
      lines <- c(r"(\begin{table}[p])", r"(\centering\footnotesize)",
                 sprintf(r"(\caption{Estimation error and computation time under %s errors at $\tau=%.1f$.})",
                         err, tau),
                 sprintf(r"(\label{tab:compute_%s_%s})", err, code),
                 r"(\begin{tabular}{llrrrr})", r"(\toprule)",
                 r"(Method & $M$ & AE & $T_{\mathrm{total}}$ (s) & $T_W$ (s) & $T_{\mathrm{QR}}$ (s) \\)",
                 r"(\midrule)")
      for (i in seq_len(nrow(panels))) {
        n <- panels$N[i]
        p <- panels$p[i]
        lines <- c(lines, sprintf(r"(\multicolumn{6}{c}{$(N,p)=(%d,%d)$} \\)", n, p))
        for (method in names(names_method)) {
          partitions <- if (method == "IP") 0L else m_set
          for (j in seq_along(partitions)) {
            m <- partitions[j]
            row <- s[s$N == n & s$p == p & s$heteroscedastic == hetero &
                       s$tau == tau & s$method == method & s$M == m, , drop = FALSE]
            stopifnot(nrow(row) == 1L)
            cells <- c(if (j == 1L) names_method[[method]] else "",
                       if (method == "IP") "--" else as.character(m))
            for (metric in c("AE", "T_total", "T_W", "T_QR")) {
              cells <- c(cells, sprintf("%.3f (%.3f)", row[[paste0(metric, "_mean")]],
                                        row[[paste0(metric, "_sd")]]))
            }
            lines <- c(lines, paste(paste(cells, collapse = " & "), r"(\\)"))
          }
        }
        if (i < nrow(panels)) lines <- c(lines, r"(\midrule)")
      }
      lines <- c(lines, r"(\bottomrule)", r"(\end{tabular})", r"(\par\smallskip)",
                 r"(\begin{minipage}{0.97\textwidth}\scriptsize)",
                 sprintf("Entries are means (standard deviations) across %d replications.", reps),
                 r"(AE is the $\ell_1$ coefficient error, including the intercept. All fits use the same estimated GMM weights. $T_W$ is the measured common calibration time. ADMM $T_{\mathrm{QR}}$ models an ideal zero-network parallel critical path. Central and renewable $T_{\mathrm{QR}}$ are measured serial wall times. $T_{\mathrm{total}}=T_W+T_{\mathrm{QR}}$.)",
                 r"(\end{minipage})", r"(\end{table})")
      writeLines(lines, file.path(out, sprintf("table_%s_tau%s.tex", err, code)))
    }
  }
}

if (mode %in% c("all", "primary")) primary_displays()
if (mode %in% c("all", "computing")) computing_displays()
cat("Simulation displays saved in ", out, "\n", sep = "")
