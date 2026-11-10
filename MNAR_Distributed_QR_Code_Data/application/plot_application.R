# Draw Figure 3 from the fixed-mask analysis.
# Rscript application/run_application.R --mode fixed
# Rscript application/plot_application.R
# Supply --input PATH --output PATH to plot a separately stored run.

plot_application <- function(input, output) {
  d <- read.csv(input, stringsAsFactors = FALSE)
  keys <- c("method", "M", "tau")
  if (nrow(d) != 63L || anyDuplicated(d[, keys]))
    stop("Figure 3 needs exactly 63 distinct method, partition and quantile rows.")
  curves <- data.frame(
    method = c("IP", rep("ADMM", 3), rep("Renewable", 3)),
    M = c(0L, 20L, 50L, 100L, 20L, 50L, 100L),
    label = c("IP", paste0("IPW-ADMM (M = ", c(20, 50, 100), ")"),
              paste0("IPW-renewable (M = ", c(20, 50, 100), ")")),
    color = c("black", "#de3434", "#19a33a", "#2246f0", "#09b7c8", "#b72ad9", "#e6a021"),
    shape = c(1, 2, 3, 4, 5, 6, 0))
  by_curve <- lapply(seq_len(nrow(curves)), function(j) {
    sub <- d[d$method == curves$method[[j]] & d$M == curves$M[[j]], , drop = FALSE]
    sub <- sub[order(sub$tau), , drop = FALSE]
    if (!isTRUE(all.equal(sub$tau, seq(.1, .9, by = .1)))) stop("Incomplete quantile grid.")
    sub
  })
  if (!all(is.finite(d$height_coefficient_lb)) ||
      !all(is.finite(d$total_figure_seconds) & d$total_figure_seconds >= 0))
    stop("Non-finite coefficient or invalid timing in the fixed-mask results.")
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  if (capabilities("cairo")) {
    cairo_pdf(output, width = 6.3, height = 7.6, family = "sans")
  } else {
    pdf(output, width = 6.3, height = 7.6, family = "Helvetica")
  }
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(2, 1), mar = c(5.5, 4.6, 2.1, 1.0), mgp = c(2.5, .7, 0),
      tcl = -.22, family = "sans")
  ylim_coef <- range(d$height_coefficient_lb)
  ylim_coef[[2L]] <- ylim_coef[[2L]] + .18 * max(diff(ylim_coef), 1)
  panel <- function(key, y_label, limits, caption) {
    plot(NA, xlim = c(.07, .93), ylim = limits, xaxt = "n",
         xlab = "Tau", ylab = y_label, main = "Maternal height (per five inches)")
    axis(1, at = c(.2, .4, .6, .8))
    for (j in seq_along(by_curve)) {
      curve <- by_curve[[j]]
      lines(curve$tau, curve[[key]], type = "b", lty = 2, col = curves$color[[j]],
            pch = curves$shape[[j]], lwd = .85, cex = .65)
    }
    legend("topleft", legend = curves$label, col = curves$color, pch = curves$shape,
           lty = 2, lwd = .85, pt.cex = .65, cex = .72, bg = "white", seg.len = 2)
    mtext(caption, side = 1, line = 4.1)
  }
  panel("height_coefficient_lb", "Estimated coefficient (lb)", ylim_coef,
        "(a) Estimation performance")
  panel("total_figure_seconds", "Time (s)", c(0, max(1, 2 * d$total_figure_seconds)),
        "(b) Time performance")
  invisible(output)
}

plot_application_main <- function(args = commandArgs(TRUE)) {
  if (any(args == "--help")) {
    cat("Usage: Rscript application/plot_application.R [--input CSV] [--output PDF]\n")
    return(invisible(NULL))
  }
  invocation <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  root <- dirname(normalizePath(sub("^--file=", "", invocation[[1L]])))
  path <- file.path(root, "results", "fixed_mask")
  settings <- list(input = file.path(path, "summary.csv"), output = file.path(path, "figure3.pdf"))
  if (length(args) %% 2L) stop("Every option requires a value.")
  for (i in seq_along(args)[seq_along(args) %% 2L == 1L]) {
    key <- sub("^--", "", args[[i]])
    if (!startsWith(args[[i]], "--") || !key %in% names(settings)) stop("Unknown option: ", args[[i]])
    settings[[key]] <- args[[i + 1L]]
  }
  do.call(plot_application, settings)
}

if (sys.nframe() == 0L) plot_application_main()
