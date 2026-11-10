# Two-step inverse-probability calibration and quantile sandwich estimates.
# Site reductions use all records, while response-dependent arrays are indexed
# only on respondents. All criteria use the original total sample size.

cal_norm <- function(x) sqrt(sum(x * x))

cal_blocks <- function(n, partitions) {
  if (length(partitions) != 1L || !is.finite(partitions) ||
      partitions != as.integer(partitions) || partitions < 1L || partitions > n)
    stop("partitions must be an integer between 1 and the number of records.")
  cut <- floor(seq(0, n, length.out = partitions + 1L))
  lapply(seq_len(partitions), function(k) seq.int(cut[k] + 1, cut[k + 1]))
}

cal_inverse <- function(a, label) {
  chol_factor <- tryCatch(chol((a + t(a)) / 2), error = function(e) NULL)
  if (is.null(chol_factor)) stop(label, " is not positive definite.")
  chol2inv(chol_factor)
}

fit_mnar <- function(t, h, r, partitions = 1L, bounds = 4,
                     tolerance = 1 / length(r), starts = NULL, offset = 0,
                     influence = FALSE,
                     initialization = c("grid", "multistart"),
                     outcome_dependent = TRUE) {
  initialization <- match.arg(initialization)
  t <- as.matrix(t)
  h <- as.matrix(h)
  r <- as.logical(r)
  n <- length(r)
  q <- ncol(t)
  ell <- ncol(h)
  if (n < 2L || anyNA(r) || !any(r) || all(r) || nrow(h) != n ||
      ell < q || q == 0L || !all(is.finite(h)) ||
      !is.finite(tolerance) || tolerance <= 0)
    stop("Invalid response indicator, moment dictionary, or gradient tolerance.")
  if (nrow(t) == n) t <- t[r, , drop = FALSE]
  if (nrow(t) != sum(r) || !all(is.finite(t)))
    stop("t must have one finite row per respondent.")
  if (length(offset) == n) offset <- offset[r]
  if (length(offset) == 1L) offset <- rep(offset, nrow(t))
  if (length(offset) != nrow(t) || !all(is.finite(offset)))
    stop("Offset must be scalar or have one finite value per respondent.")

  if (is.matrix(bounds)) {
    if (!identical(dim(bounds), c(q, 2L))) stop("Expected q by 2 bounds.")
    lower <- bounds[, 1L]
    upper <- bounds[, 2L]
  } else {
    if (!length(bounds) %in% c(1L, q)) stop("Supply one bound or q bounds.")
    upper <- rep_len(bounds, q)
    lower <- -upper
  }
  if (!all(is.finite(c(lower, upper))) || any(lower >= upper))
    stop("Bounds must define a finite nonempty box.")
  blocks <- cal_blocks(n, partitions)
  pos <- integer(n)
  pos[r] <- seq_len(sum(r))
  site_data <- lapply(blocks, function(ids) {
    respondent <- r[ids]
    ri <- ids[respondent]
    ti <- pos[ri]
    h0 <- h[ids[!respondent], , drop = FALSE]
    list(t = t[ti, , drop = FALSE], h = h[ri, , drop = FALSE],
         off = offset[ti], missing_sum = colSums(h0), missing_cross = crossprod(h0))
  })
  missing_sum <- Reduce(`+`, lapply(site_data, `[[`, "missing_sum"))
  missing_cross <- Reduce(`+`, lapply(site_data, `[[`, "missing_cross"))

  moments <- function(theta, derivative = TRUE, second = FALSE, weight_times_m = NULL) {
    m <- -missing_sum
    if (derivative) jac <- matrix(0, ell, q)
    if (second) omega <- missing_cross
    if (!is.null(weight_times_m)) extra <- matrix(0, q, q)
    for (s in site_data) {
      odds <- exp(-drop(s$t %*% theta) - s$off)
      m <- m + drop(crossprod(s$h, odds))
      if (derivative) jac <- jac - crossprod(s$h, s$t * odds)
      if (second) omega <- omega + crossprod(s$h * odds)
      if (!is.null(weight_times_m))
        extra <- extra + crossprod(s$t, s$t * (odds * drop(s$h %*% weight_times_m)))
    }
    out <- list(m = m / n)
    if (derivative) out$G <- jac / n
    if (second) out$omega <- omega / n
    if (!is.null(weight_times_m)) out$extra <- extra / n
    out
  }

  # Evaluate the same covering grid for both criteria without repeating
  # passes over the data. Batch the linear predictors within each site.
  points <- max(2L, ceiling(log(n + exp(1))^(1 / q)))
  if (points^q > 100000) stop("The covering grid exceeds 100000 candidates.")
  axes <- Map(function(a, b) seq(a, b, length.out = points), lower, upper)
  grid <- as.matrix(do.call(expand.grid, c(axes, KEEP.OUT.ATTRS = FALSE)))
  grid_moments <- matrix(-missing_sum, ell, nrow(grid))
  for (s in site_data) {
    for (i in seq.int(1L, nrow(grid), by = 8L)) {
      idx <- i:min(i + 7L, nrow(grid))
      lp <- s$t %*% base::t(grid[idx, , drop = FALSE])
      odds <- exp(-lp - s$off)
      grid_moments[, idx] <- grid_moments[, idx, drop = FALSE] + crossprod(s$h, odds)
    }
  }
  grid_moments <- grid_moments / n

  mid <- (lower + upper) / 2
  if (initialization == "multistart") mid[1L] <- qlogis(mean(r))
  if (outcome_dependent) {
    shift <- numeric(q)
    shift[q] <- if (initialization == "multistart") 1 else min(.5, (upper[q] - lower[q]) / 4)
    candidates <- list(mid, mid - shift, mid + shift)
  } else {
    candidates <- list(mid)
  }
  if (!is.null(starts)) {
    extra_starts <- if (is.matrix(starts)) lapply(seq_len(nrow(starts)), function(i) starts[i, ]) else
      if (is.list(starts)) starts else list(starts)
    candidates <- c(candidates, extra_starts)
  }
  if (!all(vapply(candidates, function(a) length(a) == q && all(is.finite(a)), logical(1))))
    stop("Each starting value must have q finite coordinates.")
  candidates <- lapply(candidates, function(a) pmin(upper, pmax(lower, a)))

  stage <- function(W, initial) {
    cached_theta <- NULL
    cached <- NULL
    evaluate <- function(theta) {
      if (!identical(theta, cached_theta)) {
        cached_theta <<- theta
        z <- moments(theta)
        wm <- drop(W %*% z$m)
        cached <<- list(value = sum(z$m * wm) / 2, grad = drop(crossprod(z$G, wm)), G = z$G,
                        wm = wm)
      }
      cached
    }
    objective <- function(theta) evaluate(theta)$value
    gradient <- function(theta) evaluate(theta)$grad
    value0 <- objective(initial)
    estimate <- tryCatch(
      optim(initial, objective, gradient, method = "L-BFGS-B", lower = lower, upper = upper,
            control = list(maxit = 1000L, factr = 1, pgtol = tolerance / 10)),
      error = function(e) NULL
    )
    theta <- initial
    if (!is.null(estimate) && is.finite(estimate$value) && estimate$value <= value0)
      theta <- estimate$par

    # Polish a local optimum to the attained-gradient condition. The exact
    # Hessian matters for the working model, where the residual moments need
    # not vanish. If necessary use the positive Gauss–Newton part.
    for (iteration in seq_len(12L)) {
      z <- evaluate(theta)
      if (cal_norm(z$grad) <= tolerance / 4) break
      gauss_newton <- crossprod(z$G, W %*% z$G)
      correction <- moments(theta, derivative = FALSE, weight_times_m = z$wm)$extra
      step <- tryCatch(drop(solve(gauss_newton + correction, z$grad)),
                       error = function(e) numeric())
      if (!length(step) || !all(is.finite(step)) || sum(step * z$grad) <= 0)
        step <- tryCatch(drop(solve(gauss_newton, z$grad)), error = function(e) numeric())
      if (!length(step) || !all(is.finite(step)) || sum(step * z$grad) <= 0) break
      moved <- FALSE
      for (a in 2^-(0:25)) {
        candidate <- pmin(upper, pmax(lower, theta - a * step))
        if (identical(candidate, theta)) next
        new <- evaluate(candidate)
        if (is.finite(new$value) && new$value <= z$value &&
            (new$value < z$value || cal_norm(new$grad) < cal_norm(z$grad))) {
          theta <- candidate
          moved <- TRUE
          break
        }
      }
      if (!moved) break
    }
    list(theta = theta, value = objective(theta), gradient = cal_norm(gradient(theta)))
  }
  fit_stage <- function(W, initials) {
    fits <- lapply(initials, function(theta) stage(W, theta))
    fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  }
  best_grid1 <- grid[which.min(colSums(grid_moments^2)), ]
  candidate1 <- c(candidates, list(best_grid1))
  if (initialization == "grid") {
    values <- vapply(candidate1, function(a) sum(moments(a, derivative = FALSE)$m^2), numeric(1))
    candidate1 <- candidate1[which.min(values)]
  }
  first <- fit_stage(diag(ell), candidate1)
  z1 <- moments(first$theta, second = TRUE)
  W <- cal_inverse(z1$omega, "The uncentered calibration second moment")

  candidate2 <- list(first$theta)
  if (initialization == "multistart") {
    best_grid2 <- grid[which.min(colSums(grid_moments * (W %*% grid_moments))), ]
    candidate2 <- c(candidate2, candidates, list(best_grid2))
  }
  second <- fit_stage(W, candidate2)
  z2 <- moments(second$theta, second = TRUE)
  left <- crossprod(z2$G, W)
  K <- solve(left %*% z2$G, left)
  at_boundary <- any(pmin(second$theta - lower, upper - second$theta) <=
                     1e-7 * pmax(1, upper - lower))
  J <- n * sum(z2$m * drop(W %*% z2$m))
  df <- ell - q
  G_weighted <- forwardsolve(base::t(chol(z2$omega)), z2$G)
  singular <- svd(G_weighted, nu = 0L, nv = 0L)$d
  fit <- list(type = "MNAR", theta = second$theta, theta1 = first$theta,
              weights = 1 + exp(-drop(t %*% second$theta) - offset),
              gradient = second$gradient, converged = second$gradient <= tolerance && !at_boundary,
              tolerance = tolerance, at_boundary = at_boundary,
              overid_statistic = J, overid_df = df,
              overid_pvalue = if (df > 0L) pchisq(J, df, lower.tail = FALSE) else NA_real_,
              whitened_condition = max(singular) / min(singular),
              objective = second$value, moment = z2$m, W = W, K = K,
              stage1 = z1, stage2 = z2,
              data = list(t = t, h = h, r = r, offset = offset, pos = pos,
                          n = n, partitions = partitions))
  if (influence) fit$influence <- fit_mnar_influence(fit, working = FALSE)
  fit
}

fit_mar <- function(x, r, influence = FALSE) {
  x <- as.matrix(x)
  r <- as.logical(r)
  n <- length(r)
  if (nrow(x) != n || anyNA(r) || !any(r) || all(r) || !all(is.finite(x)))
    stop("Invalid MAR design or response indicator.")
  fitted <- glm.fit(x, as.numeric(r), family = binomial(), control = glm.control(epsilon = 1e-12))
  p <- plogis(drop(x %*% fitted$coefficients))
  if (!all(is.finite(p)) || any(p[r] <= 0)) stop("MAR fit has nonfinite weights.")
  A_inv <- cal_inverse(crossprod(x, x * (p * (1 - p))) / n, "The MAR Hessian")
  gradient <- cal_norm(crossprod(x, r - p) / n)
  fit <- list(type = "MAR", theta = fitted$coefficients, weights = 1 / p[r],
              converged = fitted$converged && gradient <= 1e-8, gradient = gradient,
              data = list(x = x, r = r, p = p, n = n, A_inv = A_inv))
  if (influence) fit$influence <- (x * (r - p)) %*% base::t(A_inv)
  fit
}

make_mnar_influence <- function(fit, working) {
  data <- fit$data
  w2 <- fit$weights - 1
  if (!working) return(function(rows) {
    rr <- data$r[rows]
    g <- -data$h[rows, , drop = FALSE]
    g[rr, ] <- -g[rr, , drop = FALSE] * w2[data$pos[rows[rr]]]
    -g %*% t(fit$K)
  })

  # Appendix B: contraction of the derivative of Omega1 against W m2.
  w1 <- exp(-drop(data$t %*% fit$theta1) - data$offset)
  m1 <- fit$stage1$m
  m2 <- fit$stage2$m
  G1 <- fit$stage1$G
  G2 <- fit$stage2$G
  omega1 <- fit$stage1$omega
  left <- crossprod(G2, fit$W)
  c2 <- drop(fit$W %*% m2)
  A1 <- crossprod(G1)
  A2 <- left %*% G2
  contract <- matrix(0, nrow(G1), ncol(G1))
  for (rows in cal_blocks(data$n, data$partitions)) {
    ids <- rows[data$r[rows]]
    p <- data$pos[ids]
    tt <- data$t[p, , drop = FALSE]
    hh <- data$h[ids, , drop = FALSE]
    hc <- drop(hh %*% c2)
    A1 <- A1 + crossprod(tt, tt * (w1[p] * drop(hh %*% m1))) / data$n
    A2 <- A2 + crossprod(tt, tt * (w2[p] * hc)) / data$n
    contract <- contract - 2 * crossprod(hh, tt * (w1[p]^2 * hc)) / data$n
  }
  A1inv <- cal_inverse(A1, "The first-stage working calibration Hessian")
  A2inv <- cal_inverse(A2, "The second-stage working calibration Hessian")
  center1 <- drop(crossprod(G1, m1))
  center2 <- drop(crossprod(G2, c2))
  omega_c <- drop(omega1 %*% c2)
  function(rows) {
    rr <- data$r[rows]
    p <- data$pos[rows[rr]]
    tt <- data$t[p, , drop = FALSE]
    h <- data$h[rows, , drop = FALSE]
    hr <- h[rr, , drop = FALSE]
    g1 <- g2 <- -h
    g1[rr, ] <- hr * w1[p]
    g2[rr, ] <- hr * w2[p]
    phi1 <- sweep(sweep(g1, 2L, m1) %*% G1, 2L, center1)
    phi1[rr, ] <- phi1[rr, , drop = FALSE] - tt * (w1[p] * drop(hr %*% m1))
    if1 <- -phi1 %*% t(A1inv)
    phi2 <- sweep(sweep(g2, 2L, m2) %*% t(left), 2L, center2)
    phi2[rr, ] <- phi2[rr, , drop = FALSE] - tt * (w2[p] * drop(hr %*% c2))
    # IF(W) is -W IF(Omega1) W; include both stages' contributions.
    phi2 <- phi2 - (sweep(g1 * drop(g1 %*% c2), 2L, omega_c) + if1 %*% t(contract)) %*% t(left)
    -phi2 %*% t(A2inv)
  }
}

fit_mnar_influence <- function(fit, working = TRUE, rows = NULL) {
  if (!identical(fit$type, "MNAR")) stop("Supply a fit from fit_mnar().")
  evaluator <- make_mnar_influence(fit, working)
  if (!is.null(rows)) return(evaluator(rows))
  do.call(rbind, lapply(cal_blocks(fit$data$n, fit$data$partitions), evaluator))
}

# Epanechnikov density estimator for exact check loss; Gaussian density
# for the fixed-resolution projection. All reported variances are on the
# coefficient-estimate scale (divide the influence covariance by N).
qr_inference <- function(x, y, r, beta, tau, fit = NULL, weights = NULL,
                         bandwidth = NULL, smooth = FALSE, robust = FALSE,
                         partitions = 1L) {
  x <- as.matrix(x)
  r <- as.logical(r)
  n <- length(r)
  nr <- sum(r)
  d <- ncol(x)
  k <- length(tau)
  beta <- matrix(beta, nrow = d)
  if (anyNA(r) || nr == 0L || !nrow(x) %in% c(n, nr) ||
      !length(y) %in% c(n, nr) || ncol(beta) != k ||
      any(!is.finite(tau) | tau <= 0 | tau >= 1))
    stop("Inconsistent QR data, coefficients, or quantile levels.")
  if (is.null(bandwidth)) bandwidth <- if (smooth) .0125 else n^(-.2)
  if (length(bandwidth) != 1L || !is.finite(bandwidth) || bandwidth <= 0)
    stop("Bandwidth must be a positive finite scalar.")
  if (is.null(weights)) weights <- if (is.null(fit)) rep(1, nr) else fit$weights
  if (length(weights) != nr || any(!is.finite(weights) | weights <= 0))
    stop("One positive finite weight per respondent is required.")
  if (!is.null(fit) && (!identical(r, fit$data$r) || !identical(weights, fit$weights)))
    stop("Inference rows and weights must agree with the fitted response model.")
  pos <- integer(n)
  pos[r] <- seq_len(nr)
  q <- if (is.null(fit)) 0L else length(fit$theta)
  if (q > 0L) {
    nuisance <- if (identical(fit$type, "MNAR")) make_mnar_influence(fit, robust) else function(rows) {
      xx <- fit$data$x[rows, , drop = FALSE]
      (xx * (r[rows] - fit$data$p[rows])) %*% t(fit$data$A_inv)
    }
  }

  U_sum <- numeric(d * k)
  eta_sum <- numeric(q)
  UU <- matrix(0, d * k, d * k)
  U_eta <- B <- matrix(0, d * k, q)
  eta_eta <- matrix(0, q, q)
  J <- matrix(0, d * k, d * k)
  for (ids in cal_blocks(n, partitions)) {
    rr <- r[ids]
    ri <- ids[rr]
    p <- pos[ri]
    xx <- x[if (nrow(x) == n) ri else p, , drop = FALSE]
    yy <- y[if (length(y) == n) ri else p]
    if (!all(is.finite(xx)) || !all(is.finite(yy)))
      stop("Observed regression inputs must be finite.")
    U <- matrix(0, length(ids), d * k)
    w <- weights[p]
    if (q > 0L) {
      tt <- if (identical(fit$type, "MNAR")) fit$data$t[p, , drop = FALSE] else
        fit$data$x[ri, , drop = FALSE]
      eta <- nuisance(ids)
      eta_sum <- eta_sum + colSums(eta)
      eta_eta <- eta_eta + crossprod(eta)
    }
    for (j in seq_len(k)) {
      cols <- seq.int((j - 1L) * d + 1L, j * d)
      residual <- yy - drop(xx %*% beta[, j])
      u <- residual / bandwidth
      psi <- if (smooth) tau[j] - pnorm(-u) else tau[j] - (residual <= 0)
      kernel <- if (smooth) dnorm(u) else .75 * pmax(1 - u^2, 0)
      J[cols, cols] <- J[cols, cols] + crossprod(xx, xx * (w * kernel / bandwidth))
      U[rr, cols] <- xx * (w * psi)
      if (q > 0L) B[cols, ] <- B[cols, , drop = FALSE] - crossprod(xx, tt * ((w - 1) * psi))
    }
    U_sum <- U_sum + colSums(U)
    UU <- UU + crossprod(U)
    if (q > 0L) U_eta <- U_eta + crossprod(U, eta)
  }
  J <- J / n
  UU <- UU / n - tcrossprod(U_sum / n)
  V <- UU
  theta_covariance <- NULL
  if (q > 0L) {
    B <- B / n
    U_eta <- U_eta / n - tcrossprod(U_sum / n, eta_sum / n)
    eta_eta <- eta_eta / n - tcrossprod(eta_sum / n)
    mixed <- tcrossprod(U_eta, B)
    V <- UU + mixed + t(mixed) + B %*% eta_eta %*% t(B)
    theta_covariance <- eta_eta / n
  }
  inv <- cal_inverse(J, "The quantile score Jacobian")
  covariance <- inv %*% V %*% inv / n
  covariance <- (covariance + t(covariance)) / 2
  fixed_covariance <- inv %*% UU %*% inv / n
  se <- matrix(sqrt(pmax(diag(covariance), 0)), nrow = d)
  fixed_se <- matrix(sqrt(pmax(diag(fixed_covariance), 0)), nrow = d)
  if (k == 1L) {
    se <- drop(se)
    fixed_se <- drop(fixed_se)
  }
  list(se = se, covariance = covariance, fixed_se = fixed_se,
       theta_covariance = theta_covariance, condition = kappa(J, exact = TRUE))
}
