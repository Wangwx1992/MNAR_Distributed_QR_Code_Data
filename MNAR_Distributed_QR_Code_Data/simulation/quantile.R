# Independent R implementations of the weighted quantile objectives in the
# manuscript. The algorithms operate on respondent arrays and normalize
# reported losses by the number of records before response masking.

check_loss <- function(r, tau) r * (tau - (r < 0))

qr_norm <- function(x) sqrt(sum(x * x))

qr_project <- function(x, radius) {
  norm <- qr_norm(x)
  if (norm > radius) x * (radius / norm) else x
}

qr_inputs <- function(x, y, w, tau, n_total) {
  x <- as.matrix(x)
  if (!length(y) || nrow(x) != length(y) || length(w) != length(y) ||
      !ncol(x) || !all(is.finite(x)) || !all(is.finite(y)) ||
      !all(is.finite(w) & w > 0) || !is.finite(tau) || tau <= 0 ||
      tau >= 1 || !is.finite(n_total) || n_total < length(y)) {
    stop("Expected a finite respondent design, finite responses, positive weights, and 0 < tau < 1.")
  }
  x
}

qr_blocks <- function(n, partitions) {
  if (length(partitions) != 1L || !is.finite(partitions) ||
      partitions < 1L || partitions > n || partitions != as.integer(partitions)) {
    stop("The number of respondent blocks must be an integer between 1 and the respondent count.")
  }
  edges <- floor(seq(0, n, length.out = partitions + 1L))
  lapply(seq_len(partitions), function(j) seq.int(edges[j] + 1L, edges[j + 1L]))
}

qr_clock <- function() as.numeric(Sys.time())

qr_chol_solve <- function(factor, rhs) {
  drop(backsolve(factor, forwardsolve(t(factor), rhs)))
}

# Primal check loss and the compact-ball dual certificate. Each site's dual
# contribution is only a scalar and a vector of length ncol(x).
qr_ball_gap <- function(x, y, beta, dual, tau, n_total, radius) {
  d <- pmin(tau, pmax(tau - 1, dual))
  objective <- sum(check_loss(y - drop(x %*% beta), tau)) / n_total
  bound <- (sum(y * d) - radius * qr_norm(crossprod(x, d))) / n_total
  list(objective = objective, lower_bound = bound,
       gap = max(0, objective - bound))
}

# Centralized interior point algorithm on the dual quantile LP. The lower
# and upper box slacks are d-(tau-1) and tau-d. Eliminating the two multiplier
# directions yields a p-dimensional symmetric Schur complement. A
# predictor-corrector step supplies the centering parameter and second-order
# complementarity correction.
qr_dual_ip <- function(x, y, tau, n_total, radius, tolerance, maxit = 100L) {
  n <- length(y)
  p <- ncol(x)
  beta <- numeric(p)
  d <- numeric(n)
  lower <- tau - 1
  upper <- tau
  zl <- 1 + pmax(-y, 0)
  zu <- 1 + pmax(y, 0)
  diagnostics <- qr_ball_gap(x, y, beta, d, tau, n_total, radius)
  converged <- diagnostics$gap <= tolerance
  iteration <- 0L
  fraction <- function(slack, direction) {
    decreasing <- direction < 0
    if (any(decreasing)) min(-slack[decreasing] / direction[decreasing]) else Inf
  }

  while (!converged && iteration < maxit) {
    sl <- d - lower
    su <- upper - d
    if (any(sl <= 0 | su <= 0 | zl <= 0 | zu <= 0)) break
    residual <- y - drop(x %*% beta) + zl - zu
    equality <- drop(crossprod(x, d))
    diagonal <- 1 / (zl / sl + zu / su)
    schur <- crossprod(sqrt(diagonal) * x)
    factor <- tryCatch(chol(schur), error = function(e) NULL)
    if (is.null(factor)) break

    direction <- function(cl, cu) {
      h <- diagonal * (residual + cl / sl - cu / su)
      db <- qr_chol_solve(factor, drop(crossprod(x, h)) + equality)
      dd <- h - diagonal * drop(x %*% db)
      list(beta = db, d = dd, zl = (cl - zl * dd) / sl,
           zu = (cu + zu * dd) / su)
    }
    step_lengths <- function(step) {
      c(min(fraction(sl, step$d), fraction(su, -step$d)),
        min(fraction(zl, step$zl), fraction(zu, step$zu)))
    }
    affine <- direction(-sl * zl, -su * zu)
    alpha <- pmin(1, step_lengths(affine))
    mu <- mean(sl * zl + su * zu) / 2
    affine_mu <- mean(
      (sl + alpha[1] * affine$d) * (zl + alpha[2] * affine$zl) +
      (su - alpha[1] * affine$d) * (zu + alpha[2] * affine$zu)
    ) / 2
    centering <- min(1, max(0, affine_mu / mu))^3 * mu
    corrector <- direction(centering - sl * zl - affine$d * affine$zl,
                           centering - su * zu + affine$d * affine$zu)
    alpha <- pmin(1, 0.995 * step_lengths(corrector))
    d <- d + alpha[1] * corrector$d
    beta <- beta + alpha[2] * corrector$beta
    zl <- zl + alpha[2] * corrector$zl
    zu <- zu + alpha[2] * corrector$zu
    iteration <- iteration + 1L
    if (!all(is.finite(beta))) break
    diagnostics <- qr_ball_gap(x, y, qr_project(beta, radius), d, tau, n_total, radius)
    converged <- qr_norm(beta) <= radius && diagnostics$gap <= tolerance
  }

  beta <- qr_project(beta, radius)
  c(list(beta = beta, converged = converged, iterations = iteration,
         backend = "base R primal-dual interior point"), diagnostics)
}

fit_qr_central <- function(x, y, w, tau, n_total = length(y), radius = 20,
                           tolerance = min(1e-8, n_total^(-1.5))) {
  x <- qr_inputs(x, y, w, tau, n_total)
  x_weighted <- w * x
  y_weighted <- w * y

  qr_dual_ip(x_weighted, y_weighted, tau, n_total, radius, tolerance)
}

fit_qr_admm <- function(x, y, w, tau, partitions, n_total = length(y),
                        tolerance = .05, stopping = "coefficient", maxit = NULL,
                        radius = 20) {
  started <- qr_clock()
  x <- qr_inputs(x, y, w, tau, n_total)
  stopping <- match.arg(stopping, c("coefficient", "residual"))
  if (is.null(maxit)) maxit <- if (stopping == "coefficient") 5000L else 20000L
  if (missing(tolerance) && stopping == "residual") tolerance <- 1e-4
  if (maxit < 1L || tolerance <= 0 || radius <= 0) stop("Expected positive iteration limit, tolerance, and radius.")

  blocks <- qr_blocks(length(y), partitions)
  n_site <- lengths(blocks)
  p <- ncol(x)
  kappa <- if (stopping == "residual") n_site else rep(1, partitions)
  setup <- numeric(partitions)
  sites <- lapply(seq_along(blocks), function(j) {
    tick <- qr_clock()
    rows <- blocks[[j]]
    design <- w[rows] * x[rows, , drop = FALSE]
    response <- w[rows] * y[rows]
    result <- list(x = design, y = response,
                   factor = chol(crossprod(design) + diag(kappa[j], p)),
                   v = numeric(length(rows)), remainder = numeric(length(rows)))
    setup[j] <<- qr_clock() - tick
    result
  })

  ideal_elapsed <- max(setup)
  beta <- numeric(p)
  checked <- beta
  local_beta <- matrix(0, p, partitions)
  local_u <- local_beta
  converged <- FALSE
  site_elapsed <- numeric(partitions)
  criterion <- Inf

  for (iteration in seq_len(maxit)) {
    tick <- qr_clock()
    previous <- beta
    beta <- qr_project(drop((local_beta + local_u) %*% kappa) / sum(kappa), radius)
    coordination <- qr_clock() - tick
    for (j in seq_len(partitions)) {
      tick <- qr_clock()
      site <- sites[[j]]
      arg <- site$y - drop(site$x %*% local_beta[, j]) - site$v
      r <- pmax(arg - tau, 0) - pmax(-arg - (1 - tau), 0)
      rhs <- kappa[j] * (beta - local_u[, j]) +
        drop(crossprod(site$x, site$y - r - site$v))
      local_beta[, j] <- qr_chol_solve(site$factor, rhs)
      local_u[, j] <- local_u[, j] + local_beta[, j] - beta
      site$remainder <- drop(site$x %*% local_beta[, j]) + r - site$y
      site$v <- site$v + site$remainder
      sites[[j]] <- site
      site_elapsed[j] <- qr_clock() - tick
    }
    ideal_elapsed <- ideal_elapsed + coordination + max(site_elapsed)

    tick <- qr_clock()
    if (stopping == "coefficient") {
      criterion <- sum(abs(beta - previous))
      converged <- criterion > 0 && criterion <= tolerance
    } else if (iteration %% 25L == 0L) {
      consensus <- sqrt(sum(n_site * colSums((local_beta - beta)^2)) / length(y))
      diagnostics_elapsed <- numeric(partitions)
      residual_sum <- vapply(seq_along(sites), function(j) {
        site_tick <- qr_clock()
        value <- sum(sites[[j]]$remainder^2)
        diagnostics_elapsed[j] <<- qr_clock() - site_tick
        value
      }, 0)
      residual <- sqrt(sum(residual_sum) / length(y))
      criterion <- max(consensus, residual,
                       qr_norm(beta - checked) / (1 + qr_norm(beta)))
      checked <- beta
      converged <- criterion <= tolerance
      # Replace measured serial site diagnostics with their modeled critical path.
      ideal_elapsed <- ideal_elapsed + max(diagnostics_elapsed) - sum(diagnostics_elapsed)
    }
    ideal_elapsed <- ideal_elapsed + qr_clock() - tick
    if (converged) break
  }

  pieces <- lapply(seq_along(sites), function(j) {
    tick <- qr_clock()
    site <- sites[[j]]
    d <- pmin(tau, pmax(tau - 1, -site$v))
    part <- list(loss = sum(check_loss(site$y - drop(site$x %*% beta), tau)),
                 linear = sum(site$y * d), gradient = drop(crossprod(site$x, d)))
    site_elapsed[j] <<- qr_clock() - tick
    part
  })
  tick <- qr_clock()
  objective <- sum(vapply(pieces, function(part) part$loss, 0)) / n_total
  lower_bound <- (sum(vapply(pieces, function(part) part$linear, 0)) -
    radius * qr_norm(Reduce("+", lapply(pieces, function(part) part$gradient)))) / n_total
  ideal_elapsed <- ideal_elapsed + max(site_elapsed) + qr_clock() - tick
  list(beta = beta, objective = objective, lower_bound = lower_bound,
       gap = max(0, objective - lower_bound), converged = converged,
       iterations = iteration, criterion = criterion,
       parallel_seconds = ideal_elapsed, serial_seconds = qr_clock() - started)
}

qr_gaussian_loss <- function(r, tau, bandwidth) {
  r * (tau - pnorm(-r / bandwidth)) + bandwidth * dnorm(r / bandwidth)
}

fit_qr_renewable <- function(x, y, w, tau, partitions, n_total = length(y),
                             tolerance = 1e-4, maxit = 5000L) {
  started <- qr_clock()
  x <- qr_inputs(x, y, w, tau, n_total)
  blocks <- qr_blocks(length(y), partitions)
  if (length(blocks[[1L]]) < 2L) stop("The first renewable block needs at least two respondents.")
  p <- ncol(x)
  beta <- numeric(p)
  accumulated <- matrix(0, p, p)
  count <- 0L
  scores <- numeric(partitions)
  iterations <- integer(partitions)

  for (j in seq_along(blocks)) {
    rows <- blocks[[j]]
    count <- count + length(rows)
    h <- (count * log(count))^(-1 / 4)
    design <- w[rows] * x[rows, , drop = FALSE]
    response <- w[rows] * y[rows]
    previous <- beta

    evaluate <- function(candidate, derivative = TRUE) {
      residual <- response - drop(design %*% candidate)
      displacement <- candidate - previous
      penalty_gradient <- drop(accumulated %*% displacement)
      value <- (sum(qr_gaussian_loss(residual, tau, h)) +
        sum(displacement * penalty_gradient) / 2) / count
      if (!derivative) return(value)
      gradient <- (drop(crossprod(design, pnorm(-residual / h) - tau)) +
        penalty_gradient) / count
      list(value = value, gradient = gradient)
    }

    current <- evaluate(beta)
    step <- 1
    for (i in seq_len(maxit)) {
      g <- current$gradient
      if (max(abs(g)) <= tolerance) break
      for (trial in seq_len(30L)) {
        proposal <- beta - step * g
        if (evaluate(proposal, FALSE) <= current$value - 1e-4 * step * sum(g^2)) break
        if (trial < 30L) step <- step / 2
      }
      next_state <- evaluate(proposal)
      s <- proposal - beta
      delta <- next_state$gradient - g
      curvature <- sum(s * delta)
      step <- if (curvature > 0) min(sum(s^2) / curvature,
                                    curvature / sum(delta^2), 100) else 1
      current <- next_state
      beta <- proposal
      iterations[j] <- i
    }
    scores[j] <- max(abs(current$gradient))
    residual <- response - drop(design %*% beta)
    accumulated <- accumulated + crossprod(sqrt(dnorm(residual / h) / h) * design)
  }

  objective <- sum(w * check_loss(y - drop(x %*% beta), tau)) / n_total
  list(beta = beta, objective = objective, converged = all(scores <= tolerance),
       iterations = sum(iterations), block_iterations = iterations,
       block_scores = scores, serial_seconds = qr_clock() - started)
}

fit_qr_smooth <- function(x, y, w, tau, n_total = length(y), a = .0125,
                          start = NULL, radius = 20, partitions = 1L,
                          tolerance = 2e-11, maxit = 100L) {
  x <- qr_inputs(x, y, w, tau, n_total)
  if (!is.finite(a) || a <= 0) stop("The smoothing resolution must be positive.")
  p <- ncol(x)
  sites <- if (partitions == 1L) {
    list(list(x = x, y = y, w = w))
  } else {
    lapply(qr_blocks(length(y), partitions), function(rows) {
      list(x = x[rows, , drop = FALSE], y = y[rows], w = w[rows])
    })
  }

  evaluate <- function(beta, derivatives = TRUE) {
    value <- 0
    gradient <- numeric(p)
    hessian <- matrix(0, p, p)
    for (site in sites) {
      r <- site$y - drop(site$x %*% beta)
      z <- r / a
      score <- pnorm(-z) - tau
      density <- dnorm(z)
      value <- value + sum(site$w * (-r * score + a * density))
      if (derivatives) {
        gradient <- gradient + drop(crossprod(site$x, site$w * score))
        hessian <- hessian + crossprod(sqrt(site$w * density / a) * site$x)
      }
    }
    if (!derivatives) return(value / n_total)
    list(value = value / n_total, gradient = gradient / n_total,
         hessian = hessian / n_total)
  }

  if (is.null(start)) {
    cross <- Reduce("+", lapply(sites, function(site) crossprod(site$x)))
    xy <- Reduce("+", lapply(sites, function(site) crossprod(site$x, site$y)))
    start <- tryCatch(drop(solve(cross, xy)), error = function(e) numeric(p))
    if (all(abs(x[, 1L] - 1) < 1e-12)) {
      start[1L] <- start[1L] + unname(quantile(y - drop(x %*% start), tau))
    }
  }

  beta <- qr_project(start, radius)
  current <- evaluate(beta)
  converged <- FALSE
  iterations <- 0L
  for (iteration in seq_len(maxit)) {
    gradient <- current$gradient
    gap <- max(0, sum(gradient * beta) + radius * qr_norm(gradient))
    if (gap <= tolerance) {
      converged <- TRUE
      break
    }
    spectral <- eigen(current$hessian, symmetric = TRUE)
    eigenvalues <- pmax(spectral$values, .Machine$double.eps *
                         max(1, max(spectral$values)))
    rhs <- drop(crossprod(spectral$vectors, current$hessian %*% beta - gradient))
    unconstrained <- rhs / eigenvalues
    if (qr_norm(unconstrained) > radius) {
      root <- uniroot(function(lambda) qr_norm(rhs / (eigenvalues + lambda)) - radius,
                      c(0, max(1, qr_norm(rhs) / radius)), tol = 1e-13)$root
      unconstrained <- rhs / (eigenvalues + root)
    }
    direction <- qr_project(drop(spectral$vectors %*% unconstrained), radius) - beta
    slope <- sum(gradient * direction)
    if (slope >= 0) break
    step <- 1
    accepted <- FALSE
    for (trial in seq_len(40L)) {
      proposal <- beta + step * direction
      if (evaluate(proposal, FALSE) <= current$value + 1e-4 * step * slope +
          8 * .Machine$double.eps * max(1, abs(current$value))) {
        accepted <- TRUE
        break
      }
      step <- step / 2
    }
    if (!accepted) break
    beta <- proposal
    current <- evaluate(beta)
    iterations <- iteration
  }
  gap <- max(0, sum(current$gradient * beta) + radius * qr_norm(current$gradient))
  list(beta = beta, objective = current$value, gradient = current$gradient,
       hessian = current$hessian, gap = gap,
       converged = converged || gap <= tolerance, iterations = iterations)
}
