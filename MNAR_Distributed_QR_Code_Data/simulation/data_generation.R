# Simulation designs in Section 5. This file can be sourced on its own.
# R uses its own RNG, so replication indices reproduce the R sample stream.
set_design_seed <- function(key) {
  seed <- as.integer(as.double(key) %% 2147483646 + 1)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
}

primary_designs <- function(smoke = FALSE) {
  sizes <- if (smoke) c(800L, 1200L) else c(10000L, 50000L)
  missing <- if (smoke) .4 else c(.2, .4)
  grid <- expand.grid(missing = missing, error = c("normal", "t3", "lognormal"),
                      N = sizes, stringsAsFactors = FALSE)
  grid$index <- seq_len(nrow(grid)) - 1L
  grid$scenario <- sprintf("N%d_%s_m%02d", grid$N, grid$error, round(100 * grid$missing))
  grid
}

primary_seed <- function(index, replication, purpose) {
  2026103107 + 1000003 * index + 1091 * replication + purpose
}

primary_outcome <- function(n, error) {
  corr <- .35 ^ abs(outer(seq_len(4), seq_len(4), "-"))
  u <- 2 * pnorm(matrix(rnorm(n * 4), n, 4) %*% chol(corr)) - 1
  e <- switch(error, normal = rnorm(n), t3 = rt(n, df = 3) / sqrt(3),
              lognormal = (exp(rnorm(n)) - exp(.5)) / sqrt(exp(1) * expm1(1)),
              stop("Unknown error distribution: ", error))
  y <- 1 + u[, 1L] - u[, 2L] + .5 * u[, 3L] + .5 * u[, 4L] + (1 + .4 * u[, 1L]) * e
  list(x = cbind(1, u), y = y)
}

primary_truth <- function(tau, error) {
  q <- switch(error, normal = qnorm(tau), t3 = qt(tau, df = 3) / sqrt(3),
              lognormal = (exp(qnorm(tau)) - exp(.5)) / sqrt(exp(1) * expm1(1)),
              stop("Unknown error distribution: ", error))
  c(1 + q, 1 + .4 * q, -1, .5, .5)
}

primary_transform <- function(y) tanh((y - 1) / 2)

primary_moments <- function(x, z) {
  cbind(x[, 1:3, drop = FALSE], z, z ^ 2, x[, 2L] * z, x[, 3L] * z)
}

primary_reference <- function(design, smoke = FALSE) {
  set_design_seed(primary_seed(design$index, 0, 1))
  dat <- primary_outcome(if (smoke) 2^13 else 2^20, design$error)
  v <- .3 * dat$x[, 2L] - .3 * dat$x[, 3L] + primary_transform(dat$y)
  intercept <- uniroot(function(a) mean(plogis(a + v)) - (1 - design$missing),
                       interval = c(-4, 4), tol = 1e-13)$root
  set_design_seed(primary_seed(design$index, 0, 2))
  dat <- primary_outcome(if (smoke) 4000L else 100000L, design$error)
  z <- primary_transform(dat$y) + runif(length(dat$y), -sqrt(3) / 2, sqrt(3) / 2)
  h <- primary_moments(dat$x, z)
  centered <- sweep(h, 2L, colMeans(h), "-")
  scale <- sqrt(colMeans(centered ^ 2))
  scale[1L] <- 1
  list(intercept = intercept, scale = scale)
}

generate_primary <- function(design, replication, reference = primary_reference(design)) {
  set_design_seed(primary_seed(design$index, replication, 10))
  dat <- primary_outcome(design$N, design$error)
  x <- dat$x
  s <- primary_transform(dat$y)
  z <- s + runif(design$N, -sqrt(3) / 2, sqrt(3) / 2)
  pi <- plogis(reference$intercept + .3 * x[, 2L] - .3 * x[, 3L] + s)
  r <- runif(design$N) < pi
  h <- sweep(primary_moments(x, z), 2L, reference$scale, "/")
  # Randomly distribute records across the twenty calibration partitions.
  set_design_seed(primary_seed(design$index, replication, 11))
  j <- sample.int(design$N)
  list(x = x[j, , drop = FALSE], y = dat$y[j], r = r[j], pi = pi[j],
       z = z[j], h = h[j, , drop = FALSE])
}

computing_designs <- function(smoke = FALSE) {
  pairs <- if (smoke) {
    data.frame(N = c(800L, 800L, 1200L), p = c(10L, 20L, 20L))
  } else {
    data.frame(N = c(50000L, 50000L, 200000L), p = c(50L, 100L, 100L))
  }
  grid <- rbind(transform(pairs, hetero = FALSE), transform(pairs, hetero = TRUE))
  grid$scenario <- sprintf("N%d_p%d_%s", grid$N, grid$p, ifelse(grid$hetero, "hetero", "homo"))
  grid
}

computing_x1 <- function(n) {
  lower <- pnorm(-4.5)
  qnorm(lower + (1 - lower) * runif(n))
}

computing_location_scale <- function(p, hetero) {
  nc <- as.integer(.8 * p)
  nb <- p - nc
  mu_x1 <- dnorm(-4.5) / pnorm(-4.5, lower.tail = FALSE)
  var_x1 <- 1 - 4.5 * mu_x1 - mu_x1 ^ 2
  error_variance <- 25
  if (hetero) error_variance <- error_variance + var_x1 + mu_x1 ^ 2 + 10 * mu_x1
  c(center = -3 + .5 * nb + mu_x1,
    scale = sqrt(nc - 1 + var_x1 + .25 * nb + error_variance))
}

computing_reference <- function(p, hetero, smoke = FALSE) {
  set_design_seed(120261107 + 10 * p + as.integer(hetero))
  n <- if (smoke) 8192L else 262144L
  nc <- as.integer(.8 * p)
  first <- matrix(rnorm(n * 5L), n, 5L)
  first[, 1L] <- computing_x1(n)
  y <- -3 + drop(first %*% c(1, -1, 1, -1, 1)) +
    rnorm(n, sd = sqrt(nc - 5L)) + rbinom(n, p - nc, .5) +
    (if (hetero) 5 + first[, 1L] else 5) * rnorm(n)
  pars <- computing_location_scale(p, hetero)
  predictor <- drop(first %*% c(.1, -.1, .1, -.1, .1)) +
    tanh((y - pars[["center"]]) / pars[["scale"]])
  intercept <- uniroot(function(a) mean(plogis(a + predictor)) - .75,
                       interval = c(-8, 8), tol = 1e-13)$root
  c(intercept = intercept, pars)
}

generate_computing <- function(design, replication,
                               reference = computing_reference(design$p, design$hetero)) {
  n <- design$N
  p <- design$p
  hetero <- design$hetero
  key <- 220261107 + as.double(n) * 10000 + p * 100 + as.integer(hetero) * 1000000000 + replication
  set_design_seed(key)
  nc <- as.integer(.8 * p)
  continuous <- matrix(rnorm(n * nc), n, nc)
  continuous[, 1L] <- computing_x1(n)
  binary <- matrix(rbinom(n * (p - nc), 1, .5), n, p - nc)
  x <- cbind(1, continuous, binary)
  beta <- c(-3, rep_len(c(1, -1), nc), rep(1, p - nc))
  y <- drop(x %*% beta) + (if (hetero) 5 + x[, 2L] else 5) * rnorm(n)
  s <- tanh((y - reference[["center"]]) / reference[["scale"]])
  z <- s + runif(n, -sqrt(3) / 2, sqrt(3) / 2)
  v <- x[, 1:6, drop = FALSE]
  pi <- plogis(reference[["intercept"]] + drop(v[, -1L, drop = FALSE] %*%
               c(.1, -.1, .1, -.1, .1)) + s)
  r <- runif(n) < pi
  h <- cbind(v, z, z ^ 2, x[, 2L] * z, x[, 3L] * z)
  t <- cbind(v[r, , drop = FALSE], s[r])
  list(x = x[r, , drop = FALSE], y = y[r], r = r, t = t, h = h, beta = beta)
}

computing_truth <- function(base_beta, tau, hetero) {
  q <- qnorm(tau)
  base_beta[1L] <- base_beta[1L] + 5 * q
  if (hetero) base_beta[2L] <- base_beta[2L] + q
  base_beta
}
