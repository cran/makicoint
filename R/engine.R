# Internal engine for the Maki (2012) cointegration test.
# General recursive break search; reproduces the original GAUSS/tspdlib
# procedure by default, with an optional Maki (2012) paper break rule.
# All functions here are internal (not exported).

# OLS residual (column-order invariant, so the search may append candidate
# columns in any order).
.mkc_resid <- function(y, X) {
  b <- solve(crossprod(X), crossprod(X, y))
  as.vector(y - X %*% b)
}

# Fixed-break part of the cointegrating design, built once per step.
.mkc_xfixed <- function(datap, n, model, fixed) {
  k <- ncol(datap)
  R <- datap[, -1, drop = FALSE]
  X <- matrix(1, n, 1)
  if (length(fixed) > 0L) {
    for (b in fixed) X <- cbind(X, c(rep.int(0, b), rep.int(1, n - b)))
  }
  if (model == 1L || model == 3L) {
    X <- cbind(X, seq_len(n))
    if (model == 3L && length(fixed) > 0L) {
      for (b in fixed) X <- cbind(X, c(rep.int(0, b), (b + 1):n))
    }
  }
  X <- cbind(X, R)
  if ((model == 2L || model == 3L) && length(fixed) > 0L) {
    for (b in fixed) {
      X <- cbind(X, rbind(matrix(0, b, k - 1), R[(b + 1):n, , drop = FALSE]))
    }
  }
  X
}

# Candidate columns added by a trial break at position i.
.mkc_xcand <- function(datap, n, model, i) {
  k <- ncol(datap)
  dui <- c(rep.int(0, i), rep.int(1, n - i))
  if (model == 0L || model == 1L) return(matrix(dui, ncol = 1))
  R <- datap[, -1, drop = FALSE]
  dxi <- rbind(matrix(0, i, k - 1), R[(i + 1):n, , drop = FALSE])
  if (model == 2L) return(cbind(dui, dxi))
  dtri <- c(rep.int(0, i), (i + 1):n)
  cbind(dui, dtri, dxi)
}

# t-sig (general-to-specific) ADF lag, threshold 1.654.
.mkc_tsig <- function(e, maxlags) {
  n <- length(e); dy <- diff(e)
  p <- maxlags
  while (p >= 1) {
    if (n - 1 - p >= 3) {
      X <- matrix(e[(p + 1):(n - 1)], ncol = 1)
      for (i in seq_len(p)) X <- cbind(X, dy[(p + 1 - i):(n - 1 - i)])
      ly <- dy[(p + 1):(n - 1)]
      XtXi <- solve(crossprod(X))
      b <- XtXi %*% crossprod(X, ly)
      resid <- ly - X %*% b
      s2 <- sum(resid^2) / (nrow(X) - ncol(X))
      se <- sqrt(diag(s2 * XtXi))
      if (abs(b[p + 1] / se[p + 1]) > 1.654) return(p)
    }
    p <- p - 1
  }
  0L
}

# AIC/BIC ADF lag.
.mkc_ic <- function(e, maxlags, aic = TRUE) {
  n <- length(e); dy <- diff(e)
  best <- 0L; bestic <- Inf
  for (p in 0:maxlags) {
    if (n - 1 - p >= 3) {
      X <- matrix(e[(p + 1):(n - 1)], ncol = 1)
      if (p > 0) for (i in seq_len(p)) X <- cbind(X, dy[(p + 1 - i):(n - 1 - i)])
      ly <- dy[(p + 1):(n - 1)]
      b <- solve(crossprod(X), crossprod(X, ly))
      resid <- ly - X %*% b
      no <- nrow(X); np <- ncol(X)
      s2 <- sum(resid^2) / no
      ic <- if (aic) log(s2) + 2 * np / no else log(s2) + log(no) * np / no
      if (ic < bestic) { bestic <- ic; best <- p }
    }
  }
  best
}

# ADF tau on a residual; returns tau, the ADF-regression SSR, and the lag.
.mkc_adf <- function(e, lagopt, maxlags) {
  n <- length(e); dy <- diff(e)
  lag <- switch(lagopt + 1L,
                0L,                          # 0 = zero lags
                .mkc_tsig(e, maxlags),       # 1 = t-sig
                .mkc_ic(e, maxlags, TRUE),   # 2 = aic
                .mkc_ic(e, maxlags, FALSE))  # 3 = bic
  if (lagopt == 9L) lag <- maxlags           # 9 = fixed at maxlags
  r <- 2 + lag
  X <- matrix(e[(r - 1):(n - 1)], ncol = 1)
  if (lag > 0) for (q in seq_len(lag)) X <- cbind(X, dy[(r - 1 - q):(n - 1 - q)])
  ly <- dy[(r - 1):(n - 1)]
  XtXi <- solve(crossprod(X))
  b <- XtXi %*% crossprod(X, ly)
  resid <- ly - X %*% b
  ssr <- sum(resid^2)
  s2 <- ssr / (nrow(X) - ncol(X))
  se <- sqrt(diag(s2 * XtXi))
  list(tau = as.numeric(b[1] / se[1]), ssr = ssr, lag = lag)
}

# Search for one additional break given the fixed set.
# rule: "gauss" (default) chooses the min-tau segment, break at its ADF-SSR
#   minimiser; "paper" chooses the global cointegrating-SSR minimiser.
.mkc_search <- function(datap, n, model, tb, fixed, lagopt, maxlags, rule, gssr) {
  y <- datap[, 1]
  Xf <- .mkc_xfixed(datap, n, model, fixed)
  bounds <- c(0, sort(fixed), n)
  nseg <- length(bounds) - 1L
  segtau <- rep(NA_real_, nseg)
  segpos <- integer(nseg)
  seglag <- integer(nseg)
  gminssr <- Inf; gminpos <- 0L
  for (s in seq_len(nseg)) {
    lo <- bounds[s] + tb + 1L
    hi <- bounds[s + 1] - tb
    if (lo > hi) next
    sbtau <- Inf; sbssr <- Inf; sbpos <- 0L; sblag <- 0L
    for (i in lo:hi) {
      X <- cbind(Xf, .mkc_xcand(datap, n, model, i))
      e <- .mkc_resid(y, X)
      ad <- .mkc_adf(e, lagopt, maxlags)
      ssr <- if (gssr) ad$ssr else sum(e^2)
      if (ad$tau < sbtau) { sbtau <- ad$tau; sblag <- ad$lag }
      if (ssr < sbssr) { sbssr <- ssr; sbpos <- i }
      if (ssr < gminssr) { gminssr <- ssr; gminpos <- i }
    }
    segtau[s] <- sbtau; segpos[s] <- sbpos; seglag[s] <- sblag
  }
  if (all(is.na(segtau))) return(list(tau = NA_real_, bp = 0L, lag = 0L))
  smin <- which.min(segtau)
  bp <- if (identical(rule, "gauss")) segpos[smin] else gminpos
  list(tau = segtau[smin], bp = bp, lag = seglag[smin])
}

# Sequential test: accumulate the minimum tau across steps.
.mkc_test <- function(datap, m, model, tb, lagopt, maxlags, rule, gssr) {
  n <- nrow(datap)
  fixed <- integer(0)
  alltau <- rep(NA_real_, m)
  runmin <- Inf; lagused <- 0L
  for (j in seq_len(m)) {
    res <- .mkc_search(datap, n, model, tb, fixed, lagopt, maxlags, rule, gssr)
    alltau[j] <- res$tau
    if (length(res$bp) == 1L && res$bp > 0L) fixed <- sort(c(fixed, res$bp))
    if (!is.na(res$tau) && res$tau < runmin) { runmin <- res$tau; lagused <- res$lag }
  }
  list(tau = min(alltau, na.rm = TRUE), bp = fixed, lag = lagused)
}

# Maki (2012) Table 1 critical values: model 0-3, k = 1-4 regressors, m = 1-5,
# columns 1%, 5%, 10%. Returns c(NA, NA, NA) outside the tabulated range.
.mkc_cvtable <- function(k, m, model) {
  if (m < 1 || m > 5 || k < 1 || k > 4) return(c(NA_real_, NA_real_, NA_real_))
  tab <- list(
    "0" = list(
      rbind(c(-5.709,-4.602,-4.354),c(-5.416,-4.892,-4.610),c(-5.563,-5.083,-4.784),c(-5.776,-5.230,-4.982),c(-5.959,-5.426,-5.131)),
      rbind(c(-5.541,-5.004,-4.733),c(-5.717,-5.211,-4.957),c(-5.943,-5.392,-5.125),c(-6.075,-5.550,-5.297),c(-6.296,-5.760,-5.491)),
      rbind(c(-5.820,-5.341,-5.101),c(-5.984,-5.517,-5.272),c(-6.229,-5.704,-5.427),c(-6.406,-5.871,-5.603),c(-6.555,-6.038,-5.773)),
      rbind(c(-6.139,-5.650,-5.386),c(-6.303,-5.839,-5.575),c(-6.501,-5.992,-5.714),c(-6.640,-6.132,-5.892),c(-6.856,-6.306,-6.039))),
    "1" = list(
      rbind(c(-5.524,-5.038,-4.784),c(-5.708,-5.196,-4.938),c(-5.833,-5.373,-5.106),c(-6.059,-5.508,-5.245),c(-6.193,-5.699,-5.449)),
      rbind(c(-5.840,-5.359,-5.117),c(-6.011,-5.518,-5.247),c(-6.169,-5.691,-5.408),c(-6.329,-5.831,-5.558),c(-6.530,-5.993,-5.722)),
      rbind(c(-6.144,-5.645,-5.398),c(-6.271,-5.796,-5.538),c(-6.472,-5.957,-5.682),c(-6.575,-6.086,-5.820),c(-6.784,-6.250,-5.976)),
      rbind(c(-6.361,-5.913,-5.686),c(-6.556,-6.055,-5.805),c(-6.741,-6.214,-5.974),c(-6.845,-6.373,-6.096),c(-7.053,-6.494,-6.220))),
    "2" = list(
      rbind(c(-5.457,-4.895,-4.626),c(-5.863,-5.363,-5.070),c(-6.251,-5.703,-5.402),c(-6.596,-6.011,-5.723),c(-6.915,-6.357,-6.057)),
      rbind(c(-6.020,-5.558,-5.287),c(-6.628,-6.093,-5.833),c(-7.031,-6.516,-6.210),c(-7.470,-6.872,-6.563),c(-7.839,-7.288,-6.976)),
      rbind(c(-6.565,-6.035,-5.773),c(-7.232,-6.702,-6.411),c(-7.767,-7.155,-6.868),c(-8.236,-7.625,-7.329),c(-8.673,-8.110,-7.796)),
      rbind(c(-7.021,-6.520,-6.242),c(-7.756,-7.244,-6.964),c(-8.336,-7.803,-7.481),c(-8.895,-8.292,-8.004),c(-9.441,-8.869,-8.541))),
    "3" = list(
      rbind(c(-6.048,-5.541,-5.281),c(-6.620,-6.100,-5.845),c(-7.082,-6.524,-6.267),c(-7.553,-7.009,-6.712),c(-8.004,-7.414,-7.110)),
      rbind(c(-6.523,-6.055,-5.795),c(-7.153,-6.657,-6.397),c(-7.673,-7.145,-6.873),c(-8.217,-7.636,-7.341),c(-8.713,-8.129,-7.811)),
      rbind(c(-6.964,-6.464,-6.220),c(-7.737,-7.201,-6.926),c(-8.331,-7.743,-7.449),c(-8.851,-8.269,-7.960),c(-9.428,-8.800,-8.508)),
      rbind(c(-7.400,-6.911,-6.649),c(-8.167,-7.638,-7.381),c(-8.865,-8.254,-7.977),c(-9.433,-8.871,-8.574),c(-10.08,-9.482,-9.151))))
  as.numeric(tab[[as.character(model)]][[k]][m, ])
}

# Monte-Carlo critical values for any (model, k, m), Maki (2012) design:
# independent driftless random walks under the null.
.mkc_simcv <- function(kreg, m, model, trimm, lagopt, maxlags, rule, gssr,
                       simt, reps, seed) {
  set.seed(seed)
  tb <- max(1L, round(trimm * simt))
  stats <- numeric(reps)
  for (r in seq_len(reps)) {
    datap <- apply(matrix(rnorm(simt * (kreg + 1)), simt, kreg + 1), 2, cumsum)
    stats[r] <- .mkc_test(datap, m, model, tb, lagopt, maxlags, rule, gssr)$tau
  }
  stats <- sort(stats)
  idx <- function(p) max(1L, ceiling(p * reps))
  c(stats[idx(0.01)], stats[idx(0.05)], stats[idx(0.10)])
}

# Feasible maximum number of breaks for given n and trimming.
.mkc_maxfeasible <- function(n, trimm) {
  tb <- max(1L, round(trimm * n))
  floor(n / (tb + 1)) - 1L
}
