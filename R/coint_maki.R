#' Critical Values for the Maki (2012) Cointegration Test
#'
#' Returns the asymptotic critical values from Maki (2012), Table 1, for a
#' given number of regressors, number of breaks, and model.
#'
#' @param k Number of regressors (independent variables), an integer in
#'   \code{1:4}. The Maki (2012) table covers up to four regressors.
#' @param m Number of breaks, an integer in \code{1:5}.
#' @param model Model specification, an integer in \code{0:3}:
#'   \itemize{
#'     \item 0: level shift
#'     \item 1: level shift with trend
#'     \item 2: regime shift
#'     \item 3: regime shift with trend
#'   }
#' @return A numeric vector of length three with the 1\%, 5\% and 10\% critical
#'   values, or \code{c(NA, NA, NA)} outside the tabulated range.
#' @references Maki, D. (2012). Tests for cointegration allowing for an unknown
#'   number of breaks. \emph{Economic Modelling}, 29, 2011-2015.
#'   \doi{10.1016/j.econmod.2012.04.022}
#' @export
#' @examples
#' # Regime-shift model, two regressors, three breaks
#' cv_coint_maki(k = 2, m = 3, model = 2)
cv_coint_maki <- function(k, m, model) {
  .mkc_cvtable(k, m, model)
}

#' Maki (2012) Cointegration Test with Multiple Structural Breaks
#'
#' Performs the Maki (2012) residual-based test for cointegration allowing for
#' an unknown number of structural breaks, up to a user-set maximum. Breaks are
#' found by the sequential Bai-Perron procedure and the cointegrating residual
#' is tested for a unit root with an augmented Dickey-Fuller (ADF) regression;
#' the test statistic is the minimum ADF t-statistic over all candidate breaks.
#'
#' By default the break dates follow the original 'GAUSS'/'tspdlib' rule (the
#' segment with the smallest minimum ADF t-statistic, break at the ADF-regression
#' sum-of-squared-residuals minimiser). With \code{engine = "paper"} each break
#' is the global minimiser of the cointegrating-regression sum of squared
#' residuals, the rule stated in Maki (2012, Steps 2 and 4). The two engines give
#' the same test statistic but may place the breaks differently for two or more
#' breaks.
#'
#' Maki (2012, Table 1) tabulates critical values for one to five breaks. For
#' more than five breaks (his footnote 5 notes the test can use more), set
#' \code{simcv} to simulate the critical values by his Monte-Carlo design.
#'
#' @param data A matrix or data frame with the dependent variable in the first
#'   column and one to four regressors in the remaining columns.
#' @param m Maximum number of structural breaks (integer \code{>= 1}), subject to
#'   feasibility for the sample size and \code{trimm}.
#' @param model Model specification (\code{0:3}); see \code{\link{cv_coint_maki}}.
#' @param trimm Trimming fraction in \code{(0, 0.5)}; the minimum regime length as
#'   a fraction of the sample. Default \code{0.10}.
#' @param lagmethod ADF lag rule: \code{"tsig"} (default, general-to-specific
#'   t-significance), \code{"fixed"} (use \code{maxlags}), \code{"zero"},
#'   \code{"aic"} or \code{"bic"}.
#' @param maxlags Maximum ADF lag order. Default \code{12}.
#' @param engine \code{"gauss"} (default) or \code{"paper"}; see Details.
#' @param simcv Number of Monte-Carlo replications for simulated critical values.
#'   \code{0} (default) uses the Table 1 values for \code{m <= 5}. Required for
#'   inference when \code{m > 5}.
#' @param simt Series length used by the simulation. Default \code{1000}.
#' @param simseed Random seed for the simulation. Default \code{12345}.
#'
#' @return An object of class \code{"maki_test"}: a list with the test
#'   \code{statistic}, \code{breakpoints} (observation indices),
#'   \code{break_fractions}, \code{critical_values}, \code{cv_source}
#'   (\code{"table"}, \code{"simulated"} or \code{"none"}), \code{reject_1},
#'   \code{reject_5}, \code{reject_10}, \code{conclusion}, \code{m}, \code{model},
#'   \code{n}, \code{k}, \code{lags}, \code{engine} and related fields.
#'
#' @references Maki, D. (2012). Tests for cointegration allowing for an unknown
#'   number of breaks. \emph{Economic Modelling}, 29, 2011-2015.
#'   \doi{10.1016/j.econmod.2012.04.022}
#'
#'   Bai, J. and Perron, P. (1998). Estimating and testing linear models with
#'   multiple structural changes. \emph{Econometrica}, 66, 47-78.
#'   \doi{10.2307/2998540}
#'
#' @importFrom stats rnorm
#' @export
#' @examples
#' # A cointegrated pair with a level break at observation 30
#' set.seed(1)
#' n <- 60
#' x <- cumsum(rnorm(n))
#' y <- 2 + 3 * x + 4 * (seq_len(n) > 30) + rnorm(n)
#' res <- coint_maki(cbind(y, x), m = 1, model = 2)
#' res
#'
#' \donttest{
#' # Two breaks, regime-shift model
#' coint_maki(cbind(y, x), m = 2, model = 2)
#' }
coint_maki <- function(data, m = 2, model = 2, trimm = 0.10,
                       lagmethod = c("tsig", "fixed", "zero", "aic", "bic"),
                       maxlags = 12, engine = c("gauss", "paper"),
                       simcv = 0, simt = 1000, simseed = 12345) {
  if (!is.matrix(data) && !is.data.frame(data))
    stop("'data' must be a matrix or data frame.")
  datap <- as.matrix(data)
  if (!is.numeric(datap)) stop("'data' must be numeric.")
  if (anyNA(datap)) stop("'data' must not contain missing values.")
  n <- nrow(datap)
  k <- ncol(datap) - 1L
  if (k < 1L) stop("'data' must have at least two columns (dependent + regressor).")
  if (k > 4L) stop("At most four regressors are supported (Maki 2012, Table 1).")
  if (m < 1L) stop("'m' must be at least 1.")
  if (!model %in% 0:3) stop("'model' must be 0, 1, 2, or 3.")
  if (trimm <= 0 || trimm >= 0.5) stop("'trimm' must be between 0 and 0.5.")
  if (maxlags < 0) stop("'maxlags' must be non-negative.")
  if (simcv < 0) stop("'simcv' must be non-negative.")

  lagmethod <- match.arg(lagmethod)
  engine <- match.arg(engine)
  lagopt <- switch(lagmethod, tsig = 1L, fixed = 9L, zero = 0L, aic = 2L, bic = 3L)
  gssr <- engine == "gauss"

  tb <- max(1L, round(trimm * n))
  maxfeas <- .mkc_maxfeasible(n, trimm)
  if (m > maxfeas)
    stop(sprintf("m = %d exceeds the feasible maximum (%d) for n = %d and trimm = %g. Reduce m or trimm, or use a longer series.",
                 m, maxfeas, n, trimm))

  fit <- .mkc_test(datap, m, model, tb, lagopt, maxlags, engine, gssr)
  stat <- fit$tau
  bps <- fit$bp

  if (simcv > 0) {
    cv <- .mkc_simcv(k, m, model, trimm, lagopt, maxlags, engine, gssr,
                     simt, simcv, simseed)
    cv_source <- "simulated"
  } else if (m <= 5) {
    cv <- .mkc_cvtable(k, m, model)
    cv_source <- "table"
  } else {
    cv <- c(NA_real_, NA_real_, NA_real_)
    cv_source <- "none"
    warning("No tabulated critical values for m > 5; set 'simcv' to simulate them.")
  }

  reject_1 <- if (is.na(cv[1])) NA else stat < cv[1]
  reject_5 <- if (is.na(cv[2])) NA else stat < cv[2]
  reject_10 <- if (is.na(cv[3])) NA else stat < cv[3]

  if (is.na(cv[1])) {
    conclusion <- "Critical values unavailable; cannot perform inference (set 'simcv')."
  } else if (isTRUE(reject_1)) {
    conclusion <- "Reject H0 of no cointegration at the 1% level (cointegration with break(s))."
  } else if (isTRUE(reject_5)) {
    conclusion <- "Reject H0 of no cointegration at the 5% level (cointegration with break(s))."
  } else if (isTRUE(reject_10)) {
    conclusion <- "Reject H0 of no cointegration at the 10% level (cointegration with break(s))."
  } else {
    conclusion <- "Fail to reject H0: no evidence of cointegration."
  }

  vnames <- colnames(datap)
  if (is.null(vnames)) vnames <- c("y", paste0("x", seq_len(k)))

  out <- list(
    statistic = stat,
    breakpoints = bps,
    break_fractions = if (length(bps)) bps / n else numeric(0),
    critical_values = cv,
    cv_source = cv_source,
    sim_reps = if (cv_source == "simulated") simcv else NA_integer_,
    reject_1 = reject_1, reject_5 = reject_5, reject_10 = reject_10,
    conclusion = conclusion,
    m = m, model = model, n = n, k = k, trimm = trimm,
    lagmethod = lagmethod, engine = engine, lags = fit$lag,
    depvar = vnames[1], regnames = vnames[-1],
    data = datap, call = match.call()
  )
  class(out) <- "maki_test"
  out
}

#' @rdname coint_maki
#' @param x A \code{"maki_test"} object.
#' @param ... Unused.
#' @return \code{print.maki_test} returns \code{x} invisibly.
#' @export
print.maki_test <- function(x, ...) {
  model_names <- c("Level Shift", "Level Shift with Trend",
                   "Regime Shift", "Regime Shift with Trend")
  eng <- if (x$engine == "gauss") "GAUSS/tspdlib-compatible" else "Maki (2012) paper"
  cat("\nMaki (2012) Cointegration Test with Multiple Structural Breaks\n")
  cat(strrep("-", 62), "\n", sep = "")
  cat(sprintf("Model        : %s (model = %d)\n", model_names[x$model + 1L], x$model))
  cat(sprintf("Engine       : %s\n", eng))
  cat(sprintf("Observations : %d    Regressors: %d    Trimming: %.2f\n",
              x$n, x$k, x$trimm))
  cat(sprintf("Lag rule     : %s (max %d, used %d)\n", x$lagmethod,
              12L, x$lags))
  cat(sprintf("Dependent    : %s    Regressor(s): %s\n",
              x$depvar, paste(x$regnames, collapse = ", ")))
  cat(strrep("-", 62), "\n", sep = "")
  cat(sprintf("Test statistic : %10.4f\n", x$statistic))
  if (is.na(x$critical_values[1])) {
    cat("Critical values: not available (m > 5; set 'simcv' to simulate)\n")
  } else {
    cat(sprintf("Critical values: 1%% = %.3f   5%% = %.3f   10%% = %.3f  [%s]\n",
                x$critical_values[1], x$critical_values[2],
                x$critical_values[3], x$cv_source))
  }
  if (length(x$breakpoints)) {
    cat("\nEstimated breaks (obs : fraction):\n")
    for (i in seq_along(x$breakpoints))
      cat(sprintf("  break %d : %d : %.4f\n", i, x$breakpoints[i],
                  x$break_fractions[i]))
  }
  cat(strrep("-", 62), "\n", sep = "")
  cat(x$conclusion, "\n\n")
  invisible(x)
}

#' Plot a Maki Cointegration Test Result
#'
#' Draws a two-panel dashboard with 'ggplot2': the dependent variable with its
#' break-adjusted long-run fit and shaded regimes (top), and the cointegrating
#' residual about zero (bottom), with dashed lines at the estimated breaks.
#'
#' @param x A \code{"maki_test"} object from \code{\link{coint_maki}}.
#' @param ... Unused.
#' @return A \code{ggplot} object (invisibly), or \code{NULL} if 'ggplot2' is not
#'   installed.
#' @export
#' @examples
#' set.seed(1)
#' n <- 60
#' x <- cumsum(rnorm(n))
#' y <- 2 + 3 * x + 4 * (seq_len(n) > 30) + rnorm(n)
#' res <- coint_maki(cbind(y, x), m = 1, model = 2)
#' if (requireNamespace("ggplot2", quietly = TRUE)) plot(res)
plot.maki_test <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Install the 'ggplot2' package to plot a maki_test object.")
    return(invisible(NULL))
  }
  datap <- x$data
  n <- x$n
  y <- datap[, 1]
  Xf <- .mkc_xfixed(datap, n, x$model, x$breakpoints)
  b <- solve(crossprod(Xf), crossprod(Xf, y))
  fitv <- as.vector(Xf %*% b)
  resid <- y - fitv
  obs <- seq_len(n)

  pa <- "(a) Series and break-adjusted long-run relation"
  pb <- "(b) Cointegrating residual (equilibrium error)"
  df <- rbind(
    data.frame(obs = obs, value = y,     series = "observed", panel = pa),
    data.frame(obs = obs, value = fitv,  series = "fitted",   panel = pa),
    data.frame(obs = obs, value = resid, series = "residual", panel = pb)
  )
  zline <- data.frame(panel = pb, yint = 0)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$obs, y = .data$value,
                                        colour = .data$series)) +
    ggplot2::geom_hline(data = zline, ggplot2::aes(yintercept = .data$yint),
                        colour = "grey50", linewidth = 0.3, inherit.aes = FALSE) +
    ggplot2::geom_line(linewidth = 0.7, na.rm = TRUE) +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
    ggplot2::scale_colour_manual(
      values = c(observed = "navy", fitted = "firebrick", residual = "forestgreen"),
      breaks = c("observed", "fitted"), name = NULL) +
    ggplot2::labs(
      title = sprintf("Maki (2012) cointegration test - %s",
                      c("Level Shift", "Level Shift with Trend",
                        "Regime Shift", "Regime Shift with Trend")[x$model + 1L]),
      subtitle = x$conclusion, x = "Observation", y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top",
                   plot.subtitle = ggplot2::element_text(size = 9))
  if (length(x$breakpoints)) {
    p <- p + ggplot2::geom_vline(xintercept = x$breakpoints,
                                 linetype = "dashed", colour = "grey45",
                                 linewidth = 0.4)
  }
  print(p)
  invisible(p)
}
