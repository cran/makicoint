#' Main Procedure for Maki Cointegration Test
#'
#' Internal procedure function that orchestrates the break detection process.
#'
#' @param y Matrix or data frame with dependent variable in first column.
#' @param m Maximum number of breaks (0-5).
#' @param model Model specification (0-3).
#' @param trimm Trimming parameter (default 0.15).
#' @param lagoption Lag selection option (0=no lags, 1=optimal lags).
#' @return List with test statistic and breakpoints.
#' @keywords internal
pf <- function(y, m, model, trimm = 0.15, lagoption = 1) {
  n <- nrow(y)
  tb <- floor(trimm * n)
  
  if (m == 0) {
    y_mat <- as.matrix(y[, 1, drop = FALSE])
    k <- ncol(y)
    u <- matrix(1, n, 1)
    
    if (model == 0) {
      x <- cbind(u, y[, 2:k, drop = FALSE])
    } else if (model == 1) {
      tr <- matrix(1:n, ncol = 1)
      x <- cbind(u, tr, y[, 2:k, drop = FALSE])
    } else if (model == 2) {
      x <- cbind(u, y[, 2:k, drop = FALSE])
    } else if (model == 3) {
      tr <- matrix(1:n, ncol = 1)
      x <- cbind(u, tr, y[, 2:k, drop = FALSE])
    }
    
    b <- solve(t(x) %*% x) %*% t(x) %*% y_mat
    e <- as.vector(y_mat - x %*% b)
    
    if (lagoption == 0) {
      lag <- 0
    } else {
      lag <- opttlag(e)
    }
    
    dy <- diff(e)
    r <- 2 + lag
    x_adf <- matrix(e[(r - 1):(n - 1)], ncol = 1)
    
    if (lag > 0) {
      for (q in 1:lag) {
        lag_dy <- dy[(r - 1 - q):(n - 1 - q)]
        x_adf <- cbind(x_adf, lag_dy)
      }
    }
    
    result <- dftau(dy[(r - 1):(n - 1)], x_adf)
    return(list(tau = result$tau, bp = NULL))
    
  } else if (m == 1) {
    result <- mbreak1(y, n, model, tb, lagoption)
    return(result)
    
  } else if (m == 2) {
    result <- mbreak2(y, n, model, tb, lagoption)
    return(result)
    
  } else if (m == 3) {
    result <- mbreak3(y, n, model, tb, lagoption)
    return(result)
  }
  
  # For m=4 and m=5, return simplified version
  # Full implementation would require mbreak4 and mbreak5
  stop("m=4 and m=5 require additional functions. Use m <= 3.")
}

#' Get Critical Values for Maki Cointegration Test
#'
#' Returns critical values for the Maki cointegration test based on 
#' Maki (2012) Table 1.
#'
#' @param n Sample size (currently not used as values are asymptotic).
#' @param m Number of breaks (0-5).
#' @param model Model specification (0-3):
#'   \itemize{
#'     \item 0: Level shift
#'     \item 1: Level shift with trend
#'     \item 2: Regime shift
#'     \item 3: Trend and regime shift
#'   }
#' @return Numeric vector of critical values at 1\%, 5\%, and 10\% significance levels.
#' @export
#' @examples
#' # Get critical values for m=1, model=0
#' cv_coint_maki(100, m=1, model=0)
cv_coint_maki <- function(n, m, model) {
  if (m == 0) {
    if (model == 0) return(c(-3.96, -3.37, -3.07))
    if (model == 1) return(c(-4.32, -3.77, -3.47))
    if (model == 2) return(c(-4.65, -4.03, -3.75))
    if (model == 3) return(c(-4.99, -4.42, -4.13))
  }
  
  if (m == 1) {
    if (model == 0) return(c(-5.11, -4.50, -4.21))
    if (model == 1) return(c(-5.42, -4.82, -4.53))
    if (model == 2) return(c(-5.74, -5.17, -4.88))
    if (model == 3) return(c(-6.16, -5.59, -5.28))
  }
  
  if (m == 2) {
    if (model == 0) return(c(-5.71, -5.10, -4.79))
    if (model == 1) return(c(-5.97, -5.41, -5.09))
    if (model == 2) return(c(-6.41, -5.83, -5.50))
    if (model == 3) return(c(-6.82, -6.22, -5.91))
  }
  
  if (m == 3) {
    if (model == 0) return(c(-6.11, -5.49, -5.18))
    if (model == 1) return(c(-6.41, -5.83, -5.51))
    if (model == 2) return(c(-6.92, -6.32, -6.01))
    if (model == 3) return(c(-7.31, -6.71, -6.40))
  }
  
  if (m == 4) {
    if (model == 0) return(c(-6.42, -5.80, -5.48))
    if (model == 1) return(c(-6.74, -6.16, -5.84))
    if (model == 2) return(c(-7.33, -6.73, -6.41))
    if (model == 3) return(c(-7.72, -7.13, -6.81))
  }
  
  if (m == 5) {
    if (model == 0) return(c(-6.69, -6.07, -5.75))
    if (model == 1) return(c(-7.02, -6.44, -6.12))
    if (model == 2) return(c(-7.69, -7.09, -6.77))
    if (model == 3) return(c(-8.07, -7.48, -7.16))
  }
  
  return(c(NA, NA, NA))
}

#' Maki Cointegration Test with Structural Breaks
#'
#' Performs the Maki (2012) cointegration test that allows for an unknown 
#' number of structural breaks. The test detects cointegration relationships
#' in the presence of up to five structural breaks.
#'
#' @param y Matrix or data frame with dependent variable in first column and
#'   independent variable(s) in remaining columns.
#' @param m Maximum number of breaks to test (0-3). Values 4-5 require extended implementation.
#' @param model Model specification (0-3):
#'   \itemize{
#'     \item 0: Level shift
#'     \item 1: Level shift with trend  
#'     \item 2: Regime shift (changes in intercept and slope)
#'     \item 3: Trend and regime shift
#'   }
#' @param trimm Trimming parameter (default 0.15). Determines the minimum 
#'   distance between breaks as a fraction of sample size.
#' @param lagoption Lag selection (0=no lags, 1=optimal lags using t-sig criterion).
#'
#' @return A list with class "maki_test" containing:
#'   \item{statistic}{The test statistic (minimum tau)}
#'   \item{breakpoints}{Vector of detected break point locations}
#'   \item{critical_values}{Critical values at 1\%, 5\%, and 10\% levels}
#'   \item{reject_1}{Logical; reject null at 1\% level?}
#'   \item{reject_5}{Logical; reject null at 5\% level?}
#'   \item{reject_10}{Logical; reject null at 10\% level?}
#'   \item{conclusion}{Text conclusion of the test}
#'   \item{m}{Number of breaks tested}
#'   \item{model}{Model specification used}
#'   \item{n}{Sample size}
#'
#' @references
#' Maki, D. (2012). Tests for cointegration allowing for an unknown number
#' of breaks. \emph{Economic Modelling}, 29(5), 2011-2015.
#'
#' @export
#' @examples
#' # Generate cointegrated series with one break
#' set.seed(123)
#' n <- 100
#' e1 <- rnorm(n)
#' e2 <- rnorm(n)
#' x <- cumsum(e1)
#' y <- 0.5 * x + cumsum(e2)
#' y[51:100] <- y[51:100] + 2  # Add structural break
#' 
#' # Run Maki test
#' data <- cbind(y, x)
#' result <- coint_maki(data, m=1, model=0)
#' print(result)
coint_maki <- function(y, m = 1, model = 0, trimm = 0.15, lagoption = 1) {
  # Input validation
  if (!is.matrix(y) && !is.data.frame(y)) {
    stop("y must be a matrix or data frame")
  }
  
  y <- as.matrix(y)
  n <- nrow(y)
  
  if (ncol(y) < 2) {
    stop("y must have at least 2 columns (dependent and independent variable)")
  }
  
  if (!m %in% 0:5) {
    stop("m must be between 0 and 5")
  }
  
  if (!model %in% 0:3) {
    stop("model must be 0, 1, 2, or 3")
  }
  
  if (trimm <= 0 || trimm >= 0.5) {
    stop("trimm must be between 0 and 0.5")
  }
  
  # Run test
  result <- pf(y, m, model, trimm, lagoption)
  cv <- cv_coint_maki(n, m, model)
  
  reject_1 <- result$tau < cv[1]
  reject_5 <- result$tau < cv[2]
  reject_10 <- result$tau < cv[3]
  
  if (reject_5) {
    conclusion <- "Reject null hypothesis: Evidence of cointegration with structural break(s) at 5% significance level."
  } else {
    conclusion <- "Fail to reject null hypothesis: No evidence of cointegration at 5% significance level."
  }
  
  output <- list(
    statistic = result$tau,
    breakpoints = result$bp,
    critical_values = cv,
    reject_1 = reject_1,
    reject_5 = reject_5,
    reject_10 = reject_10,
    conclusion = conclusion,
    m = m,
    model = model,
    n = n
  )
  
  class(output) <- "maki_test"
  return(output)
}

#' Print Method for Maki Test Results
#'
#' Prints a formatted summary of the Maki cointegration test results,
#' including test statistic, break points, critical values, and conclusions.
#'
#' @param x An object of class "maki_test" returned by \code{\link{coint_maki}}
#' @param ... Additional arguments (currently unused)
#'
#' @return Returns the input object \code{x} invisibly. This function is 
#'   called primarily for its side effect of printing formatted test results 
#'   to the console.
#'
#' @export
print.maki_test <- function(x, ...) {
  cat("\n")
  cat("============================================================\n")
  cat("Maki Cointegration Test with Structural Breaks\n")
  cat("============================================================\n\n")
  
  cat("Model Specification:\n")
  model_names <- c("Level Shift", "Level Shift with Trend", 
                   "Regime Shift", "Trend and Regime Shift")
  cat(sprintf("  Type: %s (model=%d)\n", model_names[x$model + 1], x$model))
  cat(sprintf("  Maximum breaks tested: %d\n", x$m))
  cat(sprintf("  Sample size: %d\n\n", x$n))
  
  cat("Test Results:\n")
  cat(sprintf("  Test Statistic: %.4f\n", x$statistic))
  cat(sprintf("  Critical Values: 1%%=%.2f, 5%%=%.2f, 10%%=%.2f\n\n", 
              x$critical_values[1], x$critical_values[2], x$critical_values[3]))
  
  if (!is.null(x$breakpoints) && length(x$breakpoints) > 0) {
    cat(sprintf("  Detected Break Points: %s\n", 
                paste(x$breakpoints, collapse=", ")))
    cat(sprintf("  Break Fractions: %s\n\n", 
                paste(round(x$breakpoints/x$n, 3), collapse=", ")))
  } else {
    cat("  Detected Break Points: None\n\n")
  }
  
  cat("Hypothesis Test:\n")
  cat(sprintf("  Reject at 1%% level: %s\n", ifelse(x$reject_1, "YES", "NO")))
  cat(sprintf("  Reject at 5%% level: %s\n", ifelse(x$reject_5, "YES", "NO")))
  cat(sprintf("  Reject at 10%% level: %s\n\n", ifelse(x$reject_10, "YES", "NO")))
  
  cat("Conclusion:\n")
  cat(sprintf("  %s\n", x$conclusion))
  cat("\n")
  cat("============================================================\n")
  
  invisible(x)
}
