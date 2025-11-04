#' Compute DF Tau Statistic
#'
#' Computes the Dickey-Fuller tau statistic for unit root testing.
#'
#' @param dy Numeric vector of first differences of residuals.
#' @param x Matrix of lagged residuals and differenced lags.
#' @return A list containing:
#'   \item{tau}{The tau statistic (t-statistic for the first coefficient)}
#'   \item{s2}{Sum of squared residuals}
#' @keywords internal
dftau <- function(dy, x) {
  n1 <- nrow(x)
  k <- ncol(x)
  
  # OLS estimation
  b <- solve(t(x) %*% x) %*% t(x) %*% dy
  e <- dy - x %*% b
  s2 <- t(e) %*% e
  s2hat <- as.numeric(s2) / (n1 - k)
  
  # Standard errors
  varb <- s2hat * solve(t(x) %*% x)
  se <- sqrt(diag(varb))
  
  # t-statistic for first coefficient
  tau <- b[1] / se[1]
  
  return(list(tau = as.numeric(tau), s2 = as.numeric(s2)))
}

#' Determine Optimal Lag Using t-sig Criterion
#'
#' Determines the optimal lag length using the t-significance criterion
#' with a maximum lag of 12.
#'
#' @param e Numeric vector of residuals.
#' @return Integer indicating the optimal lag length.
#' @keywords internal
opttlag <- function(e) {
  maxlag <- 12
  p <- maxlag
  
  while (p >= 1) {
    n <- length(e)
    dy <- diff(e)
    j <- 2 + p
    
    if (j > n) {
      p <- p - 1
      next
    }
    
    # Create lagged matrix
    x <- matrix(e[(j-1):(n-1)], ncol = 1)
    
    for (i in 1:p) {
      if ((j-1-i) < 1 || (n-1-i) < 1) {
        p <- p - 1
        next
      }
      lag_dy <- dy[(j-1-i):(n-1-i)]
      x <- cbind(x, lag_dy)
    }
    
    n1 <- nrow(x)
    k <- ncol(x)
    ly <- dy[(j-1):(n-1)]
    
    if (n1 != length(ly)) {
      p <- p - 1
      next
    }
    
    # OLS
    b <- solve(t(x) %*% x) %*% t(x) %*% ly
    e2 <- ly - x %*% b
    s2hat <- as.numeric(t(e2) %*% e2) / (n1 - k)
    varb <- s2hat * solve(t(x) %*% x)
    se <- sqrt(diag(varb))
    
    taut <- b[p + 1] / se[p + 1]
    
    if (abs(taut) > 1.654) {
      break
    }
    p <- p - 1
  }
  
  lag <- max(p, 0)
  return(lag)
}
