#' Detect One Structural Break
#'
#' Internal function to detect one structural break in the cointegration relationship.
#'
#' @param datap Data matrix with dependent variable in first column.
#' @param n Sample size.
#' @param model Model specification (0-3).
#' @param tb Trimming parameter (number of observations).
#' @param lagoption Lag selection option (0 or 1).
#' @return List with tau statistic and break point location.
#' @keywords internal
mbreak1 <- function(datap, n, model, tb, lagoption) {
  y <- datap[, 1, drop = FALSE]
  k <- ncol(datap)
  
  vectau <- numeric(n)
  vecbp <- numeric(n)
  
  for (i in (tb + 1):(n - tb)) {
    u <- matrix(1, n, 1)
    du <- rbind(matrix(0, i, 1), matrix(1, n - i, 1))
    
    if (model == 0) {  # Level shift
      const <- cbind(u, du)
      x <- cbind(const, datap[, 2:k, drop = FALSE])
      
    } else if (model == 1) {  # Level shift with trend
      const <- cbind(u, du)
      tr <- matrix(1:n, ncol = 1)
      x <- cbind(const, tr, datap[, 2:k, drop = FALSE])
      
    } else if (model == 2) {  # Regime shifts
      const <- cbind(u, du)
      dx <- rbind(matrix(0, i, k - 1), datap[(i + 1):n, 2:k, drop = FALSE])
      x <- cbind(const, datap[, 2:k, drop = FALSE], dx)
      
    } else if (model == 3) {  # Trend and regime shifts
      tr <- matrix(1:n, ncol = 1)
      dtr <- rbind(matrix(0, i, 1), matrix((i + 1):n, ncol = 1))
      const <- cbind(u, du, tr, dtr)
      dx <- rbind(matrix(0, i, k - 1), datap[(i + 1):n, 2:k, drop = FALSE])
      x <- cbind(const, datap[, 2:k, drop = FALSE], dx)
    }
    
    b <- solve(t(x) %*% x) %*% t(x) %*% y
    e <- as.vector(y - x %*% b)
    
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
    vectau[i] <- result$tau
    vecbp[i] <- result$s2
  }
  
  valid_range <- (tb + 1):(n - tb)
  mintau1 <- min(vectau[valid_range])
  bp1 <- valid_range[which.min(vecbp[valid_range])]
  
  return(list(tau = mintau1, bp = bp1))
}

#' Detect Second Break Before First Break
#' @keywords internal
mbreak21 <- function(datap, n, model, tb, bp1, lagoption) {
  y <- datap[, 1, drop = FALSE]
  k <- ncol(datap)
  
  du1 <- rbind(matrix(0, bp1, 1), matrix(1, n - bp1, 1))
  if (k > 1) {
    dx1 <- rbind(matrix(0, bp1, k - 1), datap[(bp1 + 1):n, 2:k, drop = FALSE])
  }
  dtr1 <- rbind(matrix(0, bp1, 1), matrix((bp1 + 1):n, ncol = 1))
  
  vectau <- numeric(n)
  vecbp <- numeric(n)
  
  start_i <- tb + 1
  end_i <- bp1 - tb
  
  if (start_i >= end_i) {
    return(list(tau = NA, bp = NA))
  }
  
  for (i in start_i:end_i) {
    u <- matrix(1, n, 1)
    du2 <- rbind(matrix(0, i, 1), matrix(1, n - i, 1))
    
    if (model == 0) {
      const <- cbind(u, du1, du2)
      x <- cbind(const, datap[, 2:k, drop = FALSE])
      
    } else if (model == 1) {
      const <- cbind(u, du1, du2)
      tr <- matrix(1:n, ncol = 1)
      x <- cbind(const, tr, datap[, 2:k, drop = FALSE])
      
    } else if (model == 2) {
      const <- cbind(u, du1, du2)
      dx2 <- rbind(matrix(0, i, k - 1), datap[(i + 1):n, 2:k, drop = FALSE])
      x <- cbind(const, datap[, 2:k, drop = FALSE], dx1, dx2)
      
    } else if (model == 3) {
      tr <- matrix(1:n, ncol = 1)
      dtr2 <- rbind(matrix(0, i, 1), matrix((i + 1):n, ncol = 1))
      const <- cbind(u, du1, du2, tr, dtr1, dtr2)
      dx2 <- rbind(matrix(0, i, k - 1), datap[(i + 1):n, 2:k, drop = FALSE])
      x <- cbind(const, datap[, 2:k, drop = FALSE], dx1, dx2)
    }
    
    b <- solve(t(x) %*% x) %*% t(x) %*% y
    e <- as.vector(y - x %*% b)
    
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
    vectau[i] <- result$tau
    vecbp[i] <- result$s2
  }
  
  valid_range <- start_i:end_i
  mintau21 <- min(vectau[valid_range])
  bp21 <- valid_range[which.min(vecbp[valid_range])]
  
  return(list(tau = mintau21, bp = bp21))
}

#' Detect Second Break After First Break
#' @keywords internal
mbreak22 <- function(datap, n, model, tb, bp1, lagoption) {
  y <- datap[, 1, drop = FALSE]
  k <- ncol(datap)
  
  du1 <- rbind(matrix(0, bp1, 1), matrix(1, n - bp1, 1))
  if (k > 1) {
    dx1 <- rbind(matrix(0, bp1, k - 1), datap[(bp1 + 1):n, 2:k, drop = FALSE])
  }
  dtr1 <- rbind(matrix(0, bp1, 1), matrix((bp1 + 1):n, ncol = 1))
  
  vectau <- numeric(n)
  vecbp <- numeric(n)
  
  start_i <- bp1 + tb
  end_i <- n - tb
  
  if (start_i >= end_i) {
    return(list(tau = NA, bp = NA))
  }
  
  for (i in start_i:end_i) {
    u <- matrix(1, n, 1)
    du2 <- rbind(matrix(0, i, 1), matrix(1, n - i, 1))
    
    if (model == 0) {
      const <- cbind(u, du1, du2)
      x <- cbind(const, datap[, 2:k, drop = FALSE])
      
    } else if (model == 1) {
      const <- cbind(u, du1, du2)
      tr <- matrix(1:n, ncol = 1)
      x <- cbind(const, tr, datap[, 2:k, drop = FALSE])
      
    } else if (model == 2) {
      const <- cbind(u, du1, du2)
      dx2 <- rbind(matrix(0, i, k - 1), datap[(i + 1):n, 2:k, drop = FALSE])
      x <- cbind(const, datap[, 2:k, drop = FALSE], dx1, dx2)
      
    } else if (model == 3) {
      tr <- matrix(1:n, ncol = 1)
      dtr2 <- rbind(matrix(0, i, 1), matrix((i + 1):n, ncol = 1))
      const <- cbind(u, du1, du2, tr, dtr1, dtr2)
      dx2 <- rbind(matrix(0, i, k - 1), datap[(i + 1):n, 2:k, drop = FALSE])
      x <- cbind(const, datap[, 2:k, drop = FALSE], dx1, dx2)
    }
    
    b <- solve(t(x) %*% x) %*% t(x) %*% y
    e <- as.vector(y - x %*% b)
    
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
    vectau[i] <- result$tau
    vecbp[i] <- result$s2
  }
  
  valid_range <- start_i:end_i
  mintau22 <- min(vectau[valid_range])
  bp22 <- valid_range[which.min(vecbp[valid_range])]
  
  return(list(tau = mintau22, bp = bp22))
}

#' Detect Two Structural Breaks
#' @keywords internal
mbreak2 <- function(datap, n, model, tb, lagoption) {
  result1 <- mbreak1(datap, n, model, tb, lagoption)
  bp1 <- result1$bp
  
  result21 <- mbreak21(datap, n, model, tb, bp1, lagoption)
  result22 <- mbreak22(datap, n, model, tb, bp1, lagoption)
  
  if (is.na(result21$tau) && is.na(result22$tau)) {
    return(list(tau = result1$tau, bp = c(bp1)))
  }
  
  if (is.na(result21$tau)) {
    mintau2 <- result22$tau
    breaks <- sort(c(bp1, result22$bp))
  } else if (is.na(result22$tau)) {
    mintau2 <- result21$tau
    breaks <- sort(c(result21$bp, bp1))
  } else {
    if (result21$tau <= result22$tau) {
      mintau2 <- result21$tau
      breaks <- sort(c(result21$bp, bp1))
    } else {
      mintau2 <- result22$tau
      breaks <- sort(c(bp1, result22$bp))
    }
  }
  
  return(list(tau = mintau2, bp = breaks))
}
