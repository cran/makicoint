# Break detection functions for 3, 4, and 5 breaks
# These functions follow the same pattern as mbreak1 and mbreak2

#' @keywords internal
mbreak31 <- function(datap, n, model, tb, bp1, bp2, lagoption) {
  y <- datap[, 1, drop = FALSE]
  k <- ncol(datap)
  
  du1 <- rbind(matrix(0, bp1, 1), matrix(1, n - bp1, 1))
  du2 <- rbind(matrix(0, bp2, 1), matrix(1, n - bp2, 1))
  if (k > 1) {
    dx1 <- rbind(matrix(0, bp1, k - 1), datap[(bp1 + 1):n, 2:k, drop = FALSE])
    dx2 <- rbind(matrix(0, bp2, k - 1), datap[(bp2 + 1):n, 2:k, drop = FALSE])
  }
  dtr1 <- rbind(matrix(0, bp1, 1), matrix((bp1 + 1):n, ncol = 1))
  dtr2 <- rbind(matrix(0, bp2, 1), matrix((bp2 + 1):n, ncol = 1))
  
  vectau <- numeric(n)
  vecbp <- numeric(n)
  
  start_i <- tb + 1
  end_i <- bp1 - tb
  
  if (start_i >= end_i) {
    return(list(tau = NA, bp = NA))
  }
  
  for (i in start_i:end_i) {
    u <- matrix(1, n, 1)
    du3 <- rbind(matrix(0, i, 1), matrix(1, n - i, 1))
    
    if (model == 0) {
      const <- cbind(u, du1, du2, du3)
      x <- cbind(const, datap[, 2:k, drop = FALSE])
    } else if (model == 1) {
      const <- cbind(u, du1, du2, du3)
      tr <- matrix(1:n, ncol = 1)
      x <- cbind(const, tr, datap[, 2:k, drop = FALSE])
    } else if (model == 2) {
      const <- cbind(u, du1, du2, du3)
      dx3 <- rbind(matrix(0, i, k - 1), datap[(i + 1):n, 2:k, drop = FALSE])
      x <- cbind(const, datap[, 2:k, drop = FALSE], dx1, dx2, dx3)
    } else if (model == 3) {
      tr <- matrix(1:n, ncol = 1)
      dtr3 <- rbind(matrix(0, i, 1), matrix((i + 1):n, ncol = 1))
      const <- cbind(u, du1, du2, du3, tr, dtr1, dtr2, dtr3)
      dx3 <- rbind(matrix(0, i, k - 1), datap[(i + 1):n, 2:k, drop = FALSE])
      x <- cbind(const, datap[, 2:k, drop = FALSE], dx1, dx2, dx3)
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
  mintau31 <- min(vectau[valid_range])
  bp31 <- valid_range[which.min(vecbp[valid_range])]
  
  return(list(tau = mintau31, bp = bp31))
}

#' @keywords internal
mbreak32 <- function(datap, n, model, tb, bp1, bp2, lagoption) {
  y <- datap[, 1, drop = FALSE]
  k <- ncol(datap)
  
  du1 <- rbind(matrix(0, bp1, 1), matrix(1, n - bp1, 1))
  du2 <- rbind(matrix(0, bp2, 1), matrix(1, n - bp2, 1))
  if (k > 1) {
    dx1 <- rbind(matrix(0, bp1, k - 1), datap[(bp1 + 1):n, 2:k, drop = FALSE])
    dx2 <- rbind(matrix(0, bp2, k - 1), datap[(bp2 + 1):n, 2:k, drop = FALSE])
  }
  dtr1 <- rbind(matrix(0, bp1, 1), matrix((bp1 + 1):n, ncol = 1))
  dtr2 <- rbind(matrix(0, bp2, 1), matrix((bp2 + 1):n, ncol = 1))
  
  vectau <- numeric(n)
  vecbp <- numeric(n)
  
  start_i <- bp1 + tb
  end_i <- bp2 - tb
  
  if (start_i >= end_i) {
    return(list(tau = NA, bp = NA))
  }
  
  for (i in start_i:end_i) {
    u <- matrix(1, n, 1)
    du3 <- rbind(matrix(0, i, 1), matrix(1, n - i, 1))
    
    if (model == 0) {
      const <- cbind(u, du1, du2, du3)
      x <- cbind(const, datap[, 2:k, drop = FALSE])
    } else if (model == 1) {
      const <- cbind(u, du1, du2, du3)
      tr <- matrix(1:n, ncol = 1)
      x <- cbind(const, tr, datap[, 2:k, drop = FALSE])
    } else if (model == 2) {
      const <- cbind(u, du1, du2, du3)
      dx3 <- rbind(matrix(0, i, k - 1), datap[(i + 1):n, 2:k, drop = FALSE])
      x <- cbind(const, datap[, 2:k, drop = FALSE], dx1, dx2, dx3)
    } else if (model == 3) {
      tr <- matrix(1:n, ncol = 1)
      dtr3 <- rbind(matrix(0, i, 1), matrix((i + 1):n, ncol = 1))
      const <- cbind(u, du1, du2, du3, tr, dtr1, dtr2, dtr3)
      dx3 <- rbind(matrix(0, i, k - 1), datap[(i + 1):n, 2:k, drop = FALSE])
      x <- cbind(const, datap[, 2:k, drop = FALSE], dx1, dx2, dx3)
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
  mintau32 <- min(vectau[valid_range])
  bp32 <- valid_range[which.min(vecbp[valid_range])]
  
  return(list(tau = mintau32, bp = bp32))
}

#' @keywords internal
mbreak33 <- function(datap, n, model, tb, bp1, bp2, lagoption) {
  y <- datap[, 1, drop = FALSE]
  k <- ncol(datap)
  
  du1 <- rbind(matrix(0, bp1, 1), matrix(1, n - bp1, 1))
  du2 <- rbind(matrix(0, bp2, 1), matrix(1, n - bp2, 1))
  if (k > 1) {
    dx1 <- rbind(matrix(0, bp1, k - 1), datap[(bp1 + 1):n, 2:k, drop = FALSE])
    dx2 <- rbind(matrix(0, bp2, k - 1), datap[(bp2 + 1):n, 2:k, drop = FALSE])
  }
  dtr1 <- rbind(matrix(0, bp1, 1), matrix((bp1 + 1):n, ncol = 1))
  dtr2 <- rbind(matrix(0, bp2, 1), matrix((bp2 + 1):n, ncol = 1))
  
  vectau <- numeric(n)
  vecbp <- numeric(n)
  
  start_i <- bp2 + tb
  end_i <- n - tb
  
  if (start_i >= end_i) {
    return(list(tau = NA, bp = NA))
  }
  
  for (i in start_i:end_i) {
    u <- matrix(1, n, 1)
    du3 <- rbind(matrix(0, i, 1), matrix(1, n - i, 1))
    
    if (model == 0) {
      const <- cbind(u, du1, du2, du3)
      x <- cbind(const, datap[, 2:k, drop = FALSE])
    } else if (model == 1) {
      const <- cbind(u, du1, du2, du3)
      tr <- matrix(1:n, ncol = 1)
      x <- cbind(const, tr, datap[, 2:k, drop = FALSE])
    } else if (model == 2) {
      const <- cbind(u, du1, du2, du3)
      dx3 <- rbind(matrix(0, i, k - 1), datap[(i + 1):n, 2:k, drop = FALSE])
      x <- cbind(const, datap[, 2:k, drop = FALSE], dx1, dx2, dx3)
    } else if (model == 3) {
      tr <- matrix(1:n, ncol = 1)
      dtr3 <- rbind(matrix(0, i, 1), matrix((i + 1):n, ncol = 1))
      const <- cbind(u, du1, du2, du3, tr, dtr1, dtr2, dtr3)
      dx3 <- rbind(matrix(0, i, k - 1), datap[(i + 1):n, 2:k, drop = FALSE])
      x <- cbind(const, datap[, 2:k, drop = FALSE], dx1, dx2, dx3)
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
  mintau33 <- min(vectau[valid_range])
  bp33 <- valid_range[which.min(vecbp[valid_range])]
  
  return(list(tau = mintau33, bp = bp33))
}

#' @keywords internal
mbreak3 <- function(datap, n, model, tb, lagoption) {
  result2 <- mbreak2(datap, n, model, tb, lagoption)
  bp1 <- result2$bp[1]
  bp2 <- result2$bp[2]
  
  result31 <- mbreak31(datap, n, model, tb, bp1, bp2, lagoption)
  result32 <- mbreak32(datap, n, model, tb, bp1, bp2, lagoption)
  result33 <- mbreak33(datap, n, model, tb, bp1, bp2, lagoption)
  
  taus <- c(result31$tau, result32$tau, result33$tau)
  
  valid_idx <- which(!is.na(taus))
  if (length(valid_idx) == 0) {
    return(list(tau = result2$tau, bp = result2$bp))
  }
  
  min_idx <- valid_idx[which.min(taus[valid_idx])]
  mintau3 <- taus[min_idx]
  
  if (min_idx == 1) {
    breaks <- sort(c(result31$bp, bp1, bp2))
  } else if (min_idx == 2) {
    breaks <- sort(c(bp1, result32$bp, bp2))
  } else {
    breaks <- sort(c(bp1, bp2, result33$bp))
  }
  
  return(list(tau = mintau3, bp = breaks))
}
