## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----installation, eval=FALSE-------------------------------------------------
# # From CRAN
# install.packages("makicoint")
# 
# # Development version
# devtools::install_github("merwanroudane/makicoint")

## ----load---------------------------------------------------------------------
library(makicoint)

## ----example1-----------------------------------------------------------------
set.seed(123)
n <- 100
e1 <- rnorm(n)
e2 <- rnorm(n)

# Generate I(1) processes
x <- cumsum(e1)
y <- 0.5 * x + cumsum(e2)

# Add structural break at observation 50
y[51:100] <- y[51:100] + 2

# Plot the data
oldpar <- par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
plot(y, type = 'l', col = 'blue', main = 'Dependent Variable (Y) with Break',
     xlab = 'Time', ylab = 'Y')
abline(v = 50, col = 'red', lty = 2)
plot(x, type = 'l', col = 'darkgreen', main = 'Independent Variable (X)',
     xlab = 'Time', ylab = 'X')
par(oldpar)

## ----test1--------------------------------------------------------------------
# Prepare data matrix (Y first, then X)
data1 <- cbind(y, x)

# Run test with m=1 (one break), model=0 (level shift)
result1 <- coint_maki(data1, m = 1, model = 0)
print(result1)

## ----example2-----------------------------------------------------------------
set.seed(456)
n <- 150
e1 <- rnorm(n)
e2 <- rnorm(n)

x2 <- cumsum(e1)
y2 <- 0.6 * x2 + cumsum(e2)

# Add two breaks
y2[51:100] <- y2[51:100] + 1.5
y2[101:150] <- y2[101:150] + 3

# Plot
plot(y2, type = 'l', col = 'blue', main = 'Series with Two Breaks',
     xlab = 'Time', ylab = 'Y')
abline(v = c(50, 100), col = 'red', lty = 2)

# Test
data2 <- cbind(y2, x2)
result2 <- coint_maki(data2, m = 2, model = 0)
print(result2)

## ----example3-----------------------------------------------------------------
set.seed(789)
n <- 120
e1 <- rnorm(n)
e2 <- rnorm(n)

x3 <- cumsum(e1)
y3 <- 0.4 * x3 + cumsum(e2)

# Regime shift: change slope after observation 60
y3[61:120] <- y3[61:120] + 0.3 * x3[61:120]

# Plot
plot(y3, type = 'l', col = 'blue', main = 'Regime Shift',
     xlab = 'Time', ylab = 'Y')
abline(v = 60, col = 'red', lty = 2)

# Test with model=2
data3 <- cbind(y3, x3)
result3 <- coint_maki(data3, m = 1, model = 2)
print(result3)

## ----results------------------------------------------------------------------
# Access specific components
result1$statistic        # Test statistic
result1$breakpoints      # Break locations
result1$critical_values  # Critical values [1%, 5%, 10%]
result1$reject_5         # Reject at 5% level?

## ----session------------------------------------------------------------------
sessionInfo()

