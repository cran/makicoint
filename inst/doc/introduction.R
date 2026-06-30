## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 7, fig.height = 5)
has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

## ----load---------------------------------------------------------------------
library(makicoint)

## ----ex1----------------------------------------------------------------------
set.seed(123)
n <- 100
x <- cumsum(rnorm(n))
y <- 0.5 * x + cumsum(rnorm(n))
y[51:100] <- y[51:100] + 2          # a level break at observation 50

res <- coint_maki(cbind(y, x), m = 1, model = 0)
res

## ----ex1plot, eval = has_ggplot, fig.alt = "Two-panel Maki test dashboard"----
plot(res)

## ----engines------------------------------------------------------------------
g <- coint_maki(cbind(y, x), m = 2, model = 2)                 # default (GAUSS)
p <- coint_maki(cbind(y, x), m = 2, model = 2, engine = "paper")
c(gauss = g$statistic, paper = p$statistic)
rbind(gauss = g$breakpoints, paper = p$breakpoints)

## ----cv-----------------------------------------------------------------------
cv_coint_maki(k = 2, m = 3, model = 2)

## ----simcv, eval = FALSE------------------------------------------------------
# # (not run here; can take minutes)
# coint_maki(cbind(y, x), m = 7, simcv = 2000, simt = 500)

## ----results------------------------------------------------------------------
res$statistic
res$breakpoints
res$critical_values
res$cv_source

