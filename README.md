# makicoint: Maki Cointegration Test with Multiple Structural Breaks

[![CRAN status](https://www.r-pkg.org/badges/version/makicoint)](https://CRAN.R-project.org/package=makicoint)

## Overview

`makicoint` implements the Maki (2012) residual-based test for cointegration
allowing for an **unknown number of structural breaks**. It extends the
Gregory-Hansen (one break) and Hatemi-J (two breaks) tests to any feasible
number of breaks, and is useful when the long-run relationship between
non-stationary series may shift more than twice or is subject to regime changes.

## Key features

- **Any feasible number of breaks** (not capped at five), bounded by the sample
  size and the trimming parameter.
- **Four model specifications**: level shift (`0`), level shift with trend (`1`),
  regime shift (`2`), and regime shift with trend (`3`).
- **Two engines** giving the same test statistic: the default reproduces the
  original GAUSS/tspdlib implementation; `engine = "paper"` uses the Maki (2012,
  Steps 2 and 4) break rule.
- **Critical values** from Maki (2012) Table 1 (depending on the number of
  regressors, breaks and model); **simulated** critical values for more than five
  breaks via `simcv`.
- **Diagnostic plot** (`ggplot2`): the series with its break-adjusted long-run
  fit, and the cointegrating residual.
- ADF lag rules: `tsig` (default), `fixed`, `zero`, `aic`, `bic`.

## Installation

```r
install.packages("makicoint")                       # CRAN
# devtools::install_github("merwanroudane/makicoint")  # development
```

## Usage

```r
library(makicoint)

set.seed(123)
n <- 100
x <- cumsum(rnorm(n))
y <- 0.5 * x + cumsum(rnorm(n))
y[51:100] <- y[51:100] + 2          # a level break at observation 50

res <- coint_maki(cbind(y, x), m = 1, model = 0)
res
plot(res)                           # requires ggplot2
```

Two breaks, regime-shift model, and the paper engine:

```r
coint_maki(cbind(y, x), m = 2, model = 2)
coint_maki(cbind(y, x), m = 2, model = 2, engine = "paper")
```

Beyond five breaks, with simulated critical values (heavy):

```r
coint_maki(cbind(y, x), m = 7, simcv = 2000, simt = 500)
```

## Functions

- `coint_maki()` — the test.
- `cv_coint_maki(k, m, model)` — Maki (2012) Table 1 critical values.
- `print()` / `plot()` methods for the result.

## Reference

Maki, D. (2012). Tests for cointegration allowing for an unknown number of
breaks. *Economic Modelling*, 29, 2011-2015.
[doi:10.1016/j.econmod.2012.04.022](https://doi.org/10.1016/j.econmod.2012.04.022)

## Author

Dr Merwan Roudane — Independent Researcher — merwanroudane920@gmail.com —
[github.com/merwanroudane](https://github.com/merwanroudane)

## License

GPL-3
