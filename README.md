# makicoint: Maki Cointegration Test with Structural Breaks

[![CRAN status](https://www.r-pkg.org/badges/version/makicoint)](https://CRAN.R-project.org/package=makicoint)

## Overview

The `makicoint` package implements the Maki (2012) cointegration test that allows for an unknown number of structural breaks. This test is particularly useful for detecting long-run equilibrium relationships between non-stationary time series when structural changes occur over time.

## Key Features

- **Multiple Break Detection**: Tests for cointegration with up to 5 structural breaks
- **Four Model Specifications**:
  - Model 0: Level shifts
  - Model 1: Level shifts with trend
  - Model 2: Regime shifts (changes in both intercept and slope)
  - Model 3: Trend and regime shifts
- **Automatic Lag Selection**: Uses the t-sig criterion for optimal lag determination
- **Critical Values**: Based on Maki (2012) Monte Carlo simulations

## Installation

### From CRAN (once published)

```r
install.packages("makicoint")
```

### Development Version from GitHub

```r
# install.packages("devtools")
devtools::install_github("merwanroudane/makicoint")
```

## Usage

### Basic Example

```r
library(makicoint)

# Generate cointegrated series with a structural break
set.seed(123)
n <- 100
e1 <- rnorm(n)
e2 <- rnorm(n)

# Create I(1) variables
x <- cumsum(e1)
y <- 0.5 * x + cumsum(e2)

# Add a structural break at observation 50
y[51:100] <- y[51:100] + 2

# Combine into matrix (dependent variable first)
data <- cbind(y, x)

# Run Maki test with 1 break, level shift model
result <- coint_maki(data, m = 1, model = 0)
print(result)
```

### Testing for Two Breaks

```r
# Generate data with two breaks
set.seed(456)
n <- 150
x2 <- cumsum(rnorm(n))
y2 <- 0.6 * x2 + cumsum(rnorm(n))

# Add two structural breaks
y2[51:100] <- y2[51:100] + 1.5
y2[101:150] <- y2[101:150] + 3

data2 <- cbind(y2, x2)

# Test with m=2
result2 <- coint_maki(data2, m = 2, model = 0)
print(result2)
```

### Model Specifications

- **model = 0**: Level shift only
  - Tests for breaks in the intercept
  
- **model = 1**: Level shift with trend
  - Includes a deterministic time trend
  
- **model = 2**: Regime shift
  - Tests for breaks in both intercept and slope coefficients
  
- **model = 3**: Trend and regime shift
  - Most general model with trend and coefficient changes

## Mathematical Framework

The test is based on the following cointegration regression:

**Model 0 (Level Shift):**
```
y_t = μ + Σ(μ_i * D_i,t) + β'x_t + u_t
```

**Model 3 (Trend and Regime Shift):**
```
y_t = μ + Σ(μ_i * D_i,t) + βt + Σ(β_i * DT_i,t) + γ'x_t + Σ(γ_i' * D_i,t * x_t) + u_t
```

Where D_i,t are dummy variables indicating structural breaks.

## Interpretation

The test statistic is the minimum of the ADF tau statistics computed across all possible break point combinations. 

- **Reject the null hypothesis**: Evidence of cointegration with structural breaks
- **Fail to reject**: No evidence of cointegration

Lower (more negative) test statistics provide stronger evidence against the null hypothesis of no cointegration.

## Functions

- `coint_maki()`: Main function to perform the Maki cointegration test
- `cv_coint_maki()`: Get critical values for the test
- `print.maki_test()`: Print method for test results

## Reference

Maki, D. (2012). Tests for cointegration allowing for an unknown number of breaks. *Economic Modelling*, 29(5), 2011-2015. [https://doi.org/10.1016/j.econmod.2012.05.006](https://doi.org/10.1016/j.econmod.2012.05.006)

## Author

Dr. Merwan Roudane  
Independent Researcher  
Email: merwanroudane920@gmail.com

## License

GPL-3

## Contributing

Contributions, bug reports, and feature requests are welcome! Please feel free to open an issue or submit a pull request on GitHub.
