# makicoint 2.0.0

Major correctness and feature release. Results from version 1.0.0 should be
re-run, as the critical values and the break search have been corrected.

## Corrections

* **Critical values now match Maki (2012) Table 1 exactly** and depend on the
  number of regressors (1-4), as the table requires. `cv_coint_maki()` takes the
  number of regressors `k` (was the unused sample size `n`).
* The break search now runs for **all** of `m = 1, ..., 5` (and beyond);
  version 1.0.0 stopped at three breaks despite advertising five.
* The reference DOI in the DESCRIPTION is corrected to
  <doi:10.1016/j.econmod.2012.04.022>.

## New features

* **Any feasible number of breaks** via a general recursive search (no longer
  capped), bounded only by the sample size and the trimming parameter.
* **Two engines** with the same test statistic: the default reproduces the
  original 'GAUSS'/'tspdlib' implementation; `engine = "paper"` uses the Maki
  (2012, Steps 2 and 4) break rule.
* **Simulated critical values** for more than five breaks via `simcv`, using
  Maki's Monte-Carlo design.
* **`plot()` method** ('ggplot2'): a two-panel dashboard of the series with its
  break-adjusted long-run fit and the cointegrating residual.
* Additional ADF lag rules: `lagmethod = "tsig"` (default), `"fixed"`, `"zero"`,
  `"aic"`, `"bic"`.

## Interface changes

* `coint_maki()` gains `lagmethod`, `maxlags`, `engine`, `simcv`, `simt` and
  `simseed`; the default `trimm` is now `0.10`. The legacy `lagoption` argument
  is replaced by `lagmethod`.
* The result object adds `break_fractions`, `cv_source`, `lags`, `engine` and
  `k`, and stores the data to support `plot()`.

# makicoint 1.0.0

* First CRAN release: Maki (2012) cointegration test with structural breaks.
