# makicoint 1.0.0

## Initial Release

* First CRAN submission
* Implements Maki (2012) cointegration test with structural breaks
* Supports up to 3 structural breaks (m=0,1,2,3)
* Four model specifications (level shifts, trends, regime shifts)
* Automatic lag selection using t-sig criterion
* Complete documentation and examples
* Vignette with practical applications

## Features

* `coint_maki()`: Main test function
* `cv_coint_maki()`: Critical values lookup
* `print.maki_test()`: Formatted output
* Comprehensive help documentation
* Unit tests for reliability

## Reference

Maki, D. (2012). Tests for cointegration allowing for an unknown number 
of breaks. Economic Modelling, 29(5), 2011-2015.
