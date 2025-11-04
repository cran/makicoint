test_that("coint_maki works with model 0 and m=1", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + cumsum(rnorm(n))
  y[51:100] <- y[51:100] + 2
  
  data <- cbind(y, x)
  result <- coint_maki(data, m = 1, model = 0)
  
  expect_s3_class(result, "maki_test")
  expect_true(is.numeric(result$statistic))
  expect_true(!is.null(result$breakpoints))
  expect_true(result$m == 1)
  expect_true(result$model == 0)
})

test_that("coint_maki works with model 1 and m=2", {
  set.seed(456)
  n <- 150
  x <- cumsum(rnorm(n))
  y <- 0.6 * x + cumsum(rnorm(n))
  y[51:100] <- y[51:100] + 1.5
  y[101:150] <- y[101:150] + 3
  
  data <- cbind(y, x)
  result <- coint_maki(data, m = 2, model = 1)
  
  expect_s3_class(result, "maki_test")
  expect_true(is.numeric(result$statistic))
  expect_true(!is.null(result$breakpoints))
  expect_true(result$m == 2)
  expect_true(result$model == 1)
})

test_that("coint_maki validates input correctly", {
  expect_error(coint_maki(matrix(1:10, ncol=1), m = 1, model = 0))
  expect_error(coint_maki(cbind(1:10, 1:10), m = -1, model = 0))
  expect_error(coint_maki(cbind(1:10, 1:10), m = 1, model = 5))
  expect_error(coint_maki(cbind(1:10, 1:10), m = 6, model = 0))
})

test_that("cv_coint_maki returns correct critical values", {
  cv <- cv_coint_maki(100, m = 1, model = 0)
  
  expect_true(is.numeric(cv))
  expect_true(length(cv) == 3)
  expect_true(all(cv < 0))
  # Critical values should be ordered: 1% < 5% < 10% (more negative to less negative)
  expect_true(cv[1] < cv[2])  # 1% more negative than 5%
  expect_true(cv[2] < cv[3])  # 5% more negative than 10%
})

test_that("print method works", {
  set.seed(789)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + cumsum(rnorm(n))
  
  data <- cbind(y, x)
  result <- coint_maki(data, m = 1, model = 0)
  
  expect_output(print(result), "Maki Cointegration Test")
  expect_output(print(result), "Test Statistic")
  expect_output(print(result), "Critical Values")
})

test_that("Different model specifications work", {
  set.seed(321)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + cumsum(rnorm(n))
  y[51:100] <- y[51:100] + 2
  
  data <- cbind(y, x)
  
  # Model 0: Level shift
  result0 <- coint_maki(data, m = 1, model = 0)
  expect_s3_class(result0, "maki_test")
  expect_equal(result0$model, 0)
  
  # Model 1: Level shift with trend
  result1 <- coint_maki(data, m = 1, model = 1)
  expect_s3_class(result1, "maki_test")
  expect_equal(result1$model, 1)
  
  # Model 2: Regime shift
  result2 <- coint_maki(data, m = 1, model = 2)
  expect_s3_class(result2, "maki_test")
  expect_equal(result2$model, 2)
  
  # Model 3: Trend and regime shift
  result3 <- coint_maki(data, m = 1, model = 3)
  expect_s3_class(result3, "maki_test")
  expect_equal(result3$model, 3)
})

test_that("Test works with multiple explanatory variables", {
  set.seed(555)
  n <- 100
  x1 <- cumsum(rnorm(n))
  x2 <- cumsum(rnorm(n))
  y <- 0.5 * x1 + 0.3 * x2 + cumsum(rnorm(n))
  y[51:100] <- y[51:100] + 2
  
  data <- cbind(y, x1, x2)
  result <- coint_maki(data, m = 1, model = 0)
  
  expect_s3_class(result, "maki_test")
  expect_true(is.numeric(result$statistic))
})

test_that("m=0 (no breaks) works correctly", {
  set.seed(999)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + cumsum(rnorm(n))
  
  data <- cbind(y, x)
  result <- coint_maki(data, m = 0, model = 0)
  
  expect_s3_class(result, "maki_test")
  expect_true(is.numeric(result$statistic))
  expect_true(is.null(result$breakpoints))  # No breaks with m=0
  expect_equal(result$m, 0)
})

test_that("All critical value tables work", {
  # Test all combinations of m and model that should work
  for (m_val in 0:3) {
    for (model_val in 0:3) {
      cv <- cv_coint_maki(100, m = m_val, model = model_val)
      expect_true(is.numeric(cv))
      expect_equal(length(cv), 3)
      expect_true(all(cv < 0))
      expect_true(cv[1] < cv[2])  # 1% < 5%
      expect_true(cv[2] < cv[3])  # 5% < 10%
    }
  }
})
