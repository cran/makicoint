test_that("coint_maki returns a valid maki_test object", {
  set.seed(123)
  n <- 60
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + cumsum(rnorm(n))
  y[31:60] <- y[31:60] + 2
  res <- coint_maki(cbind(y, x), m = 1, model = 0)

  expect_s3_class(res, "maki_test")
  expect_true(is.numeric(res$statistic) && is.finite(res$statistic))
  expect_length(res$breakpoints, 1L)
  expect_identical(res$m, 1)
  expect_identical(res$model, 0)
  expect_identical(res$cv_source, "table")
})

test_that("number of breakpoints equals m", {
  set.seed(456)
  n <- 80
  x <- cumsum(rnorm(n))
  y <- 0.6 * x + cumsum(rnorm(n))
  res2 <- coint_maki(cbind(y, x), m = 2, model = 2)
  expect_length(res2$breakpoints, 2L)
  expect_true(all(diff(res2$breakpoints) > 0))
})

test_that("all four models run", {
  set.seed(321)
  n <- 60
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + cumsum(rnorm(n))
  y[31:60] <- y[31:60] + 2
  for (mod in 0:3) {
    res <- coint_maki(cbind(y, x), m = 1, model = mod)
    expect_s3_class(res, "maki_test")
    expect_identical(res$model, mod)
  }
})

test_that("paper engine runs and shares the statistic at m=1", {
  set.seed(11)
  n <- 60
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + cumsum(rnorm(n))
  g <- coint_maki(cbind(y, x), m = 1, model = 2, engine = "gauss")
  p <- coint_maki(cbind(y, x), m = 1, model = 2, engine = "paper")
  expect_equal(g$statistic, p$statistic, tolerance = 1e-8)
})

test_that("multiple regressors are supported", {
  set.seed(555)
  n <- 70
  x1 <- cumsum(rnorm(n)); x2 <- cumsum(rnorm(n))
  y <- 0.5 * x1 + 0.3 * x2 + cumsum(rnorm(n))
  res <- coint_maki(cbind(y, x1, x2), m = 1, model = 0)
  expect_s3_class(res, "maki_test")
  expect_identical(res$k, 2L)
})

test_that("input validation", {
  expect_error(coint_maki(matrix(1:10, ncol = 1)))          # one column
  expect_error(coint_maki(cbind(1:10, 1:10), model = 5))    # bad model
  expect_error(coint_maki(cbind(1:10, 1:10), m = 0))        # m < 1
  expect_error(coint_maki(cbind(rnorm(40), rnorm(40)), m = 30)) # infeasible
})

test_that("cv_coint_maki matches Maki (2012) Table 1", {
  expect_equal(cv_coint_maki(1, 1, 0), c(-5.709, -4.602, -4.354))
  expect_equal(cv_coint_maki(2, 3, 2), c(-7.031, -6.516, -6.210))
  expect_equal(cv_coint_maki(3, 3, 2), c(-7.767, -7.155, -6.868))
  expect_equal(cv_coint_maki(4, 5, 3), c(-10.08, -9.482, -9.151))
  expect_true(all(is.na(cv_coint_maki(5, 1, 0))))  # k out of range
  expect_true(all(is.na(cv_coint_maki(1, 6, 0))))  # m out of range
})

test_that("critical values are ordered 1% < 5% < 10%", {
  for (k in 1:4) for (m in 1:5) for (mod in 0:3) {
    cv <- cv_coint_maki(k, m, mod)
    expect_true(cv[1] < cv[2] && cv[2] < cv[3])
  }
})

test_that("print method works", {
  set.seed(789)
  n <- 60
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + cumsum(rnorm(n))
  res <- coint_maki(cbind(y, x), m = 1, model = 0)
  expect_output(print(res), "Maki")
  expect_output(print(res), "Test statistic")
})
