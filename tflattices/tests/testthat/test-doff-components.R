library(Matrix)

test_that("phipp works", {
  spdiag <- function(x) {
    n <- length(x)
    sparseMatrix(seq(n), seq(n), x = x)
  }
  # normal
  expect_equal(phipp(1L, rnorm(5)), spdiag(rep(1, 5)))
  # Gauss Variance 1 / (2x^2)
  expect_equal(phipp(2L, as.double(1:4)), spdiag(1 / (2*(1:4)^2)))
  # Poisson
  expect_equal(phipp(3L, as.double(1:4)), spdiag(exp(1:4)))
  # Exponential
  expect_equal(phipp(4L, as.double(1:6)), spdiag(1 / (1:6)^2))
})
