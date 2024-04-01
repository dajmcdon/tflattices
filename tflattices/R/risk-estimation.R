#' Unbiased estimation of the KL loss
#'
#' @param object an object of class `etfgrid` as returned by a call to [etfgrid()].
#'   The object must include an estimate of the degrees of freedom, except in 
#'   the case of Poisson. For all other cases one must use 
#'   `etfgrid(..., estimate_df = TRUE)`.
#'
#' @return A data frame with two or three columns. The first contains the values 
#'   of `lambda` used by [etfgrid()] the second contains unbiased estimates of 
#'   the expected Kullback-Leibler loss of the fit. Models with lower values
#'   are preferred. When the family is other than Poisson, the estimated
#'   degrees is the third column in the output. The returned object is an S3 
#'   object with plot and print methods
#' @export
#'
#' @examples
#' y <- rnorm(100, c(rep(0,30), rep(5,40), rep(0,30)))
#' out <- etfgrid(y, estimate_df = TRUE, lambda_max = 10)
#' risk <- kl_estimate(out)
kl_estimate <- function(object) {
  stopifnot(class(object)[1] == "etfgrid")
  fam <- object$call$exp_fam[1]
  assertthat::assert_that(
    object$call$estimate_df,
    msg = paste("An estimate of the degrees of freedom is required.",
                "Call `etfgrid` again with `estimate_df = TRUE`."))
  if (fam == "binomial") stop("binomial KL estimation is not possible.")
  y <- drop(object$y)
  n <- length(y)
  est <- switch(
    fam,
    gauss_mean = mse_estimate(object$theta, y) + 2 * object$df / n,
    gauss_var = kl_gauss_var(object$theta, y) + object$df / n,
    poisson = kl_poisson(object, y),
    exponential = kl_exponential(object$theta, y) + object$df / n,
    rep(NA, length(object$lambda)) # returns NA if none of the above families
  )
  out <- data.frame(lambda = object$lambda, kl_estimate = est, df = object$df)
  class(out) <- "etfgrid_risk"
  out
}

mse_estimate <- function(th, y) {
  colMeans((y - th)^2)
}

kl_gauss_var <- function(th, y) {
  h <- 1
  dh <- 0
  phi <- -log(th/2)
  dphi <- 1 / (2*th)
  colMeans((th + dh/h) * dphi) - colMeans(phi) 
}

kl_exponential <- function(th, y) {
  h <- 1
  dh <- 0
  phi <- -log(-th)
  dphi <- 1 / th
  colMeans((th + dh/h) * dphi) - colMeans(phi) 
}

kl_poisson <- function(object, y) {
  colMeans(abs(object$theta)) - object$df / length(y)
}


