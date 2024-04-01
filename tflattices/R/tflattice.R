

#' Exponential family trend filtering over grids
#'
#' @param y vector or array of observations
#' @param korder vector of trend filter orders, non-negative integers
#' @param exp_fam the name of the distribution/parameter, one of
#' @param dims vector of dimensions for the grid, can be specified if y is not
#'   an array but a vector of length (n1 x n2 x ... nd)
#' @param wrap logical vector indicating if any dimensions are wrapped (like
#'   a cylinder, say)
#' @param lambda non-negative, decreasing sequence of penalty parameters, will
#'   be determined automatically if omitted
#' @param n_solutions the number of solutions to use if lambda is omitted
#' @param penalize_nullspace should we include a penalty on the nullspace
#'   of the TF operator
#' @param nullspace_mult the nullspace penalty is lambda scaled by this value
#' @param return_mean_param logical, is the returned estimate the natural
#'   parameter (TRUE returns the mean parameter)
#' @param estimate_df logical, should we calculate the degrees of freedom. 
#'   This can potentially be slow. It is very 
#'   fast for `gaussian_mean` without a nullspace penalty and slower for other 
#'   options. Larger grids increase 
#'   the computing time very quickly.
#' @param ... additional arguments to configure the ADMM algorithm.
#'   See [admm_configuration()]
#'
#' @return an object of class [etfgrid]
#' @export
#'
etfgrid <- function(y,
                    korder = rep(0L, d),
                    exp_fam = c("gauss_mean","gauss_var","poisson","exponential","binomial"),
                    dims = dim(y),
                    wrap = rep(FALSE, d),
                    lambda = numeric(0),
                    n_solutions = 50L,
                    penalize_nullspace = FALSE,
                    nullspace_mult = 0.1,
                    return_mean_param = TRUE,
                    estimate_df = FALSE,
                    ...){
  if (is.null(dims)) dims <- length(y)
  d <- length(dims)
  the_call <- c(as.list(environment()), ...)
  pars <- admm_configuration(...)
  
  
  fam <-  match.arg(exp_fam)
  fam_cd <- match(fam, c("gauss_mean","gauss_var","poisson","exponential","binomial"))
  if (fam == "binomial") stop("binomial distribution not currently implemented.")
  stopifnot(length(korder) == d, length(wrap) == d, prod(dims) == length(y))
  stopifnot(all(korder >= 0), all(korder == floor(korder)))
  stopifnot(all(is.logical(wrap)))
  stopifnot(all(dims > korder + 2))
  if (length(lambda) > 0) {
    stopifnot(all(lambda >= 0))
    lambda <- sort(lambda, decreasing = FALSE)
  }
  n_solutions <- max(floor(n_solutions), length(lambda))
  stopifnot(n_solutions > 0)
  if (penalize_nullspace && nullspace_mult < 0) {
    penalize_nullspace <- FALSE
    warning("nullspace_mult > 0 if penalizing the null space")
  }
  if (fam == "gauss_var") y <- - (y - mean(y))^2
  if (fam == "poisson" && (any(round(y) != y) | any(y < 0)))
    stop("Poisson estimation requires non-negative integer y.")
  if (fam == "exponential") {
    if(any(y < 0)) stop("Exponential estimation requires non-negative y.")
    y <- -y
  }
  if (fam == "binomial" & !all(y %in% c(0, 1)))
    stop("Binomial estimation should have all 0/1 y.")
  
  y <- as.vector(y)

  out <- etfgrid_c(y, dims, korder, wrap, fam_cd, n_solutions,
                   pars$lambda_min, pars$lambda_max,
                   lambda, pars$mu_adjust, pars$rho_adjust, pars$mu,
                   pars$rho, as.integer(penalize_nullspace), nullspace_mult,
                   as.integer(estimate_df),
                   pars$maxiter, pars$tolerance, pars$verbose)
  if (fam == "gauss_var") out$y <- sqrt(-out$y)
  if (fam == "exponential") out$y <- -out$y
  class(out) <- "etfgrid"
  if (return_mean_param) {
    out$theta <- switch(fam,
                        gauss_var = 1 / (2 * out$theta),
                        poisson = exp(out$theta),
                        exponential = 1 / (out$theta),
                        out$theta)
  }
  out$call <- the_call
  out
}


#' Configuration parameters for [etfgrid()]
#' 
#' This function is not intended for direct usage. Passing these arguments
#' to [etfgrid()] controlls the algorithmic implementation.
#'
#' @param lambda_min Smallest value for the lambda sequence. Negative values
#'   (the default) result in automatic selection  (1e-4 * lambda_max)
#' @param lambda_max Largest value for the lambda sequence. Negative values
#'   (the default) result in automatic selection (infinity norm of `Dy` divided
#'   by `n`)
#' @param mu ADMM multiplier on the proximity part of the ADMM penalty (proximity to
#'   the previous iterate). Expert use only.
#' @param mu_adjust Adjustment multiplier for `mu` (enforces proximity between 
#'   the current and previous value). Expert use only.
#' @param rho ADMM multiplier on the primal/dual part of the ADMM penalty.
#' @param rho_adjust Adjustment multiplier for `rho`.
#' @param maxiter Maximum number of ADMM iterates.
#' @param tolerance Tolerance for stopping the iteration (before maxiter).
#' @param verbose Controls the printing of progress. Larger values result in
#'   more printing.
#'
#' @return a list. 
#' @export
admm_configuration <- function(lambda_min = -1, lambda_max = -1,
                               mu = -1, mu_adjust = -1,
                               rho = -1, rho_adjust = -1,
                               maxiter = 100L,
                               tolerance = 1e-3,
                               verbose = 0L){
  stopifnot(lambda_min == -1 || lambda_min >= 0)
  stopifnot(lambda_max == -1 || lambda_max > 0)
  if (lambda_min >= 0 && lambda_max > 0) stopifnot(lambda_max > lambda_min)
  stopifnot(mu == -1 || mu > 0, mu_adjust==-1 || mu_adjust > 0)
  stopifnot(rho == -1 || rho > 0, rho_adjust==-1 || rho_adjust > 0)
  maxiter = as.integer(maxiter)
  verbose = as.integer(verbose)
  stopifnot(maxiter > 0, tolerance > 0, verbose >= 0)
  list(lambda_min=lambda_min, lambda_max=lambda_max, mu=mu, mu_adjust=mu_adjust,
       rho=rho, rho_adjust=rho_adjust, maxiter=maxiter, tolerance=tolerance,
       verbose=verbose)
}
