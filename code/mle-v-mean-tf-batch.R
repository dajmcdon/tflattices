library(batchtools)
library(tflattices)
library(tidyverse)

# functions ---------------------------------------------------------------
make_signal_v <- function(n) {
  # V shaped (tf has trouble at the boundaries, so put the hard part 
  # in the middle)
  n <- as.integer(round(n))
  n2 <- n %/% 2
  r <- seq(1, n2 + 1, length.out = n2 + 1)
  l <- rev(r)
  r <- r[-1]
  if (n %% 2 == 0L) l <- l[-(n2 + 1)]
  return(c(l,r) / n)
}

# make_signal_triangle_wave <- function(n, period = 1 / floor(sqrt(n) / 2)) {
#   # similar to the v, but we need ~\sqrt(n) teeth instead of 1
#   tt <- 1:n / n
#   x <- 2 * abs(tt / p - floor(tt / p + 1 / 2))
#   10 * (-1 * (x - 1) + 2 / sqrt(n) ) / n  
# }


log10_seq <- function(from = 1, to = 1, length.out = 1, int = TRUE) {
  stopifnot(from > 0, to > 0, length.out > 0)
  s <- 10^seq(log10(from), log10(to), length.out = length.out)
  if (int) s <- round(s)
  s
}

MSE <- function(x, y) {
  # x a matrix, y a vector
  colMeans((x - y)^2)
}

KL_div <- function(x, y, distribution = c("exponential", "poisson")) {
  # x a matrix, y a vector
  distribution <- match.arg(distribution)
  switch(
    distribution,
    # for exponential, we have the mean = 1/lambda, so we invert everything
    exponential = colMeans(ifelse(x > 0, log(x / y) + y / x - 1, Inf)), 
    poisson = colMeans(y - x + ifelse(x > 0, x * log(x / y), x))
  )
}

# parameters --------------------------------------------------------------

nrepls <- 10

# data
problem_design <- list(
  vsignal = data.table::CJ( # fully factorial design
    n = log10_seq(20, 1000, 20),
    distribution = c("exponential", "poisson"),
    control_tv = c("mean", "nat_par")
  )
)

# algo
algo_design <- list( # both tf_versions, all other params the same
  tf = tibble(
    tf_version = c("mean", "mle"), 
    maxit = 25000L, 
    korder = 1L,
    tolerance = 1e-6,
    penalize_nullspace = FALSE,
    lambda_min = -1,
    lambda_max = -1
  )
)



# data --------------------------------------------------------------------

generator <- function(data, job, # required by batchtools
                      n, 
                      distribution = c("exponential", "poisson"),
                      control_tv = c("mean","nat_par"),
                      signal = c("v", "wave")) {
  distribution <- match.arg(distribution)
  control_tv <- match.arg(control_tv)
  signal <- match.arg(signal)
  sig <- switch(signal, # has tv controlled
                v = make_signal_v(n),
                wave = make_signal_triangle_wave(n))
  if (distribution == "exponential") {
    if (control_tv == "mean") sig <- 1 / sig # lambda = 1 / theta
    y <- rexp(n, rate = sig) # sig is always the rate
    sig <- 1 / sig # now sig is the mean, for MSE
  }
  if (distribution == "poisson") {
    sig <- (0.5 - sig) + log(n)# hard for mean if this is nat param
    if (control_tv == "nat_par") sig <- exp(sig) # theta = log(lambda)
    y <- rpois(n, lambda = sig) # sig is the mean already
  }
  return(list(y = y, sig = sig, distribution = distribution))
}


# estimation --------------------------------------------------------------

estimator <- function(data, job, # required by batch tools 
                      instance, # contains the output of generator()
                      tf_version = c("mean", "mle"),
                      ...) { # other args passed to etfgrid() in algo_design
  distribution <- instance$distribution
  y <- instance$y
  sig <- instance$sig
  tf_version <- match.arg(tf_version)
  if (tf_version == "mean") {
    thhat <- etfgrid(y, exp_fam = "gauss_mean", ...)
    thhat$theta[thhat$theta < 0] = 0 # in both cases, these must be positive
  } else {
    thhat <- etfgrid(y, exp_fam = distribution, ...) # returns the mean parameter
    thhat$theta[thhat$theta < 0] = 0 # in both cases, these must be positive, shouldn't be an issue
  }
  loc_mse <- which.min(MSE(thhat$theta, sig))
  loc_kl <- which.min(KL_div(thhat$theta, sig, distribution))
  lambda_mse <- thhat$lambda[loc_mse]
  lambda_kl<- thhat$lambda[loc_kl]
  return(list(
    mse = min(MSE(thhat$theta, sig)),
    kl = min(KL_div(thhat$theta, sig, distribution)),
    lambda_mse = lambda_mse,
    lambda_kl = lambda_kl,
    thhat_mse = drop(thhat$theta[,loc_mse]),
    thhat_kl = drop(thhat$theta[,loc_kl]),
    est = thhat
  ))
}


# setup the experiment and submit -----------------------------------------


#makeExperimentRegistry("../experiments/mle-v-mean-tf", 
makeExperimentRegistry("mle-v-mean-tf", packages = "tflattices")
addProblem(name = "vsignal", data = NULL, fun = generator, seed = 12345)
addAlgorithm(name = "tf", fun = estimator)


addExperiments(problem_design, algo_design, repls = nrepls)
ids <- bind_cols(findJobs(), n = unwrap(getJobPars())$n)
ids$chunk <- lpt(ids$n, parallel::detectCores())
submitJobs(ids) 