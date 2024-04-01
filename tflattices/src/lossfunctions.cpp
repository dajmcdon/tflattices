#include <RcppArmadillo.h>
#include <boost/math/special_functions/lambert_w.hpp>
#include "lossfunctions.h"


void EntrywiseSoftThreshold(arma::vec& z, double lam) {
  z.transform( [&](double x) {
    return ((x > 0) - (x < 0)) * std::max(0.0, std::fabs(x) - lam);
  } );
}

double gaussupdate(double& vy, const double& mu) { // DJM: correct
  // this is prox_fmu(v) where f=sum f_i and f_i= phi(theta_i) - theta_i y_i
  // therefore, we need to solve phi'(x) + mu*x = v + y
  // phi(x) = x^2/2 (Gaussian)
  return vy / (1 + mu);
}

double poisupdate(double& vy, const double& mu) { // DJM: correct
  // this is prox_fmu(v) where f=sum f_i and f_i= phi(theta_i) - theta_i y_i
  // therefore, we need to solve phi'(x) + mu*x = v + y
  // phi(x) = exp(x) (Poisson)
  vy /= mu;
  if (vy < 500) { // deal with potential overflow
    vy -= boost::math::lambert_w0(exp(vy) / mu);
  } else {
    // See: https://en.wikipedia.org/wiki/Lambert_W_function#Asymptotic_expansions
    // We use the first four terms.
    double la, lb;
    la = vy - log(mu);
    lb = log(la);
    vy -= la - lb + lb/la + (lb*(lb-2))/(2*la*la);
  }
  return vy;

}

double gaussvarupdate(double& vy, const double& mu) {
  // this is prox_fmu(v) where f=sum f_i and f_i= phi(theta_i) - theta_i y_i
  // therefore, we need to solve phi'(x) + mu*x = v + y
  // phi(x) = -log(-x)/2 (N(0, sig^2))
  // always negative
  vy += sqrt(2*mu + vy*vy);
  vy /= (2*mu);
  return vy;
}

double expupdate(double& vy, const double& mu) { //DJM: correct
  // this is prox_fmu(v) where f=sum f_i and f_i= phi(theta_i) - theta_i y_i
  // therefore, we need to solve phi'(x) + mu*x = v + y
  // \phi(x) = -\ln(-x), result should be negative always
  vy += sqrt(4*mu + vy*vy);
  vy /= 2*mu;
  return vy;
}

// double bernupdate(double& vy, const double& mu){
  // this is prox_fmu(v) where f=sum f_i and f_i= phi(theta_i) - theta_i y_i
  // therefore, we need to solve phi'(x) + mu*x = v + y
  // \phi(x) = \ln(1+exp(x))


//}
//
//
// This needs a value for alpha. Add to the linear updater function in admm when
// complete.

// double gamupdate(double& vy, const double& mu){
//   // this is prox_fmu(v) where f=sum f_i and f_i= phi(theta_i) - theta_i y_i
//   // therefore, we need to solve phi'(x) + mu*x = v + y
//   // \phi(x) = -\ln(-x)
//   vy = (-sqrt(-4*alpha/mu + (vy * vy)) + vy) / (2/mu);
//   return vy;
// }

