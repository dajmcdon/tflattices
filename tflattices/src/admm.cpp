#include <RcppArmadillo.h>
#include "lossfunctions.h"
#include "admm.h"
// [[Rcpp::plugins("cpp11")]]


void linear_updater(arma::vec& x, const double& mu, int linupdate) {
  switch(linupdate) {
  case 1 : x.transform( [&](double x) { return(gaussupdate(x, mu)); } );
    break;
  case 2 : x.transform( [&](double x) { return(gaussvarupdate(x, mu)); } );
    break;
  case 3 : x.transform( [&](double x) { return(poisupdate(x, mu)); } );
    break;
  case 4 : x.transform( [&](double x) { return(expupdate(x, mu)); } );
    break;
  }
}

void initilizer(arma::vec& x, int linupdate) {
  switch(linupdate) {
  case 1 : break;
  case 2 : x.transform( [&](double z) { return(-1 / (2 * z)); });
    break;
  case 3 : break;
    // x.transform( [&](double z) { return(-1 * (z==0) + (z>0) * log(z)); });
    // break;
  case 4 : x.transform( [&](double z) { return(-1 / z); });
    break;
  }
}



// This used to return the number of iterations, but the easiest way to store
// x in the theta matrix was to change the return type and use a reference to
// store the number of iterations.
void admm(int linupdate, const arma::sp_mat& DD,
          const arma::sp_mat& D,
          arma::vec& x, const arma::vec& y, arma::vec& z, arma::vec& u,
          double& lam, double& mu, double& rho,
          int pen_null, const arma::mat& B, double& alpha,
          const double& mu_adjust, const double& rho_adjust,
          int maxiter, const double& tolerance, int& iters) {
  int niter;
  double rr, ss;
  arma::vec z_old(z.n_elem, arma::fill::zeros);
  arma::vec Dx(u.n_elem, arma::fill::zeros);
  

  for (niter = 0; niter < maxiter; niter++) {
    z_old = z; // Store previous value of z

    if (pen_null == 1) {
      x = mu*x - rho * (DD*x + D.t() * (u - z)) - lam * alpha * B*B.t() * x + y; // v+y
    } else {
      x = mu*x - rho * (DD*x + D.t() * (u - z)) + y;
    }

    // (see Parikh and Boyd 2014, p158 and Ouyang, He, and Gray 2013, p3)
    // PB2014 seems wrong. We need mu such that muI - rho DD is positive definite.
    // Thus, mu > rho * ||DD||_2. The below seems wrong but kept in case.
    //
    // _WRONG_ note: need 0< mu le lam/sig_max(DD)
    // want to send mu -> zero, start at lam/||D||^2_max < lam/sig_max(DD)
    // lambda there is 1/rho here, and mu there is 1/mu here
    // this avoids some ugliness with tiny numbers

    // our update is now to solve phi'(z) + mu*z = x for z
    linear_updater(x, mu, linupdate);
    Dx = D*x;
    z = Dx + u;
    EntrywiseSoftThreshold(z, lam/rho);

    u += Dx - z; // Dual variable update

    // Compute primal and dual residual norms
    rr = arma::norm(Dx - z, "fro"); //primal
    ss = rho * arma::norm(z - z_old, "fro"); //dual

    // Check convergence criterion
    if (rr < tolerance && ss < tolerance) {
      niter++;
      iters = niter;
      break;
    }

    // Penalty adjustment (Boyd, et al. 2011  eq 3.13)
    if (rho_adjust > 1) { // rho_adjust > 1, default is no adjustment, 2 is typical
      if (rr > 10.0 * ss) {
        rho *= rho_adjust;
        u /= rho_adjust;
        mu *= rho_adjust;
      } else if (ss > 10.0 * rr) {
        rho /= rho_adjust;
        u *= rho_adjust;
        mu /= rho_adjust;
      }
    }
    if (mu_adjust > 1) mu *= mu_adjust; // mu_adjust > 1, usually 5 or so.
  }
  iters = niter;
}
