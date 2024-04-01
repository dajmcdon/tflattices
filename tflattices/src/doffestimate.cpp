#include <RcppArmadillo.h>
#include "doffestimate.h"

// [[Rcpp::export]]
arma::mat projnull(arma::sp_mat const& D, arma::uvec S) {
  arma::mat A(D);
  A.shed_rows(S);
  int n = A.n_cols;
  arma::mat u;
  arma::vec ss;
  arma::mat v;
  arma::svd_econ(u, ss, v, A, "right");
  arma::mat Id;
  Id.eye(n, n);
  Id -= v * v.t();
  return Id;
}

// [[Rcpp::export]]
arma::sp_mat phipp(int lossfunction, arma::vec th) {
  int n = th.n_elem;
  arma::sp_mat phi(n,n);
  
  switch(lossfunction) {
    case 1 : th.transform( [&](double x) { return 1.0; } );
      break;  
    case 2 : th.transform( [&](double x) { return 1 / (2*x*x); } );
      break;
    case 3 : th.transform( [&](double x) { return exp(x);  } );
      break;
    case 4 : th.transform( [&](double x) { return 1 / (x*x);  } );
      break;
  }
  
  phi.diag() = th;
  
  return phi;
}

double dof(int lossfunction, double lambda, arma::vec theta, 
           arma::sp_mat const& D, arma::mat const& B, 
           arma::vec const& y, 
           int pen_null, double alpha, double tolerance) {
  
  double df = 0;
  arma::vec bl = D * theta;
  bl.transform([&](double x) { return ((x > 0) - (x < 0)) * x; } ); 
  arma::uvec S = arma::find(bl > tolerance);
  // z.transform( [&](double x) { return ((x > 0) - (x < 0)) * x; } ); 
  // arma::uvec S = arma::find(z > tolerance);
  
  if (lossfunction == 1 && pen_null == 0) {
    df += S.n_elem;
  } else {
    arma::mat Pn = projnull(D, S);
    arma::sp_mat ddphi = phipp(lossfunction, theta);
    arma::mat K = Pn * ddphi * Pn;
    if (pen_null == 1) {
      K += lambda * alpha * B * B.t();
    }
    arma::mat P = arma::pinv(K);
    P = Pn * P * Pn;
    arma::vec ddphiv(ddphi.diag());
    if (lossfunction == 3) {
      theta -= P.diag();
      arma::uvec good = arma::find(y);
      df += arma::as_scalar( y(good).t() * theta(good) );
    } else {
      df += arma::as_scalar( ddphiv.t() * P.diag() );
    }
  }
  
  return(df);
}
