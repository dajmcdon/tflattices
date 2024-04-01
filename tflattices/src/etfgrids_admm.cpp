#include <RcppArmadillo.h>
#include "lossfunctions.h"
#include "builddmat.h"
#include "admm.h"
#include "doffestimate.h"

// [[Rcpp::depends(BH)]]

using namespace Rcpp;



// [[Rcpp::export]]
List etfgrid_c(
    arma::vec Y, arma::ivec dimensions, arma::ivec korder, arma::ivec wrap,
    int lossfunction,
    int nsol = 50,
    double lambdamin = -1, double lambdamax = -1,
    arma::vec lambda = NumericVector::create(),
    double mu_adjust = -1, double rho_adjust = -1,
    double mu = -1, double rho = -1,
    int pen_null = 0, double alpha = .1, int dfestim = 0,
    int maxiter = 100, double tolerance = 1e-3, int verbose = 0) {


  // Sanity checks
  if (nsol < 1) stop("Expected nsol > 0");
  if (maxiter < 1) stop("Expected maxiter > 0");
  if (tolerance <= 0.0) stop("Expected tolerance > 0");
  if (pen_null==1 & alpha <= 0) stop("Null space penalty requires alpha > 0.");

  int n = Y.n_elem;
  if (lambda.size() > 0) nsol = lambda.size();

  // Placeholders for solutions
  arma::mat theta(n, nsol);
  arma::vec niter(nsol);
  arma::vec df(nsol);
  df.zeros();

  // Build D matrix
  arma::sp_mat D;
  if (dimensions.n_elem ==1) {
    int wrapdim = wrap(0)*dimensions(0)+(1-wrap(0))*(dimensions(0)-korder(0)-1);
    D = buildD(wrapdim, dimensions(0), korder(0), wrap(0));
  } else {
    D = builddmat(dimensions, korder, wrap);
  }
  arma::sp_mat DD = D.t() * D;
  double Dmax_eig = 0; // max eigenvalue of DtD
  for (int i=0; i<dimensions.n_elem; i++) Dmax_eig += std::pow(4.0, korder(i)+1);

  // Generate lambda sequence if necessary
  arma::vec _lambda;
  if (lambda.size() > 0) {
    _lambda = arma::vec(lambda.begin(), lambda.size(), false);
    lambdamin = _lambda.min();
    lambdamax = _lambda.max();
    nsol = lambda.size();
  } else {
    // Where do we start lambda
    // We want it to start at ||(D*D')^(-1) * D * (Y-\phi'(\hat\theta))||_infty
    // We could try to get an upper bound, but how to avoid the matrix inverse?
    if (lambdamax <= 0) {
      lambdamax = arma::norm(D * Y, "inf"); // (D * D')^-1 * D * Y
      lambdamax *= n;
    }
    lambdamin = (lambdamin < 0) ? .001*lambdamax : lambdamin;
    // log10-spaced sequence
    _lambda = arma::logspace(log10(lambdamin), log10(lambdamax), nsol);
  }

  // ADMM parameters
  arma::mat B;
  double tolerance_abs = std::sqrt(n) * tolerance; // adjust for number of parameters
  double _rho = (rho < 0) ? lambdamin : rho;
  double _mu = (mu < 0) ? _rho*(Dmax_eig + 1e-2) : mu;

  
  if (pen_null == 1) {
    if (mu < 0) _mu += lambdamax*alpha;
    B = nullD(dimensions, korder, wrap);
  }

  // ADMM variables
  arma::vec z(D.n_rows, arma::fill::zeros);
  arma::vec u(D.n_rows, arma::fill::zeros);
  arma::vec x(n, arma::fill::zeros);
  // arma::vec x = B * B.t() * Y; this would start in the null space
  x += Y;
  initilizer(x, lossfunction);
  
  // Outer loop to compute solution path
  for (int i = 0; i < nsol; i++) {

    if (i > 0) {
      if (rho < 0) _rho = _lambda(i);
      if (mu < 0) _mu = _rho*(Dmax_eig + 1e-2);
    }

    if (verbose > 0) Rcout << ".";

    // accumulator for iters, reset every iteration of admm
    int iters = 0;

    admm(lossfunction, DD, D,
         x, Y, z, u,
         _lambda(i), _mu, _rho,
         pen_null, B, alpha,
         mu_adjust, rho_adjust,
         maxiter, tolerance_abs, iters);

    // Store solution
    theta.col(i) = x;
    niter(i) = iters;
    
    if (dfestim == 1) {
      df(i) = dof(lossfunction, _lambda(i), x, D, 
         B, Y, pen_null, alpha, tolerance);
    }

    // Verbose handlers
    if (verbose > 1) Rcout << niter(i);
    if (verbose > 2) Rcout << "(" << _lambda(i) << ")";
    if (verbose > 0) Rcout << std::endl;
  }

  

  // Return
  List out = List::create(
    Named("y") = Y,
    Named("n") = n,
    Named("lambda") = _lambda,
    Named("theta") = theta,
    Named("df") = df, 
    Named("niter") = niter
  );
  out.attr("class") = "etfgrids";
  return out;
}
