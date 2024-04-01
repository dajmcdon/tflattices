#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
using namespace Rcpp;


// This script is for testing

double chooseC(int n, int k) {
  return Rf_choose(n, k);
}



// [[Rcpp::export]]
arma::sp_mat buildD(int m, int n, int ord, int wrap){
  arma::sp_mat D(m,n);
  if(wrap==1) D.set_size(n,n);
  int c1 = ord + 1;
  double z;
  for(int i=0; i <= c1; i++){
    z = chooseC(c1, i) * std::pow(-1, i+ord+1);
    D.diag(i) += z;
  }
  if(wrap==1){
    for(int i=ord; i>-1; i--){
      for(int j=0; j<=i; j++){
        D(n-j-1,i-j) += z;
      }
      z = chooseC(c1, i) * std::pow(-1, i+ord+1);
    }
  }
  return D;
}

// [[Rcpp::export]]
arma::sp_mat leftIdkron(arma::sp_mat const &A, int ids){
  const arma::uword A_rows = A.n_rows;
  const arma::uword A_cols = A.n_cols;

  arma::sp_mat B(A_rows*ids, A_cols*ids);

  for(arma::uword j = 0; j < ids; j++){
    B.submat(j*A_rows, j*A_cols, (j+1)*A_rows-1, (j+1)*A_cols-1) = A;
  }
  return B;
}

// [[Rcpp::export]]
arma::sp_mat rightIdkron(arma::sp_mat const &A, int ids){
  const arma::uword A_rows = A.n_rows;
  const arma::uword A_cols = A.n_cols;

  arma::sp_mat B(A_rows*ids, A_cols*ids);
  arma::sp_mat Z = arma::speye(ids,ids);

  for(arma::sp_mat::const_iterator j = A.begin(); j != A.end(); ++j){
      B.submat(j.row()*ids, j.col()*ids, (j.row()+1)*ids-1, (j.col()+1)*ids-1) = (*j)*Z;
  }
  return B;
}

// [[Rcpp::export]]
arma::sp_mat builddmat(arma::uvec lens, arma::uvec ords, arma::uvec wrap){
  // This function is only used if d>1
  int d = lens.n_elem;
  int ncols = prod(lens);
  arma::uvec ms = (1-wrap) % (lens-ords-1) + wrap % lens;
  arma::uvec somerows = ncols/lens % ms;
  arma::uvec rowidx = arma::cumsum(somerows);
  int nrows = rowidx(d-1);

  arma::sp_mat BigD(nrows, ncols);
  for(int i=0; i<d; i++){
    arma::sp_mat D = buildD(ms(i), lens(i), ords(i), wrap(i));
    if (i == 0){
      BigD.rows(0,rowidx(i)-1) = rightIdkron(D, ncols/lens(i));
    } else if (i == d-1){
      BigD.rows(rowidx(i-1), rowidx(i)-1) = leftIdkron(D, ncols/lens(i));
    } else {
      int lhs = 1; int rhs = 1;
      for(int j = i+1; j < d; j++) rhs *= lens(j);
      for(int j = 0; j < i; j++) lhs *= lens(j);
      arma::sp_mat temp = rightIdkron(D, rhs);
      BigD.rows(rowidx(i-1),rowidx(i)-1) = leftIdkron(temp, lhs);
    }
  }
  return BigD;
}


// This script is for testing
//

// [[Rcpp::export]]
arma::vec repelem(arma::vec x, int times) {
  arma::vec y = arma::repelem(x, times,1);
  return y;
}
// [[Rcpp::export]]
arma::vec repmat(arma::vec x, int times){
  arma::vec y = arma::repmat(x,times,1);
  return y;
}

// [[Rcpp::export]]
arma::mat expandgrid(arma::uvec ords){

  int d = ords.n_elem;
  ords += 1;
  int repfac = 1;
  int orep = prod(ords);
  arma::mat cargs(orep, d);
  for (int i=0; i<d; i++) {
    arma::vec x = arma::regspace(1,ords(i));
    int nx = x.n_elem;
    orep /= nx;
    arma::vec temp = arma::repelem(x, repfac, 1);
    cargs.col(i) = arma::repmat(temp, orep, 1);
    repfac *= nx;
  }
  cargs -= 1;
  return cargs;
}

//[[Rcpp::export]]
arma::mat nullD(arma::uvec lens, arma::uvec ords, arma::uvec wrap){

  int d = lens.n_elem;
  int nrows = prod(lens);
  int nullity = 1;
  // regardless of k, wrapped dimensions have nullity 1, others (k+1)
  for(int j=0; j<d; j++) nullity *= (wrap(j) == 1) ? 1 : ords(j)+1;
  arma::mat nsp(nrows, nullity);
  nsp.ones();
  arma::vec temp(nrows);

  arma::mat grid = expandgrid(ords);

  for(int j=0; j<d; j++){
    if(ords(j) < 1) continue;
    arma::vec x = arma::regspace(1, lens(j));
    x /= lens(j);
    if(j==0){
      temp = arma::repelem(x, nrows / lens(j), 1);
    } else if (j==d-1){
      temp = arma::repmat(x, nrows / lens(j), 1);
    } else {
      int lhs = 1; int rhs = 1;
      for(int z = j+1; z < d; z++) rhs *= lens(z);
      for(int z = 0; z < j; z++) lhs *= lens(z);
      arma::vec ltemp = arma::repelem(x, rhs, 1);
      temp = arma::repmat(ltemp, lhs, 1);
    }
    // Now insert temp somewhere
    int times = ords(j);
    while(times > 0){
      times--;
      for(int z = 0; z < nullity; z++){
        if(grid(z,j) < 1) continue;
        nsp.col(z) %= temp;
      }
      grid.col(j) -= 1;
    }
  }
  arma::mat U;
  arma::vec s;
  arma::mat V;

  svd_econ(U, s, V, nsp, "right");
  nsp *= V * arma::diagmat(1/s);
  return nsp;
}





/*** R
library(Matrix)
nullD(10,0,0)
x=1:3
repelem(x,4) # same as rep.int(x,rep.int(4,3))
repmat(x,4) # same as rep.int(x,4)
rep.int(x,4)
rep.int(x,rep.int(4,3))
expand.grid(0:2,0:1,0,0:3)
expandgrid(c(2L,1L,0L,3L))
*/

// etf_grid testing

/***R
set.seed(123)
n = 16
gauss = c(rnorm(n*0.75, 0, 1), rnorm(n*0.25, 10, 1))

res = etfgrid(Y = gauss,
              dimensions = n,
              korder = 0,
              wrap = 0,
              lossfunction = 1,
              nsol = 5)

cbind(res$theta, gauss)
*/

/***R
pois = c(rpois(n*0.75, 5), rpois(n*0.25,10))

res = etfgrid(Y = pois,
              dimensions = n,
              korder = 0,
              wrap = 0,
              lossfunction = 2,
              nsol = 5)

cbind(res$theta, pois)
*/

/***R
set.seed(123)
n = 12
gaussvar = c(rnorm(n*0.75, 0, 1), rnorm(n*0.25, 0, sqrt(4)))
#gaussvar = rnorm(n, 0, 10)

res = etfgrid(Y = gaussvar,
              dimensions = n,
              korder = 0,
              wrap = 0, 
              lossfunction = 3,
              nsol = 5)

cbind(res$theta, gaussvar^2)
*/