
#include <RcppArmadillo.h>
#include <cmath>
#include "builddmat.h"

using namespace Rcpp;

double chooseC(int n, int k) {
  return Rf_choose(n, k);
}

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

arma::sp_mat leftIdkron(arma::sp_mat const &A, int ids){
  const arma::uword A_rows = A.n_rows;
  const arma::uword A_cols = A.n_cols;

  arma::sp_mat B(A_rows*ids, A_cols*ids);

  for(arma::uword j = 0; j < ids; j++){
    B.submat(j*A_rows, j*A_cols, (j+1)*A_rows-1, (j+1)*A_cols-1) = A;
  }
  return B;
}

arma::sp_mat rightIdkron(arma::sp_mat const &A, int ids){
  const arma::uword A_rows = A.n_rows;
  const arma::uword A_cols = A.n_cols;

  arma::sp_mat B(A_rows*ids, A_cols*ids);
  arma::sp_mat Z = arma::speye(ids,ids);

  for(arma::sp_mat::const_iterator j = A.begin(); j != A.end(); ++j){
    B.submat(j.row()*ids, j.col()*ids, (j.row()+1)*ids-1,
             (j.col()+1)*ids-1) = (*j)*Z;
  }
  return B;
}

// [[Rcpp::export]]
arma::sp_mat builddmat(arma::ivec lens, arma::ivec ords, arma::ivec wrap){
  // This function is only used if d>1
  int d = lens.n_elem;
  int ncols = prod(lens);
  arma::ivec ms = (1-wrap) % (lens-ords-1) + wrap % lens;
  arma::ivec somerows = ncols/lens % ms;
  arma::ivec rowidx = arma::cumsum(somerows);
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

arma::mat expandgrid(arma::ivec ords){
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

arma::mat nullD(arma::ivec lens, arma::ivec ords, arma::ivec wrap){

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
