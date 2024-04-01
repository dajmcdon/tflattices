
#ifndef __BUILDDMAT_H
#define __BUILDDMAT_H

double chooseC(int n, int k);
arma::sp_mat buildD(int m, int n, int ord, int wrap);
arma::sp_mat leftIdkron(arma::sp_mat const &A, int ids);
arma::sp_mat rightIdkron(arma::sp_mat const &A, int ids);
arma::sp_mat builddmat(arma::ivec lens, arma::ivec ords, arma::ivec wrap);
arma::mat expandgrid(arma::ivec ords);
arma::mat nullD(arma::ivec lens, arma::ivec ords, arma::ivec wrap);


#endif
