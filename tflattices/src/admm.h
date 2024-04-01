#ifndef __ADMM_H
#define __ADMM_H

void linear_updater(arma::vec& x, const double& mu, int linupdate);

void initilizer(arma::vec& x, int linupdate);

void admm(int linupdate, const arma::sp_mat& DD,
          const arma::sp_mat& D,
          arma::vec& x, const arma::vec& y, arma::vec& z, arma::vec& u,
          double& lam, double& mu, double& rho,
          int pen_null, const arma::mat& B, double& alpha,
          const double& mu_adjust, const double& rho_adjust,
          int maxiter, const double& tolerance, int& iters);

#endif
