#ifndef __DOFFESTIMATE_H
#define __DOFFESTIMATE_H


arma::mat projnull(arma::sp_mat const& D, arma::uvec S);

arma::sp_mat phipp(int lossfunction, arma::vec th);

double dof(int lossfunction, double lambda, arma::vec theta, 
           arma::sp_mat const& D, arma::mat const& B, 
           arma::vec const& y, 
           int pen_null, double alpha, double tolerance);
#endif