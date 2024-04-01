#ifndef __LOSSFUNCTIONS_H
#define __LOSSFUNCTIONS_H


void EntrywiseSoftThreshold(arma::vec& z, double lam);

double gaussupdate(double& vy, const double& mu);

double poisupdate(double& vy, const double& mu);

double gaussvarupdate(double& vy, const double& mu);

double expupdate(double& vy, const double& mu);


#endif
