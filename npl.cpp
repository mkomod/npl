#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec dnorm_dMu(arma::vec y, double mu, double sig)
{
    return pow(2.0 * PI, -0.5) * (y - mu) * pow(sig, -3.0) %
	exp( - 1.0/(2.0 * sig*sig) * pow(y - mu, 2.0));
}

// [[Rcpp::export]]
arma::vec dnorm_dSig(arma::vec y, double mu, double sig)
{
    return pow(2.0 * PI, -0.5) * pow(y - mu, 2.0) * pow(sig, -4.0) %
	exp( - 1.0 / (2.0 * sig*sig) * pow(y - mu, 2)) -
	pow(2.0 * PI, -0.5) * pow(sig, -2.0) *
	exp( - 1.0 / (2.0 * sig*sig) * pow(y - mu, 2.0));
}

