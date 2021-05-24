#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec dnorm_dMu(const arma::vec &y, double mu, double sig)
{
    return pow(2.0 * PI, -0.5) * (y - mu) * pow(sig, -3.0) %
	exp( - 1.0/(2.0 * sig*sig) * pow(y - mu, 2.0));
}

// [[Rcpp::export]]
arma::vec dnorm_dSig(const arma::vec &y, double mu, double sig)
{
    return pow(2.0 * PI, -0.5) * pow(y - mu, 2.0) * pow(sig, -4.0) %
	exp( - 1.0 / (2.0 * sig*sig) * pow(y - mu, 2)) -
	pow(2.0 * PI, -0.5) * pow(sig, -2.0) *
	exp( - 1.0 / (2.0 * sig*sig) * pow(y - mu, 2.0));
}

arma::vec dnorm(const arma::vec &y, double m, double s) 
{
    arma::vec res = arma::vec(y.n_rows, arma::fill::zeros);
    for (int i = 0; i < y.n_rows; ++i) {
	res(i) = R::dnorm(y(i), m, s, 0);
    }
    return res;
}

// Gaussian mixute model likelihood
arma::vec Gmm_l(const arma::vec &y, arma::vec p, arma::vec m, arma::vec s)
{
    arma::vec res = arma::vec(y.n_rows, arma::fill::zeros);
    for (int i = 0; i < y.n_rows; ++i) {
	for (int k = 0; k < p.n_rows; ++k) {
	    res(i) += p(k) * R::dnorm(y(i), m(k), s(k), 0);
	}
    }
    return res;
}

// [[Rcpp::export]]
double J_fn(arma::vec params, const arma::vec &y, const arma::vec &w, int k)
{
    arma::vec p = params.rows(0, k-1);
    p = p / sum(p); // normalise p
    arma::vec m = params.rows(k, 2*k - 1);
    arma::vec s = params.rows(2*k, 3*k -1);
    arma::vec y_ll = -log(Gmm_l(y, p, m, s));
    double res = arma::dot(w, y_ll);
    return res;
}

// [[Rcpp::export]]
arma::vec J_gr(arma::vec params, const arma::vec &y, const arma::vec &w, int k) 
{
    arma::vec p = params.rows(0, k-1);
    p = p / sum(p); // normalise p
    arma::vec m = params.rows(k, 2*k - 1);
    arma::vec s = params.rows(2*k, 3*k -1);
    arma::vec y_l =  Gmm_l(y, p, m, s);

    arma::vec deriv = arma::vec(3 * k, arma::fill::zeros);

    for (int i = 0; i < k; ++i) {
	// derivatives of pi
	arma::vec dp = dnorm(y, m(i), s(i));
	deriv(i) = dot(w, -dp/y_l);
	
	// derivatives of mu
	arma::vec dmu = - p(i) * dnorm_dMu(y, m(i), s(i));
	deriv(k + i) = dot(w, dmu / y_l);

	// derivatives of sig
	arma::vec ds = - p(i) * dnorm_dSig(y, m(i), s(i));
	deriv(2*k + i) = dot(w, ds / y_l);
    }

    return deriv;
}


