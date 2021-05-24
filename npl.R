# NPL for Gaussian Mixture model, see: Fong et al 2019
library(mvtnorm)
library(dist)		# devtools::install_github("mkomod/dist")
Rcpp::sourceCpp("./npl.cpp")
set.seed(1)

# Test data
p <- c(0.3, 0.3, 0.4)
mu <- c(-3, 0, 3)
sig <- c(0.2, 0.2, 0.2)

n <- 100
g <- sample(1:3, n, prob=p, replace=T)
y <- rnorm(n, mu[g], sig[g])
plot(density(y))


# NPL
gmm_ll <- function(y, p, mu, sig) 
{
    res <- sapply(seq_along(p), function(k) p[k] * dnorm(y, mu[k], sig[k]))
    return(-log(sum(res)))
}


J <- function(y, w, p, mu, sig) 
{
    l <- matrix(sapply(y, function(y) gmm_ll(y, p, mu, sig)), ncol=1)
    return(as.double(w %*% l))
}

J.grad <- function(y, w, p, mu, sig) 
{
    denom <- sum(sapply(seq_along(p), function(k) p[k] * dnorm(y, mu[k], sig[k])))
    d.p <- sapply(seq_along(p), function(k) {
	return(matrix(-dnorm(y, mu[k], sig[k]), ncol=1) / denom)
    })
    d.mu <- sapply(seq_along(mu), function(k) {
	return(-p[k] * dnorm_dMu(y, mu[k], sig[k]) / denom)
    })
    d.sig <- sapply(seq_along(sig), function(k) {
	return(-p[k] * dnorm_dSig(y, mu[k], sig[k]) / denom)
    })
    return(c(w %*% d.p, w %*% d.mu, w %*% d.sig))
}

J.fn <- function(pr, y, w) J(y, w, pr[1:3]/sum(pr[1:3]), pr[4:6], pr[7:9])
J.gr <- function(pr, y, w) J.grad(y, w, pr[1:3]/sum(pr[1:3]), pr[4:6], pr[7:9])

a <- 1
a.t <- 100
B <- 5000

res <- sapply(1:B, function(i) {
    w <- rdirichlet(1, c(rep(a, n), rep(a/a.t, a.t)))
    y.c <- c(y, runif(a.t, -6, 6))
    p.init <- c(runif(3), rnorm(3, sd=3), rgamma(3, 1, 1))
    p.opt <- tryCatch({
	p.opt <- optim(p.init,fn=J.fn, gr=J.gr, method="L-BFGS-B",
	      lower=c(5e-2, 5e-2, 5e-2, -100, -100, -100, 1e-2, 1e-2, 1e-2),
	      upper=c(3, 3, 3,  100,  100,  100, 10, 10, 10), y=y.c, w=w)$par
	p.opt[1:3] <- p.opt[1:3]/sum(p.opt[1:3])
	return(p.opt)
    }, error = function(e) return(rep(NA, 9)))
    return(p.opt)
})


