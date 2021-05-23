# NPL for Gaussian Mixture model, see: Fong et al 2019
library(mvtnorm)
set.seed(1)

# Test data
p <- c(0.1, 0.5, 0.4)
mu <- c(-3, 0, 3)
sig <- c(0.2, 0.2, 0.2)

n <- 100
g <- sample(1:3, n, prob=p, replace=T)
y <- rnorm(n, mu[g], sig[g])
plot(density(y))


# NPL
gmm_ll <- function(y, p, mu, sig)
{
    # approximation of the log likelihood
    res <- sapply(seq_along(p), function(k) 
	   log(p[k]) + sum(dnorm(y, mu[k], sig[k], log=T)))
    return(-max(res))
}

gmm_ll(y, p, mu, sig)

