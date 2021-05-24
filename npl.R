# NPL for Gaussian Mixture model, see: Fong et al 2019
library(dist)		# devtools::install_github("mkomod/dist")
Rcpp::sourceCpp("./npl.cpp")
set.seed(1)

# Test data
p <- c(0.1, 0.3, 0.7)
mu <- c(0, 2, 4)
sig <- c(1, 1, 1)

n <- 1000
g <- sample(1:3, n, prob=p, replace=T)
y <- rnorm(n, mu[g], sig[g])
w <- rdirichlet(1, rep(1, n))
plot(density(y))

a <- 1
a.t <- 500
B <- 10000
res <- sapply(1:B, function(i) {
    w <- dist::rdirichlet(1, c(rep(a, n), rep(a/a.t, a.t)))
    y.c <- c(y, runif(a.t, -2, 6))
    p.init <- c(rdirichlet(1 , rep(1,3)), rnorm(3, sd=3), rgamma(3, 1, 1))
    p.opt <- tryCatch({
	p.opt <- optim(p.init,fn=J_fn, gr=J_gr, method="L-BFGS-B",
	      lower=c(5e-2, 5e-2, 5e-2, -100, -100, -100, 5e-2, 5e-2, 5e-2),
	      upper=c(  10,   10,   10,  100,  100,  100,   10,   10,   10), 
	      y=y.c, w=w, k=3)$par
	p.opt[1:3] <- p.opt[1:3]/sum(p.opt[1:3])
	return(p.opt)
    }, error = function(e) { print(e); return(rep(NA, 9)) })
    print(i)
    return(p.opt)
})

res.fil <- res[ , apply(res, 2, function(r) !any(which(is.na(r))))]
apply(res.fil, 1, mean)
f1 <- MASS::kde2d(res.fil[4, ], res.fil[5, ], lims=c(-2,6,-2,6), n = 2000)
f2 <- MASS::kde2d(res.fil[5, ], res.fil[6, ], lims=c(-2,6,-2,6), n = 2000)
png("density.png", width=600, height=600)
image(f1, xlim=c(-2, 6), ylim=c(-2, 6), col=hcl.colors(12))
dev.off()
image(f2, xlim=c(-2, 6), ylim=c(-2, 6), col=hcl.colors(12))

