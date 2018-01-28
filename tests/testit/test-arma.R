# Test data set
set.seed(123)

xi <- list()
xi[[1]] <- matrix( rnorm(225), ncol = 9, nrow = 25)
xi[[2]] <- matrix( rnorm(225), ncol = 9, nrow = 25)
xi[[3]] <- matrix( rnorm(225), ncol = 9, nrow = 25)

colnames(xi[[1]]) <- colnames(xi[[2]]) <- colnames(xi[[3]]) <- letters[1:9]


ri <- list()
ri[[1]] <- matrix( rnorm(25), ncol = 1)
ri[[2]] <- matrix( rnorm(25), ncol = 1)
ri[[3]] <- matrix( rnorm(25), ncol = 1)

theta <- colnames(xi[[1]])

x <- do.call(cbind, xi)
X <- do.call(rbind, xi); colnames(X) <- theta
r <- do.call(cbind, ri)
R <- matrix(r, ncol = 1)
n <- nrow(r)
k <- ncol(r)
w <- rep(1, n)

I <- diag(n)

s <- 1/n * crossprod(r)
qS <- qr.solve(s)

S <- s %x% I
W <- qS %x% I

XDX <- t(X) %*% W %*% X
XDY <- t(X) %*% W %*% R

# calc_ssr
rss_a <- ssr_est(r = r, s = s, w = w)

rss_m <- crossprod( S %*% R)


# reshape
mm <- matrix(c(11,21,31,41,
               12,22,32,42,
               13,23,33,43,
               14,24,34,44),
             ncol = 4)

mm_a <- arma_reshape(mm, 2)

mm_m <- matrix(t(mm), nrow = 2, byrow =T )

# wls
BB <- t(qr.coef(qr(XDX), XDY))
# Bb <- coef(lm.gls(R ~ 0 + X, W = W)) # requies MASS
bb <- wls_est(x = x, r = r, qS = qS,
              w = w, sizetheta = length(theta),
              fullreg = 1, tol = .Machine$double.eps)

XDX_a <- wls_est(x = x, r = r, qS = qS,
                 w = w, sizetheta = length(theta),
                 fullreg = 0, tol = .Machine$double.eps)

dimnames(XDX_a) <- list(theta, theta)

BB <- as.vector(BB)
# Bb <- as.vector(Bb)
bb <- as.vector(bb)

# wt_mean
wtm_r <- weighted.mean(x = X[, "a"], w = X[, "b"])
wtm_a <- wt_mean(x = X[, "a"], w = X[, "b"])

# calc_robust

data( costs )
# library(sandwich)
# res <- lm("Cost ~ Sk + Sl + Se", data = costs)
# sandwich_se <- sqrt(diag(vcovHC(res, type = "HC0")))
sandwich_se <- c(355.081244412386,
                 3161.61958708514,
                 1095.47069852517,
                 3296.21760057729)


eqs <- "Cost ~ b0 + b1 * Sk  + b2 * Sl  + b3 * Se"
nes <- nlsur(eqns = eqs, data = costs,
             type = "NLS", stata = TRUE, robust = TRUE)

nlsur_se <- summary(nes)$coefficients[,2]
names(nlsur_se) <- NULL

#### calc_ssr ####
assert("ssr_est", {
  all.equal(rss_a, rss_m)
})

#### arma_reshape ####
assert("arma_reshape", {
  all.equal(mm_a, mm_m)
})

#### calc_reg ####
assert("wls_est", {
  all.equal(bb, BB)
  # all.equal(bb, Bb)
  all.equal(XDX_a, XDX)
})

#### wt_mean ####
assert("wt_mean", {
  all.equal(wtm_r, wtm_a)
})

#### calc_robust ####
assert("cov_robust", {
  all.equal(sandwich_se, nlsur_se)
})
