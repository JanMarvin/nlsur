require(testthat)


library(MASS)

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

W <- qS %x% I

XDX <- t(X) %*% W %*% X
XDY <- t(X) %*% W %*% R

BB <- t(qr.solve(XDX, XDY))

erg <- lm.gls(R ~ 0 + X, W = W)
Bb <- coef(erg)

bb <- calc_reg(x = x, r = r, qS = qS, w = w, sizetheta = length(theta), fullreg = 1)

names(bb) <- NULL
names(Bb) <- NULL
names(BB) <- NULL

bb <- as.vector(bb)
Bb <- as.vector(Bb)
BB <- as.vector(BB)


# Test
test_that("WLS equals Matrix and lm.gls variant", {
  expect_that(bb, equals(BB))
  expect_that(bb, equals(Bb))
})
