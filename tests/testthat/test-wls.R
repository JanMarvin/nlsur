# WLS test 1 -------------------------------------------------------------------

set.seed(123)
X <- matrix(abs(rnorm(450,2)), 150, 3)
Y <- matrix(abs(rnorm(150,2)), 150, 1)
W <- diag(x = 3)
qW <- qr.solve(qr(W))
w <- matrix(rep(rep(x = 1, 150), 3), ncol = 3)

erg.1 <- lm_gls(X = X, Y = Y, W = W, neqs = 3)

X <- cbind(
  X[1:50,],
  X[51:100,],
  X[101:150,]
)

Y <- cbind(
  Y[1:50],
  Y[51:100],
  Y[101:150]
)


erg.2 <- wls_est(x = X, r = Y, qS = qW, w = w,
                 sizetheta = 3, fullreg = 1, tol = 1e-10)

#### wls_est 1 ####
test_that("wls_est 1", {
  expect_equal(as.numeric(erg.1), as.numeric(erg.2))
})


# WLS test 2 -------------------------------------------------------------------

set.seed(123)
X <- matrix(abs(rnorm(450,2)), 150, 3)
Y <- matrix(abs(rnorm(150,2)), 150, 1)
W <- matrix(c(2, 1, 0, 1, 2, 1, 0, 1, 2), 3, 3)
qW <- qr.solve(qr(W))
w <- matrix(rep(rep(x = 1, 150), 3), ncol = 3)

erg.1 <- lm_gls(X = X, Y = Y, W = W, neqs = 3)

X <- cbind(
  X[1:50,],
  X[51:100,],
  X[101:150,]
)

Y <- cbind(
  Y[1:50],
  Y[51:100],
  Y[101:150]
)

erg.2 <- wls_est(x = X, r = Y, qS = qW, w = w,
                 sizetheta = 3, fullreg = 1, tol = 1e-10)

#### calc_robust ####
test_that("wls_est 2", {
  expect_equal(as.numeric(erg.1), as.numeric(erg.2))
})


# Cov test 1 -------------------------------------------------------------------

set.seed(123)
X <- matrix(abs(rnorm(450,2)), 150, 3)
Y <- matrix(abs(rnorm(150,2)), 150, 1)
W <- diag(x = 3)
qW <- qr.solve(qr(W))
w <- matrix(rep(rep(x = 1, 150), 3), ncol = 3)

erg.1 <- lm_gls(X = X, Y = Y, W = W, neqs = 3, covb = 1)

X <- cbind(
  X[1:50,],
  X[51:100,],
  X[101:150,]
)

Y <- cbind(
  Y[1:50],
  Y[51:100],
  Y[101:150]
)

erg.2 <- wls_est(x = X, r = Y, qS = qW, w = w,
                 sizetheta = 3, fullreg = 0, tol = 1e-10)

#### cov test 1 ####
test_that("cov_est 1", {
  expect_equal(as.numeric(erg.1), as.numeric(erg.2))
})
