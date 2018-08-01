set.seed(123)
dat <- data.frame(X = rnorm(n = 100, mean = 100),
                  Y = rnorm(n = 100, mean = 100),
                  W = abs(rnorm(n = 100, mean = 0))
)

model <- list()
model[[1]] <- Y ~ a + b * X


# nls

m1 <- nls(Y ~ a + b * X, data = dat, start = c(a=0, b=0))
m2 <- nls("Y ~ a + b * X", data = dat, start = c(a=0, b=0))
m3 <- nls(model[[1]], data = dat, start = c(a=0, b=0))

a <- coef(m1)
b <- coef(m2)
c <- coef(m3)


M1 <- nlsur(Y ~ a + b * X, data = dat)
M2 <- nlsur("Y ~ a + b * X", data = dat)
M3 <- nlsur(model, data = dat, startvalues = c(a = 0, b=0))

A <- coef(M1)
B <- coef(M2)
C <- coef(M3)


# weighted nls
aw <- coef(nls(Y ~ a + b * X, data = dat, start = c(a=0, b=0),
               weights = W))
bw <- coef(nls("Y ~ a + b * X", data = dat, start = c(a=0, b=0),
               weights = W))
cw <- coef(nls(model[[1]], data = dat, start = c(a=0, b=0),
               weights = W))

Aw <- coef(nlsur(Y ~ a + b * X, data = dat,
                 weights = W))
Bw <- coef(nlsur("Y ~ a + b * X", data = dat,
                 weights = W))
Cw <- coef(nlsur(model, data = dat, startvalues = c(a = 0, b=0),
                 weights = W))


# predict
pm1 <- predict(m1)
pM1 <- predict(M1)

#### nls ####
test_that("nls", {
  expect_equal(a, A)
  expect_equal(b, B)
  expect_equal(c, C)
})

#### weighted nls ####
test_that("weighted nls", {
  expect_equal(aw, Aw)
  expect_equal(bw, Bw)
  expect_equal(cw, Cw)
})

#### predict ####
test_that("predict", {
  expect_equal(pm1, c(pM1)[["Y"]])
})

