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
wm1 <- nls(Y ~ a + b * X, data = dat, start = c(a=0, b=0),
           weights = W)
wm2 <- nls("Y ~ a + b * X", data = dat, start = c(a=0, b=0),
           weights = W)
wm3 <- nls(model[[1]], data = dat, start = c(a=0, b=0),
           weights = W)

aw <- coef(wm1)
bw <- coef(wm2)
cw <- coef(wm3)

wM1 <- nlsur(Y ~ a + b * X, data = dat,
             weights = "W")
wM2 <- nlsur("Y ~ a + b * X", data = dat,
             weights = "W")
wM3 <- nlsur(model, data = dat, startvalues = c(a = 0, b=0),
             weights = "W")

Aw <- coef(wM1)
Bw <- coef(wM2)
Cw <- coef(wM3)


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


#### logLik ####
test_that("logLik", {
  expect_equal(as.numeric(logLik(m1)), as.numeric(logLik(M1)))
  expect_equal(as.numeric(logLik(m2)), as.numeric(logLik(M2)))
  expect_equal(as.numeric(logLik(m3)), as.numeric(logLik(M3)))
  expect_equal(as.numeric(logLik(wm1)), as.numeric(logLik(wM1)))
  expect_equal(as.numeric(logLik(wm2)), as.numeric(logLik(wM2)))
  expect_equal(as.numeric(logLik(wm3)), as.numeric(logLik(wM3)))
})
