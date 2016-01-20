set.seed(123)
dat <- data.frame(X = rnorm(n = 100, mean = 100),
                  Y = rnorm(n = 100, mean = 100),
                  W = abs(rnorm(n = 100, mean = 0))
)

model <- list()
model[[1]] <- Y ~ a + b * X


# nls
a <- coef(nls(Y ~ a + b * X, data = dat, start = c(a=0, b=0)))
b <- coef(nls("Y ~ a + b * X", data = dat, start = c(a=0, b=0)))
c <- coef(nls(model[[1]], data = dat, start = c(a=0, b=0)))


A <- coef(nlsur(Y ~ a + b * X, data = dat))
B <- coef(nlsur("Y ~ a + b * X", data = dat))
C <- coef(nlsur(model, data = dat, startvalues = c(a = 0, b=0)))


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
