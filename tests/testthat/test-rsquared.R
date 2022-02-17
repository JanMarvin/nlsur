#
# Test r-squared values of nlsur and lm
#
#


set.seed(123)

dd <- data.frame(y = rnorm(n = 100, mean = 5, sd = 5),
                 x = rnorm(n = 100, mean = 2, sd = 5),
                 w = sample(x = seq(0.1, 1, 0.1), size = 100, replace = TRUE))


model <- "y ~ b * x"

# res
res_a <- summary(lm("y~0+x", data = dd))
res_1 <- summary(nls(formula = model, data = dd, start = c(b = 0)))
res_2 <- summary(nlsur(eqns = model, data = dd))

x1.1 <- res_a$r.squared
x1.2 <- res_2$zi["R-squared"]

# res weighted
res_a <- summary(lm("y~0+x", data = dd, weights = w))
res_1 <- summary(nls(formula = model, data = dd, weights = w, start = c(b = 0)))
res_2 <- summary(nlsur(eqns = model, data = dd, weights = w, type = 1, stata = FALSE))

x1.3 <- res_a$r.squared
x1.4 <- res_2$zi["R-squared"]



#### model with constant ####
model <- "y ~ a + b * x"

# res
res_a <- summary(lm("y~x", data = dd))
res_1 <- summary(nls(formula = model, data = dd, start = c(a = 0, b = 0)))
res_2 <- summary(nlsur(eqns = model, data = dd))

x2.1 <- res_a$r.squared
x2.2 <- res_2$zi["R-squared"]

# res weighted
res_a <- summary(lm("y~x", data = dd, weights = w))
res_1 <- summary(nls(formula = model, data = dd, weights = w, start = c(a = 0, b = 0)))
res_2 <- summary(nlsur(eqns = model, data = dd, weights = w, type = 1, stata = FALSE))

x2.3 <- res_a$r.squared
x2.4 <- res_2$zi["R-squared"]



#### run test ####

test_that("rsquared", {
  expect_equal(as.numeric(x1.1), as.numeric(x1.2))
  expect_equal(as.numeric(x1.3), as.numeric(x1.4))
  expect_equal(as.numeric(x2.1), as.numeric(x2.2))
  expect_equal(as.numeric(x2.3), as.numeric(x2.4))
})
