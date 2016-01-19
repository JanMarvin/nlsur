
set.seed(123)

X <- rnorm(n = 100, mean = 0)
Y <- rnorm(n = 100, mean = 100)


dat <- data.frame(X = X, Y = Y)

model <- list()
model[[1]] <- Y ~ a + b * X


nls("Y ~ a + b * X")

nls(Y ~ a + b * X)


nlsur(Y ~ a + b * X, data = dat,
      startvalues = c(a = 0, b=0))

nlsur("Y ~ a + b * X", data = dat,
      startvalues = c(a = 0, b=0))


nlsur(model, data = dat,
      startvalues = c(a = 0, b=0))
