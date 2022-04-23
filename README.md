# NLSUR ![R-CMD-check](https://github.com/JanMarvin/nlsur/workflows/R-CMD-check/badge.svg)

`nlsur` is a package to estimate a nonlinear least squares for single equations
or systems of equations.
The function to interact with is `nlsur()`. This function can estimate Nonlinear-Least Squares (NLS), Feasible Generalized NLS (FGNLS) and Iterative FGNLS (IFGNLS).

The packages supports a variety of functions like `print()`, `coef()`, `summary()`, `logLik()`, `vcov()` and `predict()`.

## Installation

With `drat`:
```R
drat::addRepo("JanMarvin")
install.packages("nlsur")
```

With `devtools`:
```R
devtools::install_github("JanMarvin/nlsur")
```

## Application

With `nlsur()` it is rather straight forward to estimate nonlinear demand systems. As example the following Translog demand system can be estimated.

```R
data(costs)

dat <- costs
# apply a patch to create Greene Ed. 7 Data
dat$Sm[dat$Year == 1958] <- 0.61886
dat$Pe[dat$Year == 1950] <- 1.12442
dat$Pm[dat$Year == 1949] <- 1.06625

# model
model <- list(
  Sk ~ bk + dkk * log(Pk/Pm) + dkl * log(Pl/Pm) + dke * log(Pe/Pm),
  Sl ~ bl + dkl * log(Pk/Pm) + dll * log(Pl/Pm) + dle * log(Pe/Pm),
  Se ~ be + dke * log(Pk/Pm) + dle * log(Pl/Pm) + dee * log(Pe/Pm)
)

erg <- nlsur(eqns = model, data = dat, type = "FGNLS")
erg
```

Additional parameters may be obtained using `nlcom()` a wrapper for the delta method.

```R
# indirect estimation of translog parameters
bm  <- nlcom(object = erg, form = "1 -be -bk -bl",  rname = "bm")
dkm <- nlcom(object = erg, form = "-dkk -dkl -dke", rname = "dkm")
dlm <- nlcom(object = erg, form = "-dkl -dll -dle", rname = "dlm")
dem <- nlcom(object = erg, form = "-dke -dle -dee", rname = "dem")

# and now dmm (nlcom can search for parameters)
dmm <- nlcom(object = erg, form = "-dkm -dlm -dem", rname = "dmm")
```

