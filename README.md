# NLSUR

`nlsur` ist a package to estimate a nonlinear least squares for single equations
or systems of equations.
The function to interact with is `nlsur()`. This function can estimate Nonlinear-Least Squares (NLS), Feasible Generalized NLS (FGNLS) and Iterative FGNLS (IFGNLS).

# Installation

```{r}
devtools::install_github("JanMarvin/nlsur")
```

# Application

With `nlsur` the following Translog demand system can be estimated.

```{r}
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

# Status

[![Build Status](https://travis-ci.org/JanMarvin/nlsur.svg?branch=master)](https://travis-ci.org/JanMarvin/nlsur)
