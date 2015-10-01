library("nlsur")
dat  <- NULL
erg1 <- NULL
erg2 <- NULL
erg3 <- NULL

# inst/extdata/dat is petridish data from Stata 13.
source("inst/extdata/petridish.R")


model <- list(p1 ~ b1 * b2 ^ time,
              p2 ~ g1 * g2 ^ time)

startvalues <- c(b1=1e-01, b2=1e-01,
                g1=1e-01, g2=1e-01)

startvalues <- c(b1=0,b2=0,g1=0,g2=0)
erg1 <- nlsur(eqns = model, startvalues = startvalues, data = dat, nls = TRUE,
              debug = TRUE)

erg1 <- ifgnls(eqns = model, startvalues = startvalues, data = dat, type = 1,
               trace = T)
erg2 <- ifgnls(eqns = model, startvalues = startvalues, data = dat, type = 2,
               trace = T)
erg3 <- ifgnls(eqns = model, startvalues = startvalues, data = dat, type = 3,
               trace = T)

# startvalues <- c(b1=0,b2=0,g1=0,g2=0)
# startvalues <- c(b1=1e-03, b2=1e-03,
#                  g1=1e-03, g2=1e-03)
# startvalues <- c(b1=1e-02, b2=1e-02,
#                  g1=1e-02, g2=1e-02)

sequenz <- seq(0.01, 1, 0.1)
for (a in sequenz) {
  for (b in sequenz) {
    for (c in sequenz) {
      for (d in sequenz) {

        startvalues <- c(b1 = a, b2 = b, g1 = c, g2 = d)

        erg1 <- ifgnls(eqns = model, startvalues = startvalues, data = dat,
                       type = 1)
        cat(coef(erg1), " || ", startvalues, "\n")

      }
    }
  }
}

