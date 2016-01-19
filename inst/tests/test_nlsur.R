# library("nlsur")
# dat  <- NULL
# erg1 <- NULL
# erg2 <- NULL
# erg3 <- NULL
#
# # inst/extdata/dat is petridish data from Stata 13.
# source( system.file("extdata", "petridish.R", package = "nlsur") )
#
#
# model <- list(p1 ~ b1 * b2 ^ time,
#               p2 ~ g1 * g2 ^ time)
#
# startvalues <- c(b1=1e-01, b2=1e-01,
#                  g1=1e-01, g2=1e-01)
#
# erg1 <- nlsur(eqns = model, startvalues = startvalues, data = dat, type = 1,
#               eps = .Machine$double.eps, trace = TRUE)
#
# erg2 <- nlsur(eqns = model, startvalues = startvalues, data = dat, type = 2,
#               eps = .Machine$double.eps, trace = TRUE)
#
# erg3 <- nlsur(eqns = model, startvalues = startvalues, data = dat, type = 3,
#               eps = .Machine$double.eps, trace = TRUE)
