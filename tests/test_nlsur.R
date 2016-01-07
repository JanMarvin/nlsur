library("nlsur")
dat  <- NULL
erg1 <- NULL
erg2 <- NULL
erg3 <- NULL

# inst/extdata/dat is petridish data from Stata 13.
source( system.file("extdata", "petridish.R", package = "nlsur") )


model <- list(p1 ~ b1 * b2 ^ time,
              p2 ~ g1 * g2 ^ time)

startvalues <- c(b1=1e-01, b2=1e-01,
                 g1=1e-01, g2=1e-01)

erg1 <- nlsur(eqns = model, startvalues = startvalues, data = dat, type = 1,
              eps = .Machine$double.eps, trace = TRUE)

erg2 <- nlsur(eqns = model, startvalues = startvalues, data = dat, type = 2,
              eps = .Machine$double.eps, trace = TRUE)

erg3 <- nlsur(eqns = model, startvalues = startvalues, data = dat, type = 3,
              eps = .Machine$double.eps, trace = TRUE)

# startvalues <- c(b1=0,b2=0,g1=0,g2=0)
# startvalues <- c(b1=1e-03, b2=1e-03,
#                  g1=1e-03, g2=1e-03)
# startvalues <- c(b1=1e-02, b2=1e-02,
#                  g1=1e-02, g2=1e-02)

# sequenz <- seq(0.01, 1, 0.1)
# for (a in sequenz) {
#   for (b in sequenz) {
#     for (c in sequenz) {
#       for (d in sequenz) {
#
#         startvalues <- c(b1 = a, b2 = b, g1 = c, g2 = d)
#
#         erg1 <- nlsur(eqns = model, startvalues = startvalues, data = dat,
#                       type = 1)
#         erg2 <- nlsur(eqns = model, startvalues = startvalues, data = dat,
#                       type = 2)
#         erg3 <- nlsur(eqns = model, startvalues = startvalues, data = dat,
#                       type = 3)
#         cat(coef(erg1), " || ", coef(erg2), " || ", coef(erg3), " || ", startvalues, "\n")
#
#       }
#     }
#   }
# }
#
#
# model1 <- list(p1 ~ b1 * b2 ^ time)
# model2 <- list(p2 ~ g1 * g2 ^ time)
#
# startvalues1 <- c(b1=1e-01, b2=1e-01)
# startvalues2 <- c(g1=1e-01, g2=1e-01)
#
# e1 <- nlsur(eqns = model1, startvalues = startvalues1, data = dat, type = 1)
# e2 <- nlsur(eqns = model2, startvalues = startvalues2, data = dat, type = 1)
#
# ea <- nls(model1[[1]], data = dat, start = startvalues1)
# eb <- nls(model2[[1]], data = dat, start = startvalues2)
