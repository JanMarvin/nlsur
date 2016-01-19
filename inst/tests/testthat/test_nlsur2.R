# library("nlsur")
# dat  <- NULL
# erg1 <- NULL
# erg2 <- NULL
# erg3 <- NULL
#
# # inst/extdata/dat is mfgcost data from Stata 13.
# source( system.file("extdata", "mfgcoast.R", package = "nlsur") )
#
# # model
# model <- list(
#   (s_k ~ bk + dkk*log(pk/pm) + dkl*log(pl/pm) + dke*log(pe/pm)),
#   (s_l ~ bl + dkl*log(pk/pm) + dll*log(pl/pm) + dle*log(pe/pm)),
#   (s_e ~ be + dke*log(pk/pm) + dle*log(pl/pm) + dee*log(pe/pm))
# )
#
# # Function to create 0 values for and model
# getstartvals <- function(model, data) {
#   # automatic creation of start values
#   modelparameters <- unlist(lapply(model, all.vars))
#   svals <- modelparameters[which(!modelparameters %in% names(data))]
#   svals <- unique(svals[order(svals)])
#   # svals
#   strtvls <- rep(0, length(svals))
#   names(strtvls) <- svals
#   return(strtvls)
# }
#
# # startvalues
# startvalues <- getstartvals(model = model, data = dat)
# print(startvalues)
#
# erg1 <- nlsur(eqns = model, startvalues = startvalues, data = dat, type = 1,
#                eps = .Machine$double.eps, trace = TRUE)
# erg1
#
# erg2 <- nlsur(eqns = model, startvalues = startvalues, data = dat, type = 2,
#                eps = .Machine$double.eps, trace = TRUE)
# erg2
#
# erg3 <- nlsur(eqns = model, startvalues = startvalues, data = dat, type = 3,
#                eps = .Machine$double.eps, trace = TRUE)
# erg3
#
# # print summary
# summary(erg1)
# summary(erg2)
# summary(erg3)
