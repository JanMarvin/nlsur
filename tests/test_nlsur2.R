library("nlsur")
dat  <- NULL
erg1 <- NULL
erg2 <- NULL
erg3 <- NULL

# inst/extdata/dat is mfgcost data from Stata 13.
source("~/Source/nlsur/inst/extdata/mfgcoast.R")

# model
model <- list(
  (s_k ~ bk + dkk*log(pk/pm) + dkl*log(pl/pm) + dke*log(pe/pm)),
  (s_l ~ bl + dkl*log(pk/pm) + dll*log(pl/pm) + dle*log(pe/pm)),
  (s_e ~ be + dke*log(pk/pm) + dle*log(pl/pm) + dee*log(pe/pm))
)

# Function to create 0 values for and model
getstartvals <- function(model, data) {
  # automatic creation of start values
  modelparameters <- unlist(lapply(model, all.vars))
  svals <- modelparameters[which(!modelparameters %in% names(data))]
  svals <- unique(svals[order(svals)])
  # svals
  strtvls <- rep(0, length(svals))
  names(strtvls) <- svals
  return(strtvls)
}

# startvalues
startvalues <- getstartvals(model = model, data = dat)
print(startvalues)

erg1 <- ifgnls(eqns = model, startvalues = startvalues, data = dat,
               trace = FALSE)
erg1

erg2 <- ifgnls(eqns = model, startvalues = startvalues, data = dat, type = 2,
               trace = FALSE)
erg2

erg3 <- ifgnls(eqns = model, startvalues = startvalues, data = dat, type = 3,
               trace = FALSE)
erg3

# print summary
summary(erg1)
summary(erg2)
summary(erg3)
