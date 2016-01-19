greene <- read.table(header = TRUE, sep  = ",", text = "
name,Coef,SE
bk,0.05682,0.00131
dkm,-0.02169,0.00963
bl,0.25355,0.001987
dll,0.07488,0.00639
be,0.04383,0.00105
dle,-0.00321,0.00275
bm,0.64580,0.00299
dlm,-0.07169,0.00941
dkk,0.02987,0.00575
dee,0.02938,0.00741
dkl,0.0000221,0.00367
dem,-0.01797,0.01075
dke,-0.00820,0.00406
dmm,0.11134,0.02239",
row.names = 1)

greene <- greene[order(row.names(greene)),]


greene

dat  <- NULL
erg1 <- NULL
erg2 <- NULL
erg3 <- NULL

# inst/extdata/dat is mfgcost data from Stata 13.
source( system.file("extdata", "mfgcoast.R", package = "nlsur") )

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

# erg1 <- nlsur(eqns = model, startvalues = startvalues, data = dat, type = 1,
#                eps = .Machine$double.eps, trace = TRUE)
# erg1

erg2 <- nlsur(eqns = model, startvalues = startvalues, data = dat, type = 2,
               eps = .Machine$double.eps, trace = TRUE)
erg2

# erg3 <- nlsur(eqns = model, startvalues = startvalues, data = dat, type = 3,
#                eps = .Machine$double.eps, trace = TRUE)
# erg3
#
# # print summary
# summary(erg1)
# summary(erg2)
# summary(erg3)

# indirect estimation of translog parameters
bm <- nlcom(object = erg2, form = "1 -be -bk -bl", rname= "bm")
# bm

dkm <- nlcom(object = erg2, form = "-dkk -dkl -dke", rname = "dkm")
# dkm

dlm <- nlcom(object = erg2, form = "-dkl -dll -dle", rname = "dlm")
# dlm

dem <- nlcom(object = erg2, form = "-dke -dle -dee", rname = "dem")
# dem

dmm <- nlcom(object = erg2, form = "-dkm -dlm -dem", rname = "dmm")
# dmm


# last one is equivalent to the longer form of:
# dmm <- nlcom(object = erg,
#  form = "-(-dkk -dkl -dke) -(-dkl -dll -dle) -(-dke -dle -dee)")
# dmm

est <- summary(erg2)$coefficients
ind <- rbind(bm, dkm, dlm, dem, dmm)

res <- rbind(est, ind)
res <- res[order(rownames(res)),]
res[,1:2]

round(unlist(res[, "Estimate"]), digits = 6) - greene$Coef
round(unlist(res[, "Std. Error"]), digits = 5) - greene$SE

