
# set enviroment
env <- new.env()

# compare greene with nlsur
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

dat  <- NULL
erg2 <- NULL

bm <- dkm <- dlm <- dem <- dmm <- NULL

# use costs data
data(costs)

dat <- costs
# apply a patch to create Greenes Ed. 7 Data
dat$Sm[dat$Year == 1958] <- 0.61886
dat$Pe[dat$Year == 1950] <- 1.12442
dat$Pm[dat$Year == 1949] <- 1.06625

# model

model <- list(
  Sk ~ bk + dkk * log(Pk/Pm) + dkl * log(Pl/Pm) + dke * log(Pe/Pm),
  Sl ~ bl + dkl * log(Pk/Pm) + dll * log(Pl/Pm) + dle * log(Pe/Pm),
  Se ~ be + dke * log(Pk/Pm) + dle * log(Pl/Pm) + dee * log(Pe/Pm)
)

e1   <- nlsur(eqns = model, data = dat, type = "NLS", stata = FALSE,
              qrsolve = FALSE)
e2   <- nlsur(eqns = model, data = dat, type = "NLS", stata = FALSE,
              qrsolve = TRUE)


E1   <- nlsur(eqns = model, data = dat, type = "IFGNLS",
              MASS = FALSE)
E2   <- nlsur(eqns = model, data = dat, type = "IFGNLS",
              MASS = TRUE)


erg2 <- nlsur(eqns = model, data = dat, type = "FGNLS")
erg2

# indirect estimation of translog parameters
bm <- nlcom(object = erg2, form = "1 -be -bk -bl", rname= "bm")
# bm

dkm <- nlcom(object = erg2, form = "-dkk -dkl -dke", rname = "dkm")
# dkm

dlm <- nlcom(object = erg2, form = "-dkl -dll -dle", rname = "dlm")
# dlm

dem <- nlcom(object = erg2, form = "-dke -dle -dee", rname = "dem")
# dem

# # testthat does not like my get.
dmm <- nlcom(object = erg2, form = "-dkm -dlm -dem", rname = "dmm", envir = env)
# dmm <- nlcom(object = erg2,
#              form = "-(-dkk -dkl -dke) -(-dkl -dll -dle) -(-dke -dle -dee)")
# dmm

est <- summary(erg2)$coefficients
ind <- rbind(bm, dkm, dlm, dem, dmm)

res <- rbind(est, ind)
res <- res[order(rownames(res)),]


res_n <- round(unlist(res[, "Estimate"]), digits = 5)
res_g <- round(greene$Coef, digits = 5)

names(res_g) <- rownames(greene)


test_that("Compare nlsur options", {
  expect_equal(coef(e1), coef(e2))
  # not exact due to numericDeriv
  expect_equal(round(coef(E1),5), round(coef(E2),5))
})

test_that("Compare Greene to nlsur", {
  expect_equal(res_n, res_g)
})
