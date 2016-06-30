
#' Estimate nonlinear combinations of nlsur estimates
#'
#' @param object of class nlsur
#' @param form formula e.g. "be/bk". May contain names(coef(object)) or prior
#' nlcom estimations.
#' @param alpha value for conf. interval
#' @param rname optional rowname for result
#' @importFrom car deltaMethod
#' @importFrom stats coef formula pt qnorm
#' @seealso \link{deltaMethod}
#' @examples
#' \dontrun{
#' dkm <- nlcom(object = erg, form = "-dkk -dkl -dke", rname = "dkm")
#' dkm
#'
#' dlm <- nlcom(object = erg, form = "-dkl -dll -dle", rname = "dlm")
#' dlm
#'
#' dem <- nlcom(object = erg, form = "-dke -dle -dee", rname = "dem")
#' dem
#'
#' dmm <- nlcom(object = erg, form = "-dkm -dlm -dem", rname = "dmm")
#' dmm
#'
#' # last one is equivalent to the longer form of:
#' dmm <- nlcom(object = erg,
#'  form = "-(-dkk -dkl -dke) -(-dkl -dll -dle) -(-dke -dle -dee)")
#' dmm
#' }
#'
#' @export
nlcom <- function(object, form, alpha, rname) {

  # store original form
  oform <- form

  # check if all form vars are part of object
  vars <- all.vars(formula(paste(form, "~ 0")))
  vars <- vars[!(vars %in% names(coef(object)))]

  # if form contains other nlcom estimates, fetch their formula and update
  # the current form with the older formula.
  if( !identical(vars, character(0)) ){
    for (i in vars){
      tmp <- get0(x = i)

      if(!(class(tmp) == "nlcom"))
        stop(paste(i, "is not of class nlcom."))

      fname <- paste("(", attr(tmp, "form"), ")")

      form <- gsub(pattern = i, replacement = fname, x = form, fixed = TRUE)

      }
  }

  z <- deltaMethod(object, form)

  tval <- z$Estimate / z$SE

  nE <- sum(object$n)
  kE <- sum(object$k)

  prob <- 2 * (1 - pt(abs(tval), (nE * kE )))

  z   <- cbind(z, tval, prob)

  colnames(z) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)" )

  lk <- length(coef(object))
  neqs <- length(object$n)

  if (missing (alpha))
    alpha <- 0.05

  alphaz <- c( alpha/2, 1-alpha/2)

  # add conficence intervals
  cinterv <- z$Estimate  + qnorm(alphaz) * ( z$`Std. Error`/ sqrt (lk) * neqs )
  names(cinterv) <- as.character(alphaz)

  if(missing(rname))
    attr(z, "rname")   <- oform
  else
    attr(z, "rname")   <- rname

  attr(z, "form")    <- form
  attr(z, "oform")   <- oform
  attr(z, "Confint") <- cinterv
  class(z)           <- "nlcom"

  z
}

#' @export
print.nlcom <- function(x, ...) {

  mat <- matrix(unlist(x), nrow = 1)
  dimnames(mat) <- list(attr(x, "rname"),
                        names(unlist(x)))

  cat("nlcom: ", attr(x, "oform"), "\n\n")

  printCoefmat(mat, signif.legend = FALSE , ...)

  if (attr(x, "oform") != attr(x, "form"))
    cat("\nnlcom object estimated with prior nlcom estimate.\n")

}
