
#' Estimate nonlinear combinations of nlsur estimates
#'
#' @param object of class nlsur
#' @param form formula e.g. "be/bk"
#' @param alpha value for conf. interval
#' @importFrom car deltaMethod
#'
#' @export
nlcom <- function(object, form, alpha) {
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


  attr(z, "Confint") <- cinterv
  attr(z, "rname")   <- form
  class(z)           <- "nlcom"

  z
}

#' @export
print.nlcom <- function(object, ...) {

  mat <- matrix(unlist(object), nrow = 1)
  dimnames(mat) <- list(attr(object, "rname"),
                        names(unlist(object)))

  printCoefmat(mat, signif.legend = FALSE , ...)

}
