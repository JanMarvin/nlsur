# Copyright (c) 2017 Jan Marvin Garbuszus
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

#' Estimate nonlinear combinations of nlsur estimates
#'
#' @param object of class nlsur
#' @param form formula e.g. "be/bk". May contain names(coef(object)) or prior
#' nlcom estimations.
#' @param alpha value for conf. interval default is 0.05
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
nlcom <- function(object, form, alpha = 0.05, rname) {

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

      if(class(tmp) == "nlcom") {
        fname <- paste("(", attr(tmp, "form"), ")")

        form <- gsub(pattern = i, replacement = fname, x = form, fixed = TRUE)
      }
    }
  }

  alevel <- 1 - alpha

  z <- deltaMethod(object, form, level = alevel)

  # separate the convinterval from the z output
  cinterv <- z[-c(1:2)]

  # only keep est and se
  z <- z[c("Estimate", "SE")]

  # calculate t-value
  tval <- z$Estimate / z$SE

  # calculate prob
  nE <- sum(object$n); kE <- sum(object$k)
  prob <- 2 * (1 - pt(abs(tval), (nE * kE )))

  # create output
  z   <- cbind(z, tval, prob)
  colnames(z) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)" )

  lk <- length(coef(object))
  neqs <- length(object$n)

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

#' Print nlcom results
#'
#' @param x nlcom object
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
