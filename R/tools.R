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

#' Function to create startvalues for nlsur models
#' @param model nlsur model
#' @param data the data frame used for evaluation
#' @param val value
#'
#' @export
getstartvals <- function(model, data, val) {
  # automatic creation of start values
  modelparameters <- unlist(lapply(model, all.vars))
  svals <- modelparameters[which(!modelparameters %in% names(data))]
  svals <- unique(svals[order(svals)])
  # svals
  strtvls <- rep(val, length(svals))
  names(strtvls) <- svals

  strtvls
}


#' Check if object is of class formula
#' @param x object
#' @export
is.formula <- function(x) {
  isTRUE(class(x) == "formula")
}


#' Check if formula contains constant
#' @param x formula
#' @details Primitive function to check a formula for a constant part.
#' Function checks first and last term on rhs for constant variables
#' at front and back position.
#' @examples
#' \dontrun{
#' constant(y ~ x + a * z) # x
#' constant(y ~ x * b + 1) # 1
#' constant(y  ~ 0 + x) # NULL
#' constant(y  ~ x) # x
#' constant(y ~ x1 * b1 + b0 + x2 * b2) # wont find b0
#' constant( y  ~ (x*b +k) + a*y + b*z ) # wont find k
#' constant( y  ~ (k+ x*b) + a*y + b*z ) # k
#' constant( y  ~  a*y + b*z + (k + x*b) ) # wont find k
#' constant( y  ~  a*y + b*z + (x*b + k) ) # k
#' }
#' @export
constant <- function(x) {

  if(!is.formula(x))
    # try if call can be written as formula
    x <- as.formula(x)
  # try again
  if (!is.formula(x))
    stop("requires formula")

  rhs <- x[[3L]]

  # check front side of formula
  while ( !(is.symbol(rhs) | is.numeric(rhs)) ) {
    # go down one level and check if formula contains a new part
    # not multiplied containing a rhs
    if ( (rhs[[1]] == "+") | (rhs[[1]] == "-") | rhs[[1]] == "(") {
      rhs <- rhs[[2L]]
    } else {
      # reset rhs no constant found
      rhs <- x[[3L]]
      break
    }
  }

  # check back side of formula
  while ( !(is.symbol(rhs) | is.numeric(rhs)) ) {
    # go down one level and check if formula contains a constant
    if ( (rhs[[1]] == "+") | (rhs[[1]] == "-") ) {
      rhs <- rhs[[3L]]
    } else {
      if (rhs[[1]] == "(") {
        rhs <- rhs[[2L]]
      } else {
        break
      }
    }
  }

  z <- NULL
  # symbol as character
  if (is.symbol(rhs))
    z <- as.character(rhs)
  # numeric as numeric 0 is NULL
  if (is.numeric(rhs) & !identical(rhs, 0))
    z <- rhs

  z
}
