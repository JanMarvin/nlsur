#' Estimate nonlinear combinations of nlsur estimates
#'
#' @param object of class nlsur
#' @param form formula e.g. "be/bk". May contain names(coef(object)) or prior
#' nlcom estimations.
#' @param alpha value for conf. interval default is 0.05
#' @param rname optional rowname for result
#' @param envir optional name of environment to search for additional parameters
#' @importFrom stats coef formula pt qnorm
#' @seealso \link[car]{deltaMethod}
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
nlcom <- function(object, form, alpha = 0.05, rname, envir) {


  if (!inherits(object, "nlsur"))
    stop("no nlsur object")

  # store original form
  oform <- form

  # check if all form vars are part of object
  vars <- all.vars(formula(paste(form, "~ 0")))
  vars <- vars[!(vars %in% names(coef(object)))]

  if (missing(envir))
    envir <- environment()

  # if form contains other nlcom estimates, fetch their formula and update
  # the current form with the older formula.
  if( !identical(vars, character(0)) ){
    for (i in vars){
      tmp <- get0(x = i, envir = envir)

      if(class(tmp) == "nlcom") {
        fname <- paste("(", attr(tmp, "form"), ")")

        form <- gsub(pattern = i, replacement = fname, x = form, fixed = TRUE)
      } else {
        if (class(tmp) == "numeric" || class(tmp) == "integer")
          form <- gsub(pattern = i, replacement = tmp, x = form, fixed = TRUE)
      }
    }
  }

  alevel <- 1 - alpha

  z <- dm(object, form, level = alevel)

  # separate the convinterval from the z output
  cinterv <- z[c("lwr", "upr")]

  # only keep est and se
  z <- z[c("est", "se")]

  # calculate t-value
  tval <- z$est / z$se

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

#' @method print nlcom
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

#' dm simple delta method implementation
#'
#' @param object of class nlsur
#' @param form formula e.g. "be/bk".
#' @param level value for conf. interval default is 0.05
#' @importFrom stats D vcov
dm <- function (object, form, level=0.05) {

  df   <- as.data.frame(t(coef(object)))
  nams <- names(df)
  fml  <- parse(text = form)
  cvs  <- seq_along(df)

  est <- eval(fml, envir = df)
  der <- sapply(cvs, FUN = function(x) eval(D(fml, nams[x]), envir = df))

  se <- as.vector(sqrt(
    t(der) %*% vcov(object) %*% der
  ))

  p <- (1 - level)/2
  z <- - qnorm(p)

  lwr <- est - z * se
  upr <- est + z * se

  z <- data.frame(est, se, lwr, upr, row.names = form)

  z
}
