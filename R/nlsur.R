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

#' Non-Linear Seemingly Unrelated Regression
#'
#' \code{.nlsur()} is a function for estimation of a non-linear seemingly
#'  unrelated regression model in R.
#'
#' @param eqns is can be a single equation or a quation system. If eqns is a
#' single equation it will internaly be converted to a list. Estimation of a
#' single equation might as well be done using \code{nls()}.
#' @param data is the data set on which the equation is applied. This can be of
#' every type \code{eval()} can handle.
#' @param startvalues is a vector of intial start values. For
#' @param S is a weighing matrix used for estimation in Feasible Generalized
#' Non-Linear Least Squares (FGNLS) and Iterative FGNLS. For \code{nlsur()}
#' this is assumed to be the identity matrix. Hence, it is not included. If
#' included S is expected to be a matrix.
#' @param nls is a logical and default if estimation is done for NLSUR or NLS.
#' @param fgnls is a logical and must be set, if estimation is done for FGNLS.
#' This is called in a function called \code{fgnls()} and should not be set by
#' the user.
#' @param ifgnls is a logical and must be set, if estimation is done for ifgnls.
#' This is called in a function called \code{nlsur()} and should not be set by
#' the user.
#' @param qrsolve is a logica, if TRUE \code{qr.coef(qr(x), r)} is called which
#' should be the most robust way for estimation of nls. For this all equations
#' will be rbinded, which might lead to memory bottlenecks.
#' @param MASS is a logical, if TRUE \code{lm.gls()} is called for estimation of
#' a linear regression with a weighting matrix. Otherwise Rs matrix functions
#' will be used. Estimation results with MASS might be more stable, but it is
#' unable to handle a \code{Matrix::Diagonal()}. This will be converted with
#' \code{as.matrix()} which can be quite RAM consuming.
#' @param trace is a logical. If TRUE the current iterations SSR is called.
#' @param eps the epislon used for convergence in nlsur(). Default is 1e-5.
#' @param tau is another convergence variable. Default is 1e-3.
#' @param tol qr.solves tolerance for detecting linear dependencies.
#'
#' @details nlsur is a function for estimation of a non-linear least squares
#' (NLS). In addition to \code{nls()} it is capable of estimation of system of
#' equations. This estimation is done in a non-linear seemingly unrelated
#' regression approach.
#'
#' @references
#' Bates, D. M. and Watts, D. G. (1988) Nonlinear Regression Analysis and Its
#'  Applications, Wiley
#' @references
#' Gallant, A. Ronald (1987): Nonlinear Statistical Models. Wiley: New York
#' @importFrom Matrix diag kronecker rankMatrix
#' @importFrom parallel mcmapply mclapply
#' @importFrom stats as.formula coef deriv
#' @useDynLib nlsur, .registration=TRUE
#' @export .nlsur
.nlsur <- function(eqns, data, startvalues, S = NULL, robust = robust,
                   nls = FALSE, fgnls = FALSE, ifgnls = FALSE, qrsolve = FALSE,
                   MASS = FALSE, trace = FALSE, eps = eps, tau = tau,
                   maxiter = maxiter, tol = tol, initial = initial)
{
  z    <- list()
  itr  <- 0
  conv <- FALSE

  lhs  <- rhs <- ri <- xi <- list()
  r    <-   x <- NULL

  neqs <- length(eqns)
  n    <- vector("integer", length=neqs)
  k    <- vector("integer", length=neqs)

  wts  <- data$w

  nlsur_coef <- new.env(hash = TRUE)

  # set initial theta, if it contains NA values replace them with 0
  theta <- theta.old <- startvalues
  theta[is.na(theta)] <- 0

  # check if S-matrix was provided for estimation of nls step
  if (nls){
    if (is.null(S)) {
      if (trace)
        cat("create initial weight matrix Sigma.\n")
      S <- diag(1, ncol=neqs, nrow=neqs)
      nls <- TRUE # keep nls flag.
    } else {
      if (trace)
        cat("Use diagonal of Sigma matrix.\n")
      S <- diag(diag(S), nrow = nrow(S), ncol = ncol(S))
      nls <- FALSE # reset nls flag. needed to enter weighted least squares
    }
  }

  qS <- qr.solve(qr(S, tol = tol), tol = tol)
  s  <- chol(qS)


  eqns_lhs <- mclapply(X = eqns, FUN = function(x)x[[2L]])
  eqns_rhs <- mclapply(X = eqns, FUN = function(x)x[[3L]])

  eqnames <- sapply(X = eqns_lhs, FUN = function(x)capture.output(print(x)))

  ## assign theta: make them available for eval
  for (i in 1:length(theta)) {
    name <- names(theta)[i]
    val <- theta[i]
    storage.mode(val) <- "double"
    assign(name, val, envir = nlsur_coef)
  }

  #### Initial evaluation ------------------------------------------------------

  # begin equation loop: for (i in 1:neqs) {}
  lhs <- mclapply(X = eqns_lhs, FUN = eval, envir = data, enclos = nlsur_coef)
  rhs <- mclapply(X = eqns_rhs, FUN = eval, envir = data, enclos = nlsur_coef)
  ri  <- mcmapply("-", lhs, rhs, SIMPLIFY = FALSE)

  rm(lhs, rhs)

  x   <- mclapply(X = eqns, FUN = function(x) {
    attr(eval(deriv(x, names(theta)),
              envir = data, enclos = nlsur_coef), "gradient")
  })
  # end equation loop

  n <- as.integer(lapply(X = x, FUN = nrow))

  # rankMatrix uses svd and is slow use once only
  if (initial) {
    k <- as.integer(lapply(X = x, FUN = rankMatrix))

    z$k            <- k
    z$df           <- n - k
  }

  r <- do.call(cbind, ri)

  if (qrsolve | MASS) {
    x <- do.call(rbind, x)
  } else {
    x <- do.call(cbind, x)
  }

  # eval might return NaN
  if ( any (is.nan(x))){
    stop("NA/NaN/Inf in derivation found. Most likely due to artificial data.")
  }


  # Evaluate initial ssr
  ssr.old <- calc_ssr(r, s, wts)

  # Print initial ssr
  if (trace)
    cat("Initial SSR: ", ssr.old, "\n")



  # values for initial alpha and its divisor
  alph       <- 1
  divi       <- 2
  itr        <- 0
  alpha      <- alph # stepsizeparameter

  # repeat until convergence is reached
  while (!conv) {

    # convergence was not reached given maxiter
    if (itr == maxiter){
      message(paste(itr, "nls iterations and convergence not reached."),
              paste("Last theta is: \n", theta, "\n"))
      return(0)
    }

    # If alpha < alph increase it again. Spotted in nls.c
    alpha <- min(divi*alpha, alph)
    # Alt: Stata variant, set alpha to alph
    # alpha <- alph

    # initiate while loop
    ssr <- Inf
    theta.old <- theta

    # begin regression
    # Regression of residuals on derivs
    if (nls & qrsolve) {

      r <- matrix(r, ncol = 1)

      theta.new <- qr.coef(qr(x), r)

    } else {

      # MASS for calculation of wls
      if (MASS)
      {

        r <- matrix(r, ncol = 1)

        # Use MASS::lm.gls inspired function for the weighted regression, this
        # will call eigen() to pre-weight r and x which then can be solved with
        # qr. This will be slower than blockwise wls but is numerically stable
        theta.new <- lm_gls(X = x, Y = r, W = S, wts = wts, neqs = neqs, tol = tol)

      } else {

        # Weighted regression of residuals on derivs ---
        theta.new <- calc_reg(x, r, qS, wts, length(theta), 1, tol)

      }

    }
    # end regression

    # update theta
    theta.new        <- as.vector(theta.new)
    names(theta.new) <- names(theta)
    theta            <- theta.new

    while ( ssr > ssr.old )
    { # begin iter

      # use the scalar to get a new theta
      theta.new <- theta.old + alpha * theta

      # theta.new can be NA. Change NA to 0.
      coef_na <- names(theta.new)[is.na(theta.new)]
      theta.new[is.na(theta.new)] <- 0

      ## assign new thetas thetas = makes them available to eval
      for (i in 1:length(theta.new)) {
        name <- names(theta.new)[i]
        val <- theta.new[i]
        storage.mode(val) <- "double"
        assign(name, val, envir = nlsur_coef)
      }

      # eval eqn with the new theta
      lhs <- rhs <- ri <- xi <- list()
      r <- x <- NULL

      # begin equation loop: for (i in 1:neqs) {}
      lhs <- mclapply(X = eqns_lhs, FUN = eval, envir = data, enclos = nlsur_coef)
      rhs <- mclapply(X = eqns_rhs, FUN = eval, envir = data, enclos = nlsur_coef)
      ri  <- mcmapply("-", lhs, rhs, SIMPLIFY = FALSE)

      rm(lhs, rhs)

      x   <- mclapply(X = eqns, FUN = function(x) {
        attr(eval(deriv(x, names(theta)),
                  envir = data, enclos = nlsur_coef), "gradient")
      })
      # end equation loop

      r <- do.call(cbind, ri)

      if (qrsolve | MASS) {
        x <- do.call(rbind, x)
      } else {
        x <- do.call(cbind, x)
      }


      theta_na <- names(theta)[is.na(theta)]
      # x[,colnames(x) %in% theta_na] <- NA

      # Reevaluation of ssr
      ssr <- calc_ssr(r, s, wts)

      # divide stepsizeparameter
      alpha <- alpha/divi

    } # end iter

    ssr.old <- ssr

    # Print updated ssr
    if (trace)
      cat("SSR: ", ssr, "\n")

    # Stopping rules. [Gallant (1987) p.29]

    # ssr: |ssr.old - ssr| < eps | ssr.old + tau|

    conv1 <- !isTRUE(abs(ssr.old - ssr) >
                       eps * (ssr.old + tau))

    # theta: ||theta - theta.new|| < eps (||theta|| + tau)
    # conv2 <- !isTRUE(norm(as.matrix(theta - theta.new)) <
    #                    eps * (norm(as.matrix(theta)) + tau))

    # alpha was already divided
    conv2 <- !all( (alpha*divi) * abs(theta) >
                     eps * (abs(theta.old) + tau), na.rm = TRUE )

    # this is what Stata documents what they do for nl. include alpha?
    # conv2 <- all( alpha * abs(theta.new) <= eps * (abs(theta) + tau) )

    # both convergence criteria should be TRUE [Himmelblau (1972)] according to
    # Bates and Watts (1988) p.49
    conv <- all(conv1, conv2)

    itr <- itr + 1
    theta <- theta.new
    ssr.old <- ssr
    theta.old <- theta

  }

  # replace 0 values with NA
  theta[theta_na] <- NA

  ## Create covb matrix ##

  if (qrsolve | MASS){

    r <- do.call(cbind, ri)
    r <- matrix(r, ncol = 1)

    covb <- lm_gls(X = x, Y = r, W = S, wts = wts, neqs = neqs, tol = tol, covb = TRUE)
  } else {
    # get xdx from calc_reg
    covb <- calc_reg(x, r, qS, wts, length(theta), 0, tol)
  }

  # # if singularities are detected covb will contain cols and rows with NA
  covb <- covb[!is.na(theta), !is.na(theta)]
  covb <- qr.solve(qr(covb, tol = tol), tol = tol)

  # calculate robust standard errors
  if (robust) {
    # get xdu from calc_robust
    xdu <- calc_robust(x, r, qS, wts, length(theta))

    # only good rows
    xdu <- xdu[!is.na(theta), !is.na(theta)]

    # calc robust covb
    covb <- covb %*% xdu %*% covb
  }

  coef_names <- names(theta)[!(names(theta) %in% coef_na)]
  dimnames(covb) <- list(coef_names, coef_names)

  r <- do.call(cbind, ri)

  z$coefficients <- theta
  z$residuals    <- r
  z$eqnames      <- eqnames
  z$sigma        <- 1/n * crossprod(r, wts * r)

  z$n            <- n
  z$deviance     <- as.numeric(ssr)

  z$weights      <- wts
  z$cov          <- covb

  class(z) <- "nlsur"

  z

}

#' Fitting Iterative Feasible Non-Linear Seemingly Unrelated Regression Model
#'
#' \code{nlsur()} is used to fit nonlinear regression models. It can handle the
#' feasible and iterative feasible variants.
#'
#' @param eqns is a list object containing the model as formula. This list can
#' handle contain only a single equations (although in this case nls() might be
#' a better coice) or a system of equations.
#' @param startvalues initial values for the parameters to be estimated.
#' @param data an (optional) data frame containing the variables that will be
#' evaluated in the formula.
#' @param type can be 1 Nonlinear Least Squares (NLS), 2 Feasible Generalised
#' NLS (FGNLS) or 3 Iterative FGNLS (IFGNLS) or the respective abbrevations in
#' character form.
#' @param tol qr.solves tolerance for detecting linear dependencies.
#' @param eps the epislon used for convergence in nlsur(). Default is 1e-5.
#' @param tau is another convergence variable. Default is 1e-3.
#' @param ifgnlseps is epislon for ifgnls(). Default is 1e-10.
#' @param stata is a logical. If TRUE for nls a second evaluation will be run.
#' Stata does this by default. For this second run Stata replaces the diagonal
#' of the I matrix with the coefficients.
#' @param trace logical wheather or not SSR information should be printed.
#' Default is FALSE.
#' @param robust logical if true robust standard errors are estimated.
#' @param S is a weight matrix used for evaluation. If no weight matrix is
#' provided the identity matrix I will be used.
#' @param qrsolve logical
#' @param MASS is a logical wheather the MASS::lm.gls() function should be used
#' for weighted Regression. This can cause sever RAM usage as the weight matrix
#' tend to be huge (n-equations * n-rows).
#' @param weights Additional weight vector.
#' @param maxiter Maximum number of iterations.
#' @param val If no start values supplied, create them with this start value.
#' Default is 0.
#' @param initial logical value to define if rankMatrix is calculated every
#' iteration of nlsur.
#' @param multicores number of cores used for parallel mcl-/mcmapply if no value
#' is set this defaults to n-1.
#'
#' @details nlsur() is a wrapper around .nlsur(). The function was initialy
#' inspired by the Stata Corp Function nlsur.
#' Nlsur estimates a nonlinear least squares demand system. With nls, fgnls or
#' ifgnls which is equivalent to Maximum Likelihood estimation.
#' Nonlinear least squares requires start values and nlsur requires a weighting
#' matrix for the demand system. If no weight matrix is provided, nlsur will use
#' the identity matrix I. If type = 1 or type = "nls" is added, nlsur will use
#' the matrix for an initial estimation, once the estimation is done, it will
#' swap the diagonal with the estimated results.
#'
#' Most robust regression estimates shall be returned with both qrsolve and MASS
#' TRUE, but memory consumtion is largest this way. If MASS is FALSE a memory
#' efficent RcppArmadillo solution is used for fgnls and ifgnls. If qrsolve is
#' FALSE as well, only the Armadillo function is used.
#'
#' If \code{robust} is selected Whites HC0 is used to caclulate
#' Heteroscedasticity Robust Standard Errors.
#'
#' If \code{initial} is TRUE rankMatrix will be calculated every iteration of
#' nlsur. Meaning for nls at least once, for fgnls at least twice and for ifgnls
#' at least three times. This adds a lot of overhead, since rankMatrix is used
#' to calculate k. To assure that k does not change this can be set to TRUE.
#'
#' Nlsur has methods for the generic functions \link{coef}, \link{confint},
#' \link{deviance}, \link{df.residual}, \link{fitted}, \link{predict},
#' \link{print}, \link{residuals}, \link{summary} and \link{vcov}.
#' @return The function returns a list object of class nlsur. The list includes:
#' \describe{
#'   \item{coefficients:}{estimated coefficients}
#'   \item{residuals:}{residuals}
#'   \item{xi:}{residuals of each equation in a single list}
#'   \item{eqnames:}{list of equation names}
#'   \item{sigma:}{the weight matrix}
#'   \item{ssr:}{Residual sum of squares}
#'   \item{lhs:}{Left hand side of the evaluated model}
#'   \item{rhs:}{Right hand side of the evaluated model}
#'   \item{nlsur:}{model type. "NLS", "FGNLS" or "IFGNLS"}
#'   \item{se:}{standard errors}
#'   \item{t:}{t values}
#'   \item{covb:}{asymptotic covarince matrix}
#'   \item{zi:}{equation wise estimation results of SSR, MSE, RMSE, MAE, R2 and
#' Adj-R2. As well as n, k and df.}
#'   \item{model:}{equvation or system of equations as list containing
#' formulas}
#' }
#'
#' @examples
#' \dontrun{
#' # Greene Example 10.3
#' library(nlsur)
#'
#' url <- "http://www.stern.nyu.edu/~wgreene/Text/Edition7/TableF10-2.txt"
#' dd <- read.table(url, header = T)
#'
#' names(dd) <- c("Year", "Cost", "Sk", "Sl", "Se", "Sm", "Pk", "Pl", "Pe", "Pm")
#'
#'
#' eqns <- list( Sk ~ bk + dkk * log(Pk/Pm) + dkl * log(Pl/Pm) + dke * log(Pe/Pm),
#'               Sl ~ bl + dkl * log(Pk/Pm) + dll * log(Pl/Pm) + dle * log(Pe/Pm),
#'               Se ~ be + dke * log(Pk/Pm) + dle * log(Pl/Pm) + dee * log(Pe/Pm))
#'
#' strtvls <- c(be = 0, bk = 0, bl = 0,
#'              dkk = 0, dkl = 0, dke = 0,
#'              dll = 0, dle = 0, dee = 0)
#'
#'
#' erg <- nlsur(eqns = eqns, data = dd, startvalues = strtvls, type = 2,
#'              trace = TRUE, eps = 1e-10)
#'
#' erg
#' }
#' @references Gallant, A. Ronald (1987): Nonlinear Statistical Models.
#'  Wiley: New York
#' @seealso \link{nls}
#' @importFrom parallel mclapply detectCores
#' @importFrom stats as.formula coef na.omit
#' @importFrom utils capture.output
#' @useDynLib nlsur
#'
#' @export
nlsur <- function(eqns, data, startvalues, type=NULL, S = NULL,
                  trace = FALSE, robust = FALSE, stata = TRUE, qrsolve = FALSE,
                  weights, MASS = FALSE, maxiter = 1000, val = 0,
                  tol = 1e-7, eps = 1e-5, ifgnlseps = 1e-10,
                  tau = 1e-3, initial = FALSE, multicores) {

  mc <- getOption("mc.cores")

  # set multicore process
  if (missing(multicores))
    multicores <- detectCores()-1

  options("mc.cores" = multicores )

  # Check if eqns might be a formula
  if (!is.list(eqns)){
    # check if formula provided as string
    if (!is.formula(eqns)){
      gl <- as.formula(eqns)

      if (!is.formula(gl))
        stop("No equation specified.")

    } else
      gl <- eqns

    eqns <- list()
    eqns[[1]] <- gl
  }

  # Check if original call contains weights
  if (missing(weights))
    wts <- NULL
  else
    wts <- as.character(substitute(weights))

  # If no startvalues supplied, create them.
  if (missing (startvalues)) {
    msg <- paste("startvalues created with val =", val)
    message(msg)

    startvalues <- getstartvals(model = eqns, data = data, val = val)
  }

  # Check if all variables that are not startvalues exist in data.
  vars <- unlist(mclapply(eqns, all.vars))
  vars <- vars[which(!vars %in% names(startvalues))]
  ok   <- all(vars%in%names(data))

  # if not ok bail out
  if (!ok) {
    msg <- paste("Missmatch in model and dataset. Some equation variables \n",
                 "are not in startvalues nor data object.")
    message(msg)
    return(0)
  }

  # lm.gls does not allow weights
  if ((isTRUE(qrsolve)) & !is.null(wts))
    stop("With qrsolve and MASS you can not use weights.")


  # remove observation, if observation a parameter contains NA.
  modelparameters <- c(unlist(mclapply(eqns, all.vars)), wts)
  parms <- modelparameters[which(!modelparameters %in% names(startvalues))]

  # check for equation constants
  eqconst <- lapply(X = eqns, FUN = constant)

  # Check for wts
  if ( is.null(wts) ) {
    data$w <- 1
  } else {
    wts <- as.name(wts)

    data$w <- eval(substitute(wts), data)
  }

  # include weights to assure the correct length
  # of weights if missings are excluded.
  data <- na.omit(data[unique(c(parms,"w"))])

  nls  <- fgnls <- ifgnls <- FALSE
  n    <- nrow(data)
  z    <- NULL

  # normwts
  data$w <- data$w/sum(data$w) * n

  cl <- match.call()

  # define what will be done
  if (!is.null(type)) {
    if(type == "NLS" | type == 1) {
      nls <- TRUE
    } else {
      if (type == "FGNLS" | type == 2) {
        fgnls <- TRUE
      } else {
        if (type == "IFGNLS" | type == 3) {
          fgnls  <- TRUE
          ifgnls <- TRUE
        }
      }
    }
  }

  # Estimation of NLS ##########################################################

  # print NLS
  if (trace)
    cat("-- NLS\n")

  z <- .nlsur( eqns = eqns, data = data, startvalues = startvalues, S = S,
               robust = robust, nls = TRUE, trace = trace, qrsolve = qrsolve,
               MASS = MASS, eps = eps, tau = tau, maxiter = maxiter,
               tol = tol, initial = TRUE)

  # evaluatated at initial stage
  n <- z$n; k <- z$k; df <- z$df

  # To update standard errors in nls case Stata estimates a second nls with
  # diag(S) instead of I.
  if (nls & stata) {

    # backup of theta and S, remove z
    theta.old <- coef(z); S <- z$sigma; rm(z)

    z <- .nlsur( eqns = eqns, data = data, startvalues = theta.old, S = S,
                 robust = robust, nls = nls, trace = trace, qrsolve = FALSE,
                 MASS = MASS, eps = eps, tau = tau, maxiter = maxiter,
                 tol = tol, initial = initial)

    # Stata uses the orignal sigma for covb
    z$sigma <- diag(diag(S), nrow = n, ncol = n)

  }
  z$nlsur <- "NLS"


  # Estimation of FGNLS ########################################################
  if (fgnls) {

    # print FGNLS
    if (trace)
      cat("-- FGNLS\n")

    # backup of theta and S, remove z
    theta.old <- coef(z); S <- z$sigma; rm(z)

    z <- .nlsur(eqns = eqns, data = data, startvalues = theta.old, S = S,
                robust = robust, nls = FALSE, trace = trace, qrsolve = FALSE,
                MASS = MASS, eps = eps, tau = tau, maxiter = maxiter,
                tol = tol, initial = initial)

    # Stata uses the orignal sigma for covb
    if (!ifgnls)
      z$sigma <- S

    z$nlsur <- "FGNLS"

    # Estimation of IFGNLS #####################################################
    if (ifgnls) {

      # print IFGNLS
      if (trace)
        cat("-- IFGNLS\n")

      S <- z$sigma
      conv <- FALSE
      iter <- 0

      # repeat nlsur estimation until convergence is reached
      while (!conv)
      {

        # if convergence is not reached
        if (iter == maxiter){

          msg <- paste(iter, "nls iterations and convergence not reached.\n",
                       "Last theta is: \n")

          message(msg)
          print(coef(z))
          return(0)
        }

        # backup of theta and S remove z
        theta.old <- coef(z); S.old <- S; rm(z)

        z <- .nlsur(eqns = eqns, data = data, startvalues = theta.old,
                    S = S, robust = robust, nls = FALSE, qrsolve = FALSE,
                    MASS = MASS, eps = eps, tau = tau, maxiter = maxiter,
                    tol = tol, initial = initial)

        S     <- z$sigma
        r     <- z$residuals
        theta <- coef(z)

        s   <- chol(qr.solve(S, tol = tol))
        rss <- calc_ssr(r, s, data$w)

        iter <- iter +1

        maxthetachange <- max(abs(theta.old - theta) /
                                ( abs(theta) +1),
                              na.rm = TRUE )
        maxSigmachange <- max(abs(S.old - S) /
                                (abs(S.old) + 1),
                              na.rm = TRUE)

        # conv1 <- abs(rss.old - rss) < eps * (rss.old + tau)
        # conv2 <- norm(as.matrix(z.old$coefficients - z$coefficients)) <
        #   eps * (norm(as.matrix(z.old$coefficients)) + tau)
        # conv <- any(conv1, conv2)

        conv1 <- maxthetachange < eps
        conv2 <- maxSigmachange < ifgnlseps

        conv <- all(conv1, conv2)

        # Iteration output
        if (trace)
          cat("Iteration", iter, ": SSR", rss, "\n")

      }
      message <- paste("Convergence after iteration:", iter,".")


      z$message   <- message
      z$nlsur     <- "IFGNLS"

    }
  }

  # initial is true so collect final k and df
  if(initial) {
    n <- z$n; k <- z$k; df <- z$df
  }

  # Estimate log likelihood ####################################################
  S <- z$sigma
  N <- unique(n)
  M <- nrow(S)

  LL <- ( sum(log(data$w)) -(M*N) * (log(2 * pi) +
                                       1 - log(N) +
                                       log(det(S)) / M  +
                                       log(sum(data$w))) )/2


  # Fitted values ##############################################################
  nlsur_coef <- new.env(hash = TRUE)
  eqns_lhs   <- mclapply(X = eqns, FUN = function(x)x[[2L]])
  eqns_rhs   <- mclapply(X = eqns, FUN = function(x)x[[3L]])
  eqnames    <- sapply(X = eqns_lhs, FUN = function(x)capture.output(print(x)))
  theta      <- coef(z)

  # assign theta: make them available for eval
  for (i in 1:length(theta)) {
    name <- names(theta)[i]
    val <- theta[i]
    storage.mode(val) <- "double"
    assign(name, val, envir = nlsur_coef)
  }

  fitted <- mclapply(X = eqns_rhs, FUN = eval, envir = data, enclos = nlsur_coef)
  fitted <- as.data.frame(fitted)
  names(fitted) <- eqnames

  # create output
  z$n           <- n
  z$k           <- k
  z$df          <- df
  z$df.residual <- df
  z$LL          <- LL
  z$model       <- eqns
  z$const       <- eqconst
  z$data        <- data
  z$call        <- cl
  z$start       <- startvalues
  z$nlsonly     <- all(nls & !stata)
  z$robust      <- robust
  z$fitted      <- fitted

  # if call did not contain weights: drop them
  if (is.null(wts))
    z$weights <- NULL


  options("mc.cores" = mc )

  z
}

#' @export
print.nlsur <- function(x, ...) {
  # ... is to please check()
  print(x$coefficients, ...)
}

#' Summary of nlsur objects
#'
#' @param object object of class nlsur.
#' @param noconst logical value determining if model uses a constant or not.
#' @param multicores number of cores to be used, default n - 1
#' @param ... additional parameters (currently not used)
#' @importFrom parallel mclapply
#' @importFrom stats as.formula pt residuals weights
#' @export
summary.nlsur <- function(object, noconst = TRUE, multicores, ...) {
  # ... is to please check()

  z <- object


  mc <- getOption("mc.cores")

  # set multicore process
  if (missing(multicores))
    multicores <- detectCores()-1

  options("mc.cores" = multicores )

  data    <- z$data
  eqns    <- z$model
  neqs    <- length(eqns)
  w       <- weights(z)
  n       <- z$n
  k       <- z$k
  df      <- z$df
  r       <- residuals(z)
  eqconst <- z$const
  nlsonly <- z$nlsonly
  est     <- z$coefficients

  hasconst <- sapply(eqconst, is.character)

  # check weights
  if (is.null(w)) {
    w <- rep(1, nrow(data))
  } else {
    # check maleformed weights
    if (!all(w > 0))
      stop("Negative or zero weight found.")
  }

  eqns_lhs <- mclapply(X = eqns, FUN = function(x)x[[2L]])


  #### Estimation of covariance matrix, standard errors and z/t-values ####
  lhs     <- list()
  scale   <- vector("numeric", length=neqs)      # scalefactor
  div     <- vector("numeric", length=neqs)      # divisor
  wi      <- vector("numeric", length=neqs)      # normalized wts
  ssr     <- vector("numeric", length=neqs)      # sum of squared residuals
  mss     <- vector("numeric", length=neqs)
  mse     <- vector("numeric", length=neqs)      # mean square error
  rmse    <- vector("numeric", length=neqs)      # root of mse
  mae     <- vector("numeric", length=neqs)      # mean absolute error
  r2      <- vector("numeric", length=neqs)      # R-squared value
  adjr2   <- vector("numeric", length=neqs)      # adjusted R-squared value


  nlsur_coef <- new.env(hash = TRUE)

  # Assign values for eval
  for (i in 1:length(est)) {
    name <- names(est)[i]
    val <- est[i]
    storage.mode(val) <- "double"
    assign(name, val, envir = nlsur_coef)
  }

  lhs   <- mclapply(X = eqns_lhs, FUN = eval, envir = data, enclos = nlsur_coef)
  scale <- n/sum(w)
  div   <- n -1



  # Evaluate everything required for summary printing
  for (i in 1:neqs) {

    # if lhs is a constant eg 0, size of lhs_i and w differs
    lhs_i <- lhs[[i]]

    if (length(lhs_i) < NROW(data))
      lhs_i <- rep(lhs_i, NROW(data))

    ssr[i]   <- sum( r[,i]^2 * w) * scale[i]

    # No constant found
    if (hasconst[i]) {
      wi       <- w/sum(w) * n[i]

      lhs_wm   <- wt_mean(x = lhs_i, w = wi)
      wvar     <- (1/div[i]) * sum( wi * (lhs_i - lhs_wm)^2)
      mss[i]   <- wvar * div[i] - ssr[i]
    } else{
      mss[i]   <- sum(w * lhs_i^2) * scale[i] - ssr[i]
    }

    mse[i]   <- ssr[i] / n[i]
    rmse[i]  <- sqrt(mse[i])
    mae[i]   <- sum(abs(r[, i]))/n[i]

    r2[i]  <- mss[i] / (mss[i] + ssr[i])

    # correct adjr2 if n - df = 1
    if (div[i] > 1) {
      adjr2[i] <- 1 - (div[i] / df[i]) * (1 - r2[i])
    } else {
      adjr2[i] <- 1 - (n[i] / df[i]) * (1 - r2[i])
    }
  }

  nE <- sum(n) / sum ( w/sum(w) )
  kE  <- sum(k)


  resvar <- 1
  # for single equation
  if (neqs == 1 & nlsonly)
    resvar <- ssr / df

  # Estimate se, tval and prob
  se <- sqrt(diag(z$cov) * resvar)

  # Add names of vars that lead to these se values
  names(se) <- names(est[!is.na(est)])

  # Add names of vars that have no se value
  se_na <- est[is.na(est)]
  names(se_na) <- names(est[is.na(est)])

  # Combine and change order
  se <- c(se, se_na)
  se <- se[names(est)]

  tval <- est / se

  # z vs t
  prob <- 2 * (1 - pt(abs(tval), (nE * kE )))
  # if single equation
  if (neqs == 1 & nlsonly)
    prob <- 2 * pt(abs(tval), df, lower.tail = FALSE)

  # per equation statistics
  zi <- cbind(n, k, rmse, mae, r2, adjr2)
  zi <- as.data.frame(zi)

  # replace all character(0) if a equation does not contain a constant
  for (i in 1:neqs){
    eqconst[[i]][identical(eqconst[[i]], character(0))] <- ""
  }

  # add constant variables to summary
  if (any(hasconst)){
    eqconst <- do.call(rbind, eqconst)
    neqconst <- ncol(eqconst)

    zi <- data.frame(cbind(zi, eqconst))
  }

  zi <- data.frame(as.character(z$eqnames),zi)

  cnst <- character(0)
  # if a equation contians more than one const only add it once and fill the
  # rest with blanks
  if (any(hasconst)) {
    cnst <- c("Const")
    if (neqconst>1){
      cnst <- c(cnst, rep(x = "", (neqconst-1)))
    }
  }

  colnames(zi) <- c("","n", "k", "RMSE", "MAE", "R-squared",
                    "Adj-R-sqr.", cnst)

  # ans: returned object
  ans <- NULL

  # ans$coefficients <- z[c("coefficients", "se", "t")]
  ans$coefficients <- cbind(est, se, tval, prob)

  cnames <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  # for single equation return t values
  if (neqs == 1 & nlsonly)
    cnames <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")

  dimnames(ans$coefficients) <- list(
    names(z$coefficients),
    cnames
  )

  ans$residuals    <- residuals
  ans$nlsur        <- z$nlsur
  ans$zi           <- zi
  ans$weights      <- weights(z)
  ans$cov          <- z$cov

  # for ifgnls add log likelihood
  if (ans$nlsur == "IFGNLS")
    ans$LL <- z$LL


  options("mc.cores" = mc )

  class(ans) <- "summary.nlsur"

  ans
}

#' @importFrom stats printCoefmat weights
#' @export
print.summary.nlsur <- function(x, digits, ...) {

  # ... is to please check()
  cat("NLSUR Object of type:", x$nlsur, "\n\n")

  # check if estimation contains weights
  if (!is.null(weights(x))) {
    cat("Scaled R-squared: \n\n")
  }

  # digits to be presented
  if(missing(digits))
    digits <- 4

  print(x$zi, digits = digits)

  cat("\n")
  cat("Coefficients:\n")

  # Weights again
  if (!is.null(weights(x))) {
    cat("Weighted nlsur: \n\n")
  }

  printCoefmat(x$coefficients, digits = digits, ...)

  # For ifgnls print log likelihood
  if (x$nlsur == "IFGNLS")
    cat("Log-Likelihood:", x$LL, "\n")
}

#' @export
logLik.nlsur <- function(object, ...) {
  z <- object$LL
  class(z) <- "logLik"

  z
}

#' @export
vcov.nlsur <- function(object, ...) {
  object$cov
}

#' @export
vcov.summary.nlsur <- function(object, ...) {
  object$cov
}

#' Predict for Non-Linear Seemingly Unrelated Regression Models
#'
#' \code{predict()} is a function to predict nlsur results.
#'
#' @param object is an nlsur estimation result.
#' @param newdata an optional data frame for which the prediction is evaluated.
#' @param ... further arguments for predict. At present no optional arguments
#'  are used.
#'
#' @details predict.nlsur evaluates the nlsur equation(s) given nlsurs
#' estimated parameters using either the original data.frame or newdata. Since
#' nlsur() restricts the data object only to complete cases observations with
#' missings will not be fitted.
#'
#'
#' @examples # predict(nlsurObj, dataframe)
#' @importFrom stats coef
#' @importFrom parallel mclapply
#'
#' @export
predict.nlsur <- function(object, newdata, ...) {

  eqs <- object$model

  eqns_lhs <- mclapply(X = eqs, FUN = function(x)x[[2L]])
  eqns_rhs <- mclapply(X = eqs, FUN = function(x)x[[3L]])
  vnam     <- sapply(X = eqns_lhs, FUN = as.character)

  # check for newdata
  if (missing(newdata) || is.null(newdata)) {
    data <- object$data
  } else {
    data <- newdata
  }

  data2 <- data.frame(data, as.list(coef(object)))

  # create fit: predict result
  fit <- mclapply(X = eqns_rhs, FUN = eval, envir = data2)
  fit <- data.frame(fit)
  names(fit) <- vnam

  fit
}

#' Calculate WLS using sparse matrix and qr
#'
#' @param X n x m X matrix
#' @param Y n x k matrix
#' @param W n x n
#' @param wts vector with weights
#' @param neqs k
#' @param tol tolerance for qr
#' @param covb if true covb is calculated else theta
#' @importFrom Matrix crossprod kronecker Diagonal
#' @importFrom methods as
#' @export
lm_gls <- function(X, Y, W, wts, neqs, tol = 1e-7, covb = FALSE) {

  # if (!missing(wts) & !all(wts == 1)) {
  #   X <- wts * X
  #   Y <- wts * Y
  # }


  eW <- eigen(W, TRUE)
  d <- eW$values
  if (any(d <= 0))
    stop("'W' is not positive definite")

  A <- diag(d^-0.5,
            nrow = length(d),
            ncol = length(d)) %*% t(eW$vectors)

  n <- nrow(X)/neqs

  A <- Matrix::kronecker(X = A,
                         Y = Matrix::Diagonal(n))

  X <- as(X, "sparseMatrix"); Y <- as(Y, "sparseMatrix")

  if (covb) {
    fit <- Matrix::crossprod(A %*% X)
  } else {
    fit <- qr.coef(qr(wts * (A %*% X), tol = tol), wts * (A %*% Y) )
  }

  fit
}
