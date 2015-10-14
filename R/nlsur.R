# Copyright (c) 2015 Jan Marvin Garbuszus
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
#' \code{nlsur()} is a function for estimation of a non-linear seemingly
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
#' @param debug is a logical if debug output should be included. This can be a
#' lot and might not be of any help to you.
#' @param nls is a logical and default if estimation is done for NLSUR or NLS.
#' @param fgnls is a logical and must be set, if estimation is done for FGNLS.
#' This is called in a function called \code{fgnls()} and should not be set by
#' the user.
#' @param ifgnls is a logical and must be set, if estimation is done for ifgnls.
#' This is called in a function called \code{ifgnls()} and should not be set by
#' the user.
#' @param MASS is a logical, if TRUE \code{lm.gls()} is called for estimation of
#' a linear regression with a weighting matrix. Otherwise Rs matrix functions
#' will be used. Estimation results with MASS might be more stable, but it is
#' unable to handle a \code{Matrix::Diagonal()}. This will be converted with
#' \code{as.matrix()} which can be quite RAM consuming.
#' @param trace is a logical. If TRUE the current iterations SSR is called.
#' @param solvetol is the \code{tol} value of \code{qr.solve()}. For certain
#' estimations it is required to be fairly huge. Otherwise a false singularity
#' might be reported.
#'
#' @details nlsur is a function for estimation of a non-linear least squares
#' (NLS). In addition to \code{nls()} it is capable of estimation of system of
#' equations. This estimation is done in a non-linear seemingly unrelated
#' regression approach.
#'
#' @references
#' Bates, D. M. and Watts, D. G. (1988) Nonlinear Regression Analysis and Its
#'  Applications, Wiley
#'  Gallant, A. Ronald (1987): Nonlinear Statistical Models. Wiley: New York
#' @useDynLib nlsur
#' @import RcppArmadillo

#' @export
nlsur <- function(eqns, data, startvalues, S = NULL, debug = FALSE,
                  nls = FALSE, fgnls = FALSE, ifgnls = FALSE,
                  MASS = FALSE, trace = FALSE,
                  solvetol = solvetol, eps = eps, tau = tau)
{
  z    <- list()
  itr  <- 0
  conv <- FALSE

  lhs  <- rhs <- ri <- xi <- list()
  r    <-   x <- NULL

  n    <- nrow(data)
  neqs <- length(eqns)

  # set initial theta
  theta <- startvalues

  if (nls){
    if (is.null(S)) {
      if (trace)
        cat("create initial weight matrix Sigma.\n")
      S <- diag(1, ncol=neqs, nrow=neqs)
      nls <- TRUE # keep nls flag
    } else {
      if (trace)
        cat("Use diagonal of Sigma matrix.\n")
      S <- diag(diag(S), nrow = nrow(S), ncol = ncol(S))
      nls <- FALSE # reset nls flag. needed to enter weighted least squares
    }
  }


  if (debug)
    print(S)

  qS <- qr.solve(S)
  s  <- chol(qS)

  eqnames <- NULL

  ## assign theta: make them available for eval
  for (i in 1:length(theta)) {
    name <- names(theta)[i]
    val <- theta[i]
    storage.mode(val) <- "double"
    assign(name, val)
  }

  #### Initial evaluation ------------------------------------------------------
  # Evaluate inital lhs, rhs, ri, r and xi and x
  for (i in 1:neqs) {
    eqnames <- c(eqnames, as.formula(eqns[[i]])[[2L]])
    lhs[[i]] <- eval(as.formula(eqns[[i]])[[2L]], envir = data)
    rhs[[i]] <- eval(as.formula(eqns[[i]])[[3L]], envir = data)

    ri[[i]] <- lhs[[i]] - rhs[[i]]

    xi[[i]] <- attr(with(data, with(as.list(theta),
                                    eval(deriv(eqns[[i]], names(theta)),
                                         envir = data))), "gradient")
  }

  r <- do.call(cbind, ri)
  x <- do.call(cbind, xi)

  if ( any (is.nan(x))){
    # eval might return NaN. lm.gls will complain, since r is smaller than S.
    # Fix this by changing its value to zero.
    stop("NA/NaN/Inf in derivation found. Most likely due to artificial data.")
  }

  if (debug)
    print(r)

  # Evaluate initial ssr
  ssr.old <- calc_ssr(r, s, neqs)

  if (trace)
    cat("Initial SSR: ", ssr.old, "\n")

  itr        <- 0
  alpha      <- 1 # stepsizeparameter

  while (!conv) {

    if (debug)
      cat("Iteration: ", itr , "\n")

    if (itr == 1000){
      message(paste(itr, "nls iterations and convergence not reached."),
              paste("Last theta is: \n", theta, "\n"))
      return(0)
    }

    # If alpha < 1 increase it again. Spotted in nls.c
    alpha <- min(2*alpha, 1)
    # Alt: Stata variant, set alpha to 1
    # alpha <- 1

    # initiate while loop
    ssr <- Inf
    theta.old <- theta

    # r <<- r
    # x <<- x
    # qS <<- qS

    # begin regression
    # Regression of residuals on derivs
    if (nls) {

      r <- matrix(r, ncol = 1)
      x <- do.call(rbind, xi)

      theta.new <- qr.coef(qr(x), r)

      if (any(is.na(theta.new))) {
        warning("fix NA in theta.new")
        message(
          "During nls for the following variables NA was replaced with 0."
        )
        print(names(theta)[is.na(theta.new)])
        theta.new[is.na(theta.new)] <- 0
      }

      theta.new <- as.vector(theta.new)
      names(theta.new) <- names(theta)
      theta <- theta.new

    } else {

      # r     <<- r
      # x     <<- x
      # qS    <<- qS
      # theta <<- theta
      # neqs  <<- neqs

      # Weighted regression of residuals on derivs ---
      theta_test <- calc_reg(x, r, qS, length(theta), neqs, 1)
      theta.new <- as.vector(theta_test)

      names(theta.new) <- names(theta)
      theta <- theta.new
    }
    # end regression

    if (debug)
      cat("enter while ( ssr > ssr.old ) loop\n")

    while ( ssr > ssr.old )
    { # begin iter

      if (debug)
        cat("alpha: ", alpha, "\n")

      # use the scalar to get a new theta
      theta.new <- startvalues + alpha * theta

      ## assign new thetas thetas = makes them available to eval
      for (i in 1:length(theta.new)) {
        name <- names(theta.new)[i]
        val <- theta.new[i]
        storage.mode(val) <- "double"
        assign(name, val)
      }

      # eval eqn with the new theta
      lhs <- rhs <- ri <- xi <- list()
      r <- x <- NULL

      for (i in 1:neqs) {
        lhs[[i]] <- eval(as.formula(eqns[[i]])[[2]], envir = data)
        rhs[[i]] <- eval(as.formula(eqns[[i]])[[3]], envir = data)
        ri[[i]] <- lhs[[i]] - rhs[[i]]

        xi[[i]] <- attr(with(data, with(as.list(theta.new),
                                        eval(deriv(eqns[[i]], names(theta.new)),
                                             envir = data))), "gradient")
      }

      r <- do.call(cbind, ri)
      x <- do.call(cbind, xi)

      # Evaluate initial ssr
      ssr <- calc_ssr(r, s, neqs)

      # divide stepsizeparameter
      alpha <- alpha/2


      if (debug)
        cat("SSR :", ssr, "SSR_Old:", ssr.old, "\n")

    } # end iter

    ssr.old <- ssr

    if (debug)
      print(r)

    if (trace)
      cat("SSR: ", ssr, "\n")

    if(debug){
      print(warnings())
      b <- cbind(theta.old, theta)
      print(b)
    }

    # Stopping rule. [Gallant (1987) p.29]
    # Note: R uses a different convergence criterium

    # ssr: |ssr.old - ssr| < eps | ssr.old + tau|
    # conv1 <- abs(ssr.old - ssr) < eps * (ssr.old + tau)

    # theta: ||theta - theta.new|| < eps (||theta|| + tau)
    # conv2 <- norm(as.matrix(theta - theta.new)) <
    #   eps * (norm(as.matrix(theta)) + tau)

    # no idea why, but Stata uses this
    conv1 <- !isTRUE(abs(ssr.old - ssr) > eps * (ssr.old + tau))
    # conv1 <- TRUE

    conv2 <- !isTRUE(all( alpha * abs(theta) > eps * (abs(startvalues) + tau) ))
    # conv2 <- !isTRUE( alpha * all(abs(theta - theta.new) > eps * (theta + tau)) )
    # print(theta)
    # conv2 <- FALSE

    # and this is what Stata documents what they do for nl
    # conv2 <- all( alpha * abs(theta.new) <= eps * (abs(theta) + tau) )

    # both convergence criteria must be TRUE
    conv <- all(conv1, conv2)

    if(debug)
      print(b)

    itr <- itr + 1
    theta <- theta.new
    ssr.old <- ssr
    startvalues <- theta


    if(debug)
      print(itr)

  }

  z$coefficients <- theta
  z$residuals <- r
  z$xi <- xi
  z$eqnames <- eqnames
  z$sigma <- 1/n * crossprod(r)
  z$ssr <- ssr

  z$lhs <- lhs
  z$rhs <- rhs

  attr(z, "class") <- "nlsur"

  z

}

#' @export
ifgnls <- function(eqns, data, startvalues, type=NULL, S = NULL, debug = FALSE,
                   trace = FALSE, solvetol = .Machine$double.eps, nls = nls,
                   MASS = FALSE, eps = 1e-5, ifgnlseps = 1e-10, tau = 1e-3) {

  # Check if all variables that are not startvalues exist in data.
  vars <- unlist(lapply(eqns, all.vars))
  vars <- vars[which(!vars %in% names(startvalues))]
  ok <- all(vars%in%names(data))
  if (!ok) {
    message("Missmatch in model and dataset.")
    return(0)
  }

  # remove observation, if observation a parameter contains NA.
  modelparameters <- unlist(lapply(eqns, all.vars))
  parms <- modelparameters[which(!modelparameters %in% names(startvalues))]

  data <- na.omit(data[parms])

  neqs   <- length(eqns)
  nls    <- FALSE
  fgnls  <- FALSE
  ifgnls <- FALSE
  z      <- NULL
  zi     <- NULL
  n <- nrow(data)

  #
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


  # Estimation of NLS
  # nls
  if (trace)
    cat("-- NLS\n")

  z <- nlsur( eqns = eqns, data = data, startvalues = startvalues, S = S,
              debug = debug, nls = TRUE, trace = trace, solvetol = solvetol,
              MASS = MASS, eps = eps, tau = tau)

  if (nls) {
    S <- z$sigma
    z <- nlsur( eqns = eqns, data = data, startvalues = z$coefficients, S = S,
                debug = debug, nls = nls, trace = trace, solvetol = solvetol,
                MASS = MASS, eps = eps, tau = tau)
  }
  z$nlsur <- "NLS"

  # For w/e kind of reason, Stata estimates a second nls with diag(S) instead of
  # I.

  # Estimation of FGNLS
  if (fgnls) {
    # fgnls
    if (trace)
      cat("-- FGNLS\n")

    # nlserg <<- z

    S <- z$sigma
    # print(S)

    z <- nlsur(eqns = eqns, data = data, startvalues = z$coefficients,
               S = S, debug = debug, nls = FALSE, trace = trace,
               solvetol = solvetol, MASS = MASS, eps = eps, tau = tau)


    z$nlsur <- "FGNLS"

    # Estimation of IFGNLS
    if (ifgnls) {

      if (trace)
        cat("-- IFGNLS\n")

      S <- z$sigma
      conv <- FALSE
      iter <- 0
      while (!conv)
      {

        if (iter == 1000){
          message(paste(iter, "nls iterations and convergence not reached."),
                  paste("Last theta is: \n"),
                  coef(z))
          return(0)
        }

        z.old <- z
        S.old <- S

        # s <- chol(qr.solve(S))
        # rss.old <- calc_ssr(r, s, neqs)

        # print(S)

        z <- nlsur(eqns = eqns, data = data, startvalues = z$coefficients,
                   S = S, debug = debug, nls = FALSE, solvetol = solvetol,
                   MASS = MASS, eps = eps, tau = tau)

        r <- z$residuals
        S <- z$sigma
        s <- chol(qr.solve(S))


        rss <- calc_ssr(r, s, neqs)

        # eps <- 1e-5; tau <- 1e-3;
        iter <- iter +1

        maxthetachange <- max(abs(coef(z.old) - coef(z)) /
                                ( abs(coef(z.old)) +1) )
        maxSigmachange <- max(abs(S.old - S) /
                                (abs(S.old) + 1))

        if (is.nan(maxSigmachange))
          maxSigmachange <- 0

        #         print(maxthetachange)
        #         print(maxSigmachange)

        # conv1 <- abs(rss.old - rss) < eps * (rss.old + tau)
        # conv2 <- norm(as.matrix(z.old$coefficients - z$coefficients)) <
        #   eps * (norm(as.matrix(z.old$coefficients)) + tau)

        # conv <- any(conv1, conv2)

        # print(maxthetachange)
        # print(maxSigmachange)

        conv1 <- maxthetachange < eps
        conv2 <- maxSigmachange < ifgnlseps

        conv <- any(conv1, conv2)

        if (trace)
          cat("Iteration", iter, ": SSR", rss, "\n")

      }
      message <- paste("Convergence after iteration:", iter,".")

      S <- z$sigma
      N <- n
      M <- nrow(S)

      LL <- -(M*N)/2 * (1 + log(2*pi)) - N/2 * log(det(S))

      z$message <- message
      z$LL <- LL
      z$sigma <- S
      z$nlsur <- "IFGNLS"

    }
  }


  #### 2. Estimation of covariance matrix, standard errors and t-values ####
  xi      <- list()
  ri      <- list()
  lhs     <- list()
  rhs     <- list()
  n       <- vector("integer", length=neqs)      # number of observations in each equation
  k       <- vector("integer", length=neqs)      # number of (unrestricted) coefficients/
  # regressors in each equation
  df      <- vector("integer", length=neqs)      # degrees of freedom
  ssr     <- vector("numeric", length=neqs)      # sum of squared residuals
  mse     <- vector("numeric", length=neqs)      # mean square error
  rmse    <- vector("numeric", length=neqs)      # root of mse
  mae     <- vector("numeric", length=neqs)      # mean absolute error
  r2      <- vector("numeric", length=neqs)      # R-squared value
  adjr2   <- vector("numeric", length=neqs)      # adjusted R-squared value

  # Get theta, lhs, rhs and residuals from the last estimation.
  theta   <- z$coefficients
  lhs     <- z$lhs
  rhs     <- z$rhs
  xi      <- z$xi

  for (i in 1:length(theta)) {
    name <- names(theta)[i]
    val <- theta[i]
    storage.mode(val) <- "double"
    assign(name, val)
  }

  # contains some duplicated code. ToDo: figure out a way to use this in a
  # efficient manor aka run this loop as few times as possible.
  for (i in 1:neqs) {
    ri[[i]]  <- lhs[[i]] - rhs[[i]]

    n[i]     <- length(lhs[[i]])
    k[i]     <- qr(xi[[i]])$rank
    df[i]    <- n[i] - k[i]

    ssr[i]   <- as.vector(crossprod(ri[[i]]))
    mse[i]   <- ssr[i] / n[i]
    rmse[i]  <- sqrt(mse[i])
    mae[i]   <- sum(abs(ri[[i]]))/n[i]

    r2[i]    <- 1 - ssr[i] /
      ((crossprod(lhs[[i]])) - mean(lhs[[i]]) ^ 2 * n[i])
    adjr2[i] <- 1 - ((n[i] - 1) / df[i]) * (1 - r2[i])
  }
  x <- do.call(cbind, xi)
  r <- do.call(cbind, ri)

  # Create zi for export per equation results
  zi       <- list()
  zi$ssr   <- ssr
  zi$mse   <- mse
  zi$rmse  <- rmse
  zi$mae   <- mae
  zi$n     <- n
  zi$k     <- k
  zi$df    <- df
  zi$r2    <- r2
  zi$adjr2 <- adjr2

  # Estimate covb
  sigma <- z$sigma
  qS <- qr.solve(sigma)

  # covb is solve(XDX)
  covb <- calc_reg(x, r, qS, length(theta), neqs, 0)

  # Estimate SE and t-value
  se <- sqrt(diag(covb))
  tval <- theta / se

  z$se <- se
  z$t <- tval
  z$residuals <- r
  z$covb <- covb
  z$sigma <- sigma
  z$zi <- zi
  z$model <- eqns

  z
}

# [Gallant, A. Ronald (1987): Nonlinear Statistical Models. Wiley: New York]

#' @export
print.nlsur <- function(x) {
  print(x$coefficients)
}

#' @export
summary.nlsur <- function(x) {

  z <- x

  n    <- z$zi$n
  k    <- z$zi$k
  rmse <- z$zi$rmse
  r2   <- z$zi$r2
  zi   <- cbind(n, k, rmse, r2)
  dimnames(zi) <- list(as.character(x$eqnames),
                       c("n", "k", "RMSE", "R-squared"))


  est  <- z$coefficients
  se   <- z$se
  t    <- z$t
  n    <- z$zi$n[1]
  k    <- z$zi$k[1]
  df   <- z$zi$df[1]
  r    <- z$residuals
  prob <- 2 * (1 - pt(abs(t), (n - k)))

  ans <- NULL

  # ans$coefficients <- z[c("coefficients", "se", "t")]
  ans$coefficients <- cbind(est, se, t, prob)
  dimnames(ans$coefficients) <- list(
    names(x$coefficients),
    c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  )

  ans$residuals <- r
  ans$df        <- df
  ans$nlsur     <- x$nlsur
  ans$zi        <- zi

  if (ans$nlsur == "IFGNLS")
    ans$LL <- x$LL

  attr(ans, "class") <- "summary.nlsur"

  ans
}

#' @export
print.summary.nlsur <- function(x) {
  cat("NLSUR Object of type:", x$nlsur, "\n\n")
  print(x$zi)
  cat("\n")
  cat("Coefficientients:\n")
  coefs <- x$coefficients
  printCoefmat(coefs)

  if (x$nlsur == "IFGNLS")
    cat("Log-Likelihood:", x$LL, "\n")
}

#' @export
predict.nlsur <- function(obj, data) {

  eqs <- obj$model

  ddcoef <- list()
  for(i in seq(coef(obj))){
    nam <- names(coef(obj))[i]
    ddcoef[nam] <- coef(obj)[i]
  }
  ddcoef <- data.frame(ddcoef)

  data2 <- data.frame(data, ddcoef)

  fit <- list()
  vnam <- NULL
  for (i in seq(length(eqs))){
    fit[[i]] <- eval(eqs[[i]][[3]], envir = data2)
    vnam.i <- eqs[[i]][[2]]
    vnam <- c(vnam, vnam.i)
  }

  fit <- data.frame(fit)
  names(fit) <- vnam

  fit
}
