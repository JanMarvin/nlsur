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

#' @export
nlsur <- function(eqns, data, startvalues, S = NULL, debug = FALSE,
                  nls = FALSE, fgnls = FALSE, ifgnls = FALSE,
                  MASS = FALSE, trace = FALSE,
                  solvetol = .Machine$double.eps)
{

  #   # Allow passing a single formula to nlsur
  #   if (length(eqns)==1 & formula(eqns))
  #     eqns <- list(eqns)

  if ((MASS & is.null(S)) | is.null(S)) {
    S <- Matrix::kronecker(
      Matrix::Diagonal(length(eqns)),
      Matrix::Diagonal(nrow(data)))
  }

  z      <- list()
  itr    <- 0
  conv   <- FALSE
  nobs   <- nrow(data)

  lhs    <- list()
  rhs    <- list()
  residi <- NULL       # residuals equation wise
  r      <- NULL       # stacked residuals
  # partial derivatives of the residuals with respect to the parameters
  dResidTheta  <- NULL
  dResidThetai <- list()
  eps <- sqrt(.Machine$double.eps)
  # eps    <- 1e-10
  tau    <- 1e-3
  eqnames <- NULL

  # set initial theta
  theta <- startvalues

  ## assign thetas = makes them available for eval
  for (i in 1:length(theta)) {
    name <- names(theta)[i]
    val <- theta[i]
    storage.mode(val) <- "double"
    assign(name, val)
  }

  # Evaluate inital lhs, rhs, residi, r and dResidThetai and dResidTheta
  for (i in 1:length(eqns)) {
    eqnames <- c(eqnames, as.formula(eqns[[i]])[[2L]])
    lhs[[i]] <- as.matrix(eval(as.formula(eqns[[i]])[[2L]], envir = data))
    rhs[[i]] <- as.matrix(eval(as.formula(eqns[[i]])[[3L]], envir = data))

    residi[[i]] <- lhs[[i]] - rhs[[i]]
    r <- rbind(r, residi[[i]])

    dResidThetai[[i]] <-
      attributes(with(data, with(as.list(theta),
                                 eval(deriv(eqns[[i]], names(theta)),
                                      envir = data))))$gradient

    dResidTheta <- rbind(dResidTheta, dResidThetai[[i]])
  }

  if ( any (is.nan(dResidTheta))){
    # eval might return NaN. lm.gls will complain, since r is smaller than S.
    # Fix this by changing its value to zero.
    dResidTheta[is.nan(dResidTheta)] <- 0
    warning("1. Fix NaN value in dResidTheta!")
  }

  # Estimate SSR
  # t(r)%*%S%*%r
  if ( nls & is.null(S) ){
    ssr.old <- as.vector( Matrix::crossprod(r))
  } else {
    ssr.old <- as.vector( Matrix::crossprod(
      Matrix::t(Matrix::crossprod(r, S)), r) )
  }

  theta.new  <- 1
  itr     <- 0
  alpha   <- 1       # stepsizeparameter

  while (!conv) {

    if (itr == 1000){
      message(paste(itr, "nls iterations and convergence not reached."))
      # stop(itr, " nls iterations and convergence not reached.")

      return(0)
    }

    iter  <- FALSE    # gauss-newton iterative loop

    # If alpha < 1 increase it again. Spotted in nls.c
    # if (alpha < 1) alpha  <- alpha*2
    alpha <- min(2*alpha, 1)

    ssr <- ssr.old + 1

    while ( ssr > ssr.old )
      # while ( !iter )
    { # begin iter

      if (debug)
        cat("alpha: ", alpha, "\n")

      # A few random fixes for values of the function. It these are applied,
      # results are biased. These fixes allow estimation of the Sigma Matrix
      # and IFGNL() can be estimated.
      if ( any (is.na(dResidTheta))){
        dResidTheta[is.na(dResidTheta)] <- 0
        warning("2. Fix NA value in dResidTheta!")
      }
      if ( any (is.nan(dResidTheta))){
        dResidTheta[is.nan(dResidTheta)] <- 0
        warning("2. Fix NaN value in dResidTheta!")
      }
      if(any(is.infinite(dResidTheta))){
        # dResidTheta[is.infinite(dResidTheta)] <- 1
        dResidTheta[dResidTheta==Inf] <- +2^1022
        dResidTheta[dResidTheta==-Inf] <- -2^1022
        warning("Fix Inf value in dResidTheta!")
      }
      if(any(is.infinite(r))){
        # r[is.infinite(r)] <- 1
        r[r==Inf] <- +2^1022
        r[r==-Inf] <- -2^1022
        warning("Fix Inf value in r!")
      }
      if(any(is.nan(r))){
        r[is.nan(r)] <- 0
        warning("Fix NaN value in r!")
      }
      if(any(is.na(r))){
        r[is.na(r)] <- 0
        warning("Fix NA value in r!")
      }

      # weighted regression r ~ gradient*hessian
      if (!MASS){
        if ( nls & is.null(S) )
        {
          # gH <- as.matrix(coef(lm(r~dResidTheta+0)))
          gH <- qr.coef(qr(dResidTheta), r)
        } else {
          gH <- qr.solve(
            # t(dResidTheta)%*%S%*%dResidTheta
            Matrix::crossprod(
              Matrix::t(Matrix::crossprod(dResidTheta, S)),
              dResidTheta),
            tol = solvetol
          ) %*% (
            # t(dResidTheta)%*%S%*%r
            Matrix::crossprod(
              Matrix::t(Matrix::crossprod(dResidTheta, S)), r)
          )
        }
      } else{
        gH <- as.matrix(coef(MASS::lm.gls(r ~ 0 + dResidTheta, W = S)))
        # Note: Das funktioniert nicht, mit meiner gigantischen S Matrix. Damit
        # die Gewichte beruecksichtigt werden, wird zunächst mit der Matrix ge-
        # wichtet. Dabei wird der eigen() aufgerufen. Das produziert eine normale
        # Matrix. Meine Matrix::Diagonal() wird mit as.matrix() behandelt und das
        # RAM Problem taucht wieder auf.
      }

      # Sometimes gH will return a NA value. To get resonable results when
      # estimating the new theta, NA will be replaced by a Zero. theta.new will at
      # least contain the value of theta.
      if ( any(is.na(gH)) ){
        gH[is.na(gH)] <- 0
        warning("Fix NA value in gh.")
      }

      # d <- dResidTheta

      #       S2 <- Diagonal(nrow(S))
      #       diag(S2) <- sqrt(diag(S))

      #       X <- dResidTheta
      #       y <- r
      #       eW <- eigen(S) #eigs(A=S, k = nrow(X), which="LM")
      #
      #       d <- eW$values
      #       A <- diag(d^ifelse(FALSE, -0.5, 0.5)) %*% t(eW$vector)
      # Ainv <- eW$vector %*% diag(d^ifelse(inverse, 0.5, -0.5))


      # gH <- qr.coef(qr(A%*%X), A%*%y)
      # gH <- as.matrix(coef(lm(as.matrix(y)~as.matrix(X)+0)))
      # gH <- as.matrix(coef(MASS::lm.gls(r~dResidTheta+0, W = S)))

      # estimate a new theta
      # old theta + scaling-parameter * gH
      theta.new <- as.vector( theta + alpha * gH )
      names(theta.new) <- names(theta)

      if (debug){
        b <- cbind(theta, theta.new)
        print(b)
      }

      ## assign new thetas thetas = makes them available to eval
      for (i in 1:length(theta.new)) {
        name <- names(theta.new)[i]
        val <- theta.new[i]
        storage.mode(val) <- "double"
        assign(name, val)
      }

      # eval eqn with the new theta
      lhs    <- list()
      rhs    <- list()
      residi <- NULL
      r      <- NULL
      dResidTheta  <- NULL
      dResidThetai <- list()

      for (i in 1:length(eqns)) {
        lhs[[i]] <- as.matrix(eval(as.formula(eqns[[i]])[[2]], envir = data))
        rhs[[i]] <- as.matrix(eval(as.formula(eqns[[i]])[[3]], envir = data))

        residi[[i]] <- lhs[[i]] - rhs[[i]]
        r <- rbind(r, residi[[i]])

        dResidThetai[[i]] <-
          attributes(with(data, with(as.list(theta.new),
                                     eval(deriv(eqns[[i]], names(theta.new)),
                                          envir = data))))$gradient

        dResidTheta <- rbind(dResidTheta, dResidThetai[[i]])
      }

      if ( any (is.na(dResidTheta))){
        dResidTheta[is.na(dResidTheta)] <- 0
        warning("2. Fix NA value in dResidTheta!")
      }
      if ( any (is.nan(dResidTheta))){
        dResidTheta[is.nan(dResidTheta)] <- 0
        warning("3. Fix NaN value in dResidTheta!")
      }
      if(any(is.infinite(dResidTheta))){
        # dResidTheta[is.infinite(dResidTheta)] <- 1
        dResidTheta[dResidTheta==Inf] <- +2^1022
        dResidTheta[dResidTheta==-Inf] <- -2^1022
        warning("Fix Inf value in dResidTheta!")
      }
      if(any(is.infinite(r))){
        # r[is.infinite(r)] <- 1
        r[r==Inf] <- +2^1022
        r[r==-Inf] <- -2^1022
        warning("Fix Inf value in r!")
      }
      if(any(is.nan(r))){
        r[is.nan(r)] <- 0
        warning("Fix NaN value in r!")
      }
      if(any(is.na(r))){
        r[is.na(r)] <- 0
        warning("Fix NA value in r!")
      }

      # Estimate a new SSR
      # t(r)%*%S%*%r
      if ( nls & is.null(S) ){
        ssr <- as.vector(crossprod(r))
      } else {
        ssr <- as.vector( Matrix::crossprod(
          Matrix::t(Matrix::crossprod(r, S)), r) )
      }

      # divide stepsizeparameter
      alpha <- alpha/2

      # for nls iteration stops if SSR(betaN) < SSR(beta)
      # else the alogrithm tries to maximize ssr

      #       if ( ifgnls )
      #         iter <- !isTRUE(ssr > ssr.old)
      #       else
      #         iter <- !isTRUE(ssr >= ssr.old)


    } # end iter

    if (trace)
      cat("SSR: ", ssr, "\n")

    # cat (alpha, "\n")
    # cat("SSR: ", ssr, " SSR Old: ", ssr.old, "\n")
    # cat("\nSSR: ", ssr, " SSR Old: ", ssr.old, "\n")

    if(debug){
      print(warnings())
      b <- cbind(theta, theta.new)
      print(b)
    }

    # Stopping rule. [Gallant (1987) p.29]
    # Note: R uses a different convergence criterium

    #     # ssr: |ssr.old - ssr| < eps | ssr.old + tau|
    #     conv1 <- abs(ssr.old - ssr) < eps * (ssr.old + tau)
    #
    #     # theta: ||theta - theta.new|| < eps (||theta|| + tau)
    #         conv2 <- norm(as.matrix(theta - theta.new)) <
    #           eps * (norm(as.matrix(theta)) + tau)

    # no idea why, but Stata uses this
    # Stata uses this
    conv1 <- !isTRUE(abs(ssr.old - ssr) > eps * (ssr.old + tau))
    # conv1 <- TRUE

    conv2 <- !isTRUE(all( alpha * abs(gH) > eps * (abs(theta) + tau) ))
    # conv2 <- TRUE

    # and this is what Stata documents what they do for nl
    # conv2 <- all( alpha * abs(theta.new) <= eps * (abs(theta) + tau) )

    # both convergence criteria must be TRUE
    conv <- all(conv1, conv2)

    if(debug)
      print(b)

    itr <- itr + 1
    theta <- theta.new
    ssr.old <- ssr



    if(debug)
      print(itr)

    # not needed?
    if ( conv ) {
      residi <- matrix(unlist(residi), ncol = length(eqns))
    }

  }

  z$coefficients <- theta
  z$residuals <- residi
  z$lhs <- eqnames

  attr(z, "class") <- "nlsur"

  z

}

#' #### fgnls ####
#' #' @export
#' fgnls <- function(eqns, data, startvalues, S = NULL, debug = FALSE,
#'                   trace = FALSE, solvetol = .Machine$double.eps) {
#'
#'   erg <- nlsur( eqns = eqns, data = data, startvalues = startvalues,
#'                 debug = debug, nls = TRUE, trace = trace,
#'                 solvetol = solvetol)
#'
#'   S <- Matrix::kronecker((qr.solve(1/nrow(data) *
#'                                      Matrix::crossprod(erg$residuals),
#'                            tol = solvetol)),
#'                  Matrix::Diagonal(nrow(data)))
#'
#'   erg <- nlsur(eqns = eqns, data = data, startvalues = erg$coefficients, S = S,
#'                debug = debug, fgnls = TRUE, trace = trace,
#'                solvetol = solvetol)
#'
#'   erg$nlsur <- "FGNLS"
#'
#'   return(erg)
#' }

#' @export
ifgnls <- function(eqns, data, startvalues, type=NULL, S = NULL, debug = FALSE,
                   trace = FALSE, solvetol = .Machine$double.eps, nls = nls,
                   MASS = FALSE) {

  fgnls  <- FALSE
  ifgnls <- FALSE
  z      <- NULL
  zi     <- NULL

  #
  if (!is.null(type)) {
    if(type == "NLS" | type == 1) {
      type <- NULL
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
              MASS = MASS)

  z$nlsur <- "NLS"

  # Estimation of FGNLS
  if (fgnls) {
    # fgnls
    if (trace)
      cat("-- FGNLS\n")

    S <- Matrix::kronecker( qr.solve(1/nrow(data) *
                                       Matrix::crossprod(z$residuals),
                                     tol = solvetol),
                            Matrix::Diagonal(nrow(data)) )

    z <- nlsur(eqns = eqns, data = data, startvalues = z$coefficients,
               S = S, debug = debug, nls = FALSE, trace = trace,
               solvetol = solvetol, MASS = MASS)


    z$nlsur <- "FGNLS"

    # Estimation of IFGNLS
    if (ifgnls) {

      if (trace)
        cat("-- IFGNLS\n")

      conv <- FALSE
      iter <- 0
      while (!conv)
      {

        r <- matrix(z$residuals, ncol=1)

        z.old <- z
        rss.old <- as.vector(Matrix::crossprod(
          Matrix::t(Matrix::crossprod(r, S)), r))
        # rss.old <- sum(S)

        S <- Matrix::kronecker(qr.solve(1/nrow(data) *
                                          Matrix::crossprod(z$residuals),
                                        tol = solvetol),
                               Matrix::Diagonal(nrow(data)))

        z <- nlsur(eqns = eqns, data = data, startvalues = z$coefficients,
                   S = S, debug = debug, nls = FALSE, solvetol = solvetol,
                   MASS = MASS)

        r <- matrix(z$residuals, ncol=1)
        rss <- as.vector(Matrix::crossprod(Matrix::t(Matrix::crossprod(r, S)), r))
        # rss <- sum(S)

        eps <- 1e-5; tau <- 1e-3; iter <- iter +1

        conv1 <- abs(rss.old - rss) < eps * (rss.old + tau)
        conv2 <- norm(as.matrix(z.old$coefficients - z$coefficients)) <
          eps * (norm(as.matrix(z.old$coefficients)) + tau)

        conv <- all(conv1, conv2)

        if (trace)
          cat("Iteration", iter, ": SSR", rss, "\n")

      }
      message <- paste("Convergence after iteration:", iter,".")
      # nlsur3 <<- erg

      S <- 1/nrow(data) * crossprod(z$residuals)
      N <- nrow(data)
      M <- nrow(S)

      LL <- -(M*N)/2 * (1 + log(2*pi)) - N/2 * log(det(S))

      z$message <- message
      z$LL <- LL
      z$sigma <- S
      z$nlsur <- "IFGNLS"

    }
  }


  #### 2. Estimation of covariance matrix, standard errors and t-values ####
  theta <- z$coefficients

  # get the values of the parameters
  for (i in 1:length(theta)) {
    name <- names(theta)[i]
    val <- theta[i]
    storage.mode(val) <- "double"
    assign(name, val)
  }

  ## get the rank for the eqns, compute the first-stage
  ## cov matrix to finish the SUR and 3SLS methods
  X       <- NULL
  r       <- NULL
  residi  <- list()
  derivs  <- list()
  lhs     <- list()
  rhs     <- list()
  G       <- length(eqns)
  n       <- array(0, c(G))      # number of observations in each equation
  k       <- array(0, c(G))      # number of (unrestricted) coefficients/regressors in each equation
  df      <- array(0, c(G))      # degrees of freedom in each equation
  ssr     <- array(0, c(G))      # sum of squared residuals of each equation
  mse     <- array(0, c(G))      # mean square error (residuals) of each equation
  rmse    <- array(0, c(G))      # root of mse
  r2      <- array(0, c(G))      # R-squared value
  adjr2   <- array(0, c(G))      # adjusted R-squared value

  # you're working on parsing out the parameters and the estiamtes for the return structure...
  for (i in 1:length(eqns)) {
    lhs[[i]] <- as.matrix(eval(as.formula(eqns[[i]])[[2]], envir = data))
    rhs[[i]] <- as.matrix(eval(as.formula(eqns[[i]])[[3]], envir = data))
    residi[[i]] <- lhs[[i]] - rhs[[i]]
    derivs[[i]] <- deriv(as.formula(eqns[[i]]), names(theta))

    # computing the jacobian for OLS
    jacobian <- attr(eval(derivs[[i]], envir = data), "gradient")

    n[i]     <-  length(lhs[[i]])
    k[i]     <- qr(jacobian)$rank
    df[i]    <- n[i] - k[i]

    ssr[i]   <- as.vector(crossprod(residi[[i]]))
    mse[i]   <- ssr[i] / (n[i] - k[i])
    rmse[i]  <- sqrt(mse[i])

    r2[i]    <- 1 - ssr[i] / ((crossprod(lhs[[i]])) - mean(lhs[[i]]) ^ 2 * n[i])
    adjr2[i] <- 1 - ((n[i] - 1) / df[i]) * (1 - r2[i])

    X        <- rbind(X, jacobian)
    r        <- cbind(r, residi[[i]])
  }

  zi <- list()
  zi$ssr <- ssr
  zi$mse <- mse
  zi$rmse <- rmse
  zi$n <- n
  zi$k <- k
  zi$df <- df
  zi$r2 <- r2
  zi$adjr2 <- adjr2

  sigma <- 1/nrow(data) * crossprod(r)
  S <- Matrix::kronecker(
    qr.solve(sigma), Matrix::Diagonal(n[1])
  )

  # Estimate covb
  covb <- qr.solve(
    Matrix::crossprod(Matrix::t(Matrix::crossprod(X, S)),X),
    tol = solvetol
  )
  colnames(covb) <- rownames(covb)

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
  dimnames(zi) <- list(as.character(x$lhs),
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
  dimnames(ans$coefficients) <- list(names(x$coefficients),
                                     c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))

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

