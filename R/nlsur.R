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
#' @param debug is a logical if debug output should be included. This can be a
#' lot and might not be of any help to you.
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
#' @importFrom Matrix kronecker diag
#' @importFrom MASS lm.gls
#' @import RcppArmadillo
#' @useDynLib nlsur
#' @export
.nlsur <- function(eqns, data, startvalues, S = NULL, debug = FALSE,
                   nls = FALSE, fgnls = FALSE, ifgnls = FALSE, qrsolve = FALSE,
                   MASS = FALSE, trace = FALSE, eps = eps, tau = tau,
                   weights = weights)
{
  z    <- list()
  itr  <- 0
  conv <- FALSE

  lhs  <- rhs <- ri <- xi <- list()
  r    <-   x <- NULL

  neqs <- length(eqns)
  n    <- vector("integer", length=neqs)
  k    <- vector("integer", length=neqs)
  df   <- vector("integer", length=neqs)

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
    n[[i]] <- length(lhs[[i]])
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
  ssr.old <- calc_ssr(r, s, weights)

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

    # begin regression
    # Regression of residuals on derivs
    if (nls & qrsolve & is.null(weights)) {

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

      if (MASS & is.null(weights))
      {

        Sigma <- Matrix::kronecker(X = qr.solve(1/n * crossprod(r)),
                                   Y = Matrix::diag(n) )
        r <- matrix(r, ncol = 1)
        x <- do.call(rbind, xi)

        # Use MASS::lm.gls for the weighted regression, this will call eigen()
        # to pre-weight r and x which then can be solved with qr. This can
        # cause huge matrices as lm.gls() is not able to handle sparse Matrices.
        theta.new <- coef(MASS::lm.gls(r ~ 0 + x, W = Sigma))
        names(theta.new) <- names(theta)
        theta <- theta.new
      } else {

        # Weighted regression of residuals on derivs ---
        theta_test <- calc_reg(x, r, qS, weights, length(theta), 1)
        theta.new <- as.vector(theta_test)

        names(theta.new) <- names(theta)
        theta <- theta.new
      }


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
        lhs[[i]] <- eval(as.formula(eqns[[i]])[[2L]], envir = data)
        rhs[[i]] <- eval(as.formula(eqns[[i]])[[3L]], envir = data)
        ri[[i]] <- lhs[[i]] - rhs[[i]]

        xi[[i]] <- attr(with(data, with(as.list(theta.new),
                                        eval(deriv(eqns[[i]], names(theta.new)),
                                             envir = data))), "gradient")

        n[i]     <- length(rhs[[i]])
        k[i]     <- qr(xi[[i]])$rank
        df[i]    <- n[i] - k[i]
      }

      df <- unique(df)
      k  <- unique(k)
      n  <- unique(n)


      if (length(n)>1 | length(k)>1)
        stop ("unequal n or k")

      r <- do.call(cbind, ri)
      x <- do.call(cbind, xi)

      # Evaluate initial ssr
      ssr <- calc_ssr(r, s, weights)

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

  # fitted
  fitted <- as.data.frame(rhs)
  names(fitted) <- eqnames

  z$fitted       <- fitted
  z$coefficients <- theta
  z$residuals    <- r
  z$eqnames      <- eqnames
  z$sigma        <- 1/n * crossprod(r, weights * r)

  z$deviance     <- as.numeric(ssr)
  z$df.residual  <- df

  z$weights      <- weights

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
#' @param eps the epislon used for convergence in nlsur(). Default is 1e-5.
#' @param tau is another convergence variable. Default is 1e-3.
#' @param ifgnlseps is epislon for ifgnls(). Default is 1e-10.
#' @param stata is a logical. If TRUE for nls a second evaluation will be run.
#' Stata does this by default. For this second run Stata replaces the diagonal
#' of the I matrix with the coefficients.
#' @param trace logical wheather or not SSR information should be printed.
#' Default is FALSE.
#' @param debug logical wheater or not debug information will be printed.
#' Default is FALSE.
#' @param S is a weight matrix used for evaluation. If no weight matrix is
#' provided the identity matrix I will be used.
#' @param qrsolve logical
#' @param MASS is a logical wheather the MASS::lm.gls() function should be used
#' for weighted Regression. This can cause sever RAM usage as the weight matrix
#' tend to be huge (n-equations * n-rows).
#' @param weights Additional weight vector.
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
#' @examples # predict(nlsurObj, dataframe)
#' @seealso nls
#' @import RcppArmadillo
#' @useDynLib nlsur
#'
#' @export
nlsur <- function(eqns, data, startvalues, type=NULL, S = NULL, debug = FALSE,
                  trace = FALSE, stata = FALSE, qrsolve = FALSE,
                  weights, MASS = FALSE,
                  eps = 1e-5, ifgnlseps = 1e-10, tau = 1e-4) {

  if (missing(weights))
    weights <- NULL
  else
    weights <- as.character(substitute(weights))


  # Check if all variables that are not startvalues exist in data.
  vars <- unlist(lapply(eqns, all.vars))
  vars <- vars[which(!vars %in% names(startvalues))]
  ok <- all(vars%in%names(data))
  if (!ok) {
    message("Missmatch in model and dataset.")
    return(0)
  }


  # remove observation, if observation a parameter contains NA.
  modelparameters <- c(unlist(lapply(eqns, all.vars)), weights)
  parms <- modelparameters[which(!modelparameters %in% names(startvalues))]

  # Check for weights
  if ( is.null(weights) )
    w <- rep(x = 1, nrow(data))
  else {
    weights <- as.name(weights)

    w <- eval(substitute(weights), data)
  }

  data <- na.omit(data[unique(parms)])

  nls    <- FALSE
  fgnls  <- FALSE
  ifgnls <- FALSE
  z      <- NULL
  n      <- nrow(data)

  # normweights
  w <- w/sum(w) * n

  cl <- match.call()

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

  z <- .nlsur( eqns = eqns, data = data, startvalues = startvalues, S = S,
               debug = debug, nls = TRUE, trace = trace, qrsolve = qrsolve,
               MASS = MASS, eps = eps, tau = tau, weights = w)

  if (nls & stata) {

    # For w/e kind of reason, Stata estimates a second nls with diag(S) instead
    #  of I.
    S <- z$sigma

    z <- .nlsur( eqns = eqns, data = data, startvalues = z$coefficients, S = S,
                 debug = debug, nls = nls, trace = trace, qrsolve = qrsolve,
                 MASS = MASS, eps = eps, tau = tau, weights = w)

    # FixMe: Stata uses this sigma for covb, not the updated?
    z$sigma <- diag(diag(S))
  }
  z$nlsur <- "NLS"


  # Estimation of FGNLS
  if (fgnls) {
    # fgnls
    if (trace)
      cat("-- FGNLS\n")

    # nlserg <<- z

    S <- z$sigma

    z <- .nlsur(eqns = eqns, data = data, startvalues = z$coefficients,
                S = S, debug = debug, nls = FALSE, trace = trace,
                qrsolve = qrsolve, MASS = MASS, eps = eps, tau = tau,
                weights = w)

    # FixMe: Stata uses this sigma for covb, not the updated?
    if (!ifgnls)
      z$sigma <- S

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

        z <- .nlsur(eqns = eqns, data = data, startvalues = z$coefficients,
                    S = S, debug = debug, nls = FALSE,
                    qrsolve = qrsolve, MASS = MASS, eps = eps, tau = tau,
                    weights = w)

        r <- z$residuals
        S <- z$sigma
        s <- chol(qr.solve(S))

        rss <- calc_ssr(r, s, w)

        iter <- iter +1

        maxthetachange <- max(abs(coef(z.old) - coef(z)) /
                                ( abs(coef(z.old)) +1) )
        maxSigmachange <- max(abs(S.old - S) /
                                (abs(S.old) + 1))

        if (is.nan(maxSigmachange))
          maxSigmachange <- 0

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



      z$message <- message
      z$sigma <- S
      z$residuals <- r
      z$nlsur <- "IFGNLS"

    }
  }

  S <- z$sigma
  N <- n
  M <- nrow(S)

  LL <- -(M*N)/2 * (1 + log(2*pi)) - N/2 * log(det(S))
  z$LL <- LL

  z$model <- eqns
  z$data <- data
  z$call <- cl

  if (is.null(weights))
    z$weights <- NULL

  z
}

# [Gallant, A. Ronald (1987): Nonlinear Statistical Models. Wiley: New York]

#' @export
print.nlsur <- function(x, ...) {
  # ... is to please check()
  print(x$coefficients, ...)
}

#' @export
summary.nlsur <- function(object, ...) {
  # ... is to please check()

  # z is shorter
  z <- object

  data <- z$data
  eqns <- z$model
  neqs <- length(eqns)
  w    <- weights(z)

  if (is.null(w))
    w <- rep(1, nrow(data))


  #### 2. Estimation of covariance matrix, standard errors and t-values ####
  xi      <- list()
  ri      <- list()
  lhs     <- list()
  rhs     <- list()
  n       <- vector("integer", length=neqs)      # number of observations in each equation
  k       <- vector("integer", length=neqs)      # number of (unrestricted) coefficients/
  scale   <- vector("numeric", length=neqs)      # scalefactor
  div     <- vector("numeric", length=neqs)      # divisor
  wi      <- vector("numeric", length=neqs)      # normalized weights
  # regressors in each equation
  df      <- vector("integer", length=neqs)      # degrees of freedom
  ssr     <- vector("numeric", length=neqs)      # sum of squared residuals
  mss     <- vector("numeric", length=neqs)
  mse     <- vector("numeric", length=neqs)      # mean square error
  rmse    <- vector("numeric", length=neqs)      # root of mse
  mae     <- vector("numeric", length=neqs)      # mean absolute error
  r2      <- vector("numeric", length=neqs)      # R-squared value
  adjr2   <- vector("numeric", length=neqs)      # adjusted R-squared value

  # Get coefficients from the last estimation.
  est     <- z$coefficients

  for (i in 1:length(est)) {
    name <- names(est)[i]
    val <- est[i]
    storage.mode(val) <- "double"
    assign(name, val)
  }

  # contains some duplicated code.
  for (i in 1:neqs) {
    lhs[[i]]  <- eval(as.formula(eqns[[i]])[[2L]], envir = data)
    rhs[[i]]  <- eval(as.formula(eqns[[i]])[[3L]], envir = data)
    ri[[i]]   <- lhs[[i]] - rhs[[i]]

    xi[[i]] <- attr(with(data, with(as.list(est),
                                    eval(deriv(eqns[[i]], names(est)),
                                         envir = data))), "gradient")

    n[i]     <- length(lhs[[i]])
    k[i]     <- qr(xi[[i]])$rank
    df[i]    <- n[i] - k[i]

    scale[i] <- n[i]/sum(w)
    div[i]   <- n[i] - 1
    wi       <- w/sum(w) * n[i]

    ssr[i]   <- sum( ri[[i]]^2 * w) * scale[i]

    # if (const)
    lhs_wm   <- weighted.mean(x = lhs[[i]], w = wi)
    wvar     <- (1/(n[i] - 1)) * sum( wi * (lhs[[i]] - lhs_wm)^2)
    mss[i] <- wvar * div[i] - ssr[i]
    # else

    mse[i]   <- ssr[i] / n[i]
    rmse[i]  <- sqrt(mse[i])
    mae[i]   <- sum(abs(ri[[i]]))/n[i]

    r2[i]  <- mss[i] / (mss[i] + ssr[i])
    adjr2[i] <- 1 - ((n[i] - 1) / df[i]) * (1 - r2[i])
  }

  x  <- do.call(cbind, xi)
  r  <- do.call(cbind, ri)

  nE <- sum(n) / sum ( w/sum(w) )
  kE  <- sum(k)


  # Estimate covb
  sigma <- z$sigma
  qS <- qr.solve(sigma)

  if (neqs == 1){
    # single eqs: covb is s *(XX)-1 for single equations
    covb <- 1/(n-k) * sum(r^2* w) * scale * qr.solve(Matrix::crossprod(x, w*x))

  } else {
    # covb is solve(XDX)
    covb <- calc_reg(x, r, qS, w, length(est), 0)
  }
  dimnames(covb) <- list(names(est), names(est))

  # Estimate se, tval and prob
  se <- sqrt(diag(covb))
  tval <- est / se
  prob <- 2 * (1 - pt(abs(tval), (nE * kE )))

  # per equation statistics
  zi   <- cbind(n, k, rmse, mae, r2, adjr2)
  dimnames(zi) <- list(as.character(z$eqnames),
                       c("n", "k", "RMSE", "MAE", "R-squared",
                         "Adj-R-sqr."))

  # ans: returned object
  ans <- NULL

  # ans$coefficients <- z[c("coefficients", "se", "t")]
  ans$coefficients <- cbind(est, se, tval, prob)

  cnames <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  if (neqs == 1)
    cnames <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")

  dimnames(ans$coefficients) <- list(
    names(z$coefficients),
    cnames
  )

  ans$residuals    <- residuals
  ans$nlsur        <- z$nlsur
  ans$sigma        <- sigma
  ans$zi           <- zi
  ans$cov          <- covb
  ans$weights      <- weights(z)

  if (ans$nlsur == "IFGNLS")
    ans$LL <- z$LL

  class(ans) <- "summary.nlsur"

  ans
}

#' @export
print.summary.nlsur <- function(x, ...) {
  # ... is to please check()
  cat("NLSUR Object of type:", x$nlsur, "\n\n")

  if (!is.null(weights(x))) {
    cat("Scaled R-squared: \n\n")
  }

  print(x$zi, digits = 4)

  cat("\n")
  cat("Coefficients:\n")

  if (!is.null(weights(x))) {
    cat("Weighted nlsur: \n\n")
  }

  printCoefmat(x$coefficients, digits = 4)

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
  summary(object)$cov
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
#' @details In contrast to other regression objects nlsur does not store the
#' complete model in the resulting object. This requires a data object for
#' predict. A limitation due to the fact that nlsur estimations tend to be
#' estimated on larger data.frames.
#'
#' @examples # predict(nlsurObj, dataframe)
#'
#' @export
predict.nlsur <- function(object, newdata, ...) {

  eqs <- object$model

  if (missing(newdata) || is.null(newdata)) {
    data <- object$data
  } else {
    data <- newdata
  }

  data2 <- data.frame(data, as.list(coef(object)))

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
