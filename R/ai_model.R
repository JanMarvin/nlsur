#' Q/AI model function
#'
#' Modelfunction to create Deatons and Muellbauers (1980) famous
#' Almost-Ideal Demand System or the Quadratic Almost-Ideal Demand System by
#' Banks et al. (1997).
#'
#' @param w character vector of m budgetshares used for estimation.
#' @param p character vector of m prices.
#' @param exp single character vector of total expenditure.
#' @param alph0 start value for translog price index.
#' @param logp logical if prices are log prices.
#' @param logexp logical if expenditure is log expenditure.
#' @param priceindex character either "translog" or "S" for the stone price
#'  index.
#' @param modeltype character either "AI" or "QAI" for AI or QAI model.
#' @param ray logical if Ray (1983) scaling should be included. If TRUE requires
#' demographic vector.
#' @param demogr character vector of k demographic variables.
#'
#' @references Deaton, Angus S., Muellbauer, John: An Almost Ideal Demand
#'  System, The American Economic Review 70(3), American Economic Association,
#'  312-326, 1980
#' @references Banks, James, Blundell, Richard, Lewbel, Arthur: Quadratic Engel
#'  Curves and Consumer Demand, The Review of Economics and Statistics 79(4),
#'  The MIT Press, 527-539, 1997
#' @references Ray, Ranjan: Measuring the costs of children: An alternative
#'  approach, Journal of Public Economics 22(1), 89-102, 1983
#'
#' @seealso ai and qai
#'
#' @export
ai.model <- function(w, p, exp, alph0 = 10, logp = TRUE, logexp = TRUE,
                     priceindex = "translog", modeltype = "AI",
                     ray = FALSE, demogr) {

  if (missing(demogr) & ray)
    stop("Demogr. scaling selected, but no demographic variables specified.")

  if (missing(demogr))
    demogr <- character(0)


  #### Create AI equation system ####
  # w_i = alpha + \sum gamma log(p) + beta log(exp/P)

  neqs <- length(w)

  #### alpha ####
  # restricted
  alpha <- NULL
  alpha <- paste0("a", sprintf("%02d",1:neqs))
  if (priceindex == "translog" | priceindex == "S")
    alpha[neqs] <- paste0("(1-",paste(alpha[-neqs], collapse="-"),")")

  #### gammas ####
  # gammas are restricted.
  # gammas create a matrix.
  # row sums add to zero. gammas are positive semidefinit.
  # g11 g12 -g11-g12
  # g12 g22 -g12-g22
  gamma <- NULL
  gamma <- matrix(NA, neqs, neqs)
  for( i in 1:neqs) {
    for (j in 1:neqs) {
      gamma[i,j] <- paste0("g",sprintf("%02d",i),sprintf("%02d",j))
    }
  }

  # positive semidefinit
  gamma[lower.tri(gamma)] <- t(gamma)[lower.tri(gamma)]

  # make gammas add up to zero
  for (i in 1:neqs)
    gamma[i,neqs] <- paste0("(-",paste(gamma[i,-neqs], collapse="-"),")")
  for (i in 1:neqs)
    gamma[neqs,i] <- paste0("(-",paste(gamma[-neqs,i], collapse="-"),")")

  #### logp ####
  if ( !isTRUE(logp) ) {
    logp <- paste0("log(",p,")")
  } else {
    logp <- p
  }

  #### beta ####
  # restricted to zero
  beta <- NULL
  beta <- paste0("b", sprintf("%02d",1:neqs))
  beta[neqs] <- paste0("(-",paste(beta[-neqs], collapse="-"),")")

  if (ray) {

    # rho is the scaling parameter
    rho <- paste0("rho_",demogr)

    m0 <- paste("( 1 + ", paste(rho, "*", demogr, collapse = " + "), ")")

    # log(m0)+translog
    m0 <- paste("log(", m0, ")")

    # sum ( eta_i ) = 0.
    eta <- matrix(NA, length(demogr), neqs )
    for( i in seq_along(demogr)) {
      for (j in 1:neqs) {
        eta[i,j] <- paste0("eta_",demogr[i],"_",sprintf("%02d",j))
      }
    }
    for (i in 1:length(demogr))
      eta[i,neqs] <- paste0("- (", paste(eta[i,-neqs], collapse= " + "), ")")

    # eta(z)
    etaz <- matrix(NA, length(demogr), neqs)
    # eta_i * z
    for (i in 1:neqs)
      etaz[,i] <- paste0(eta[,i], "*", demogr)

    # c(p,z)
    cpz <- paste("(",
                 paste( paste("exp(", logp, "*", t(etaz), ")"),
                        collapse = " * "),
                 ")")


    # beta* is (beta_i + eta_i * z)
    beta_star <- rep(NA, neqs)

    for (j in 1:neqs) {
      beta_star[j] <-  paste( beta[j], " + ",
                              paste( etaz[,j], collapse = " + "),
                              collapse = "+")
      beta_star[j] <- paste("(", beta_star[j], ")")
    }

  }

  #### log(a) ####
  alpha0 <- alph0
  alogp <- paste(paste(alpha,"*", logp), collapse=" + ")

  #### lambda ####
  # restricted to sum lambda = 0
  lambda <- NULL
  lambda <- paste0("l", sprintf("%02d",1:neqs))
  lambda[neqs] <- paste0("(-",paste(lambda[-neqs], collapse="-"),")")

  #### log(exp) ####
  if ( !isTRUE(logexp) ) {
    logexp <- paste0("log(",exp,")")
  } else {
    logexp <- exp
  }

  #### translog price index ####
  pmatrix <- matrix(NA, neqs, neqs)

  for(i in 1:neqs) {
    pmatrix[i,] <- paste0(logp[i],"*",logp)
  }

  translogmatrix <- matrix(NA, neqs, neqs)
  for(i in 1:neqs)
    translogmatrix[i,] <- paste0(gamma[i,], "*", pmatrix[i,])

  translog <- paste("(", alpha0, "+", alogp, "+ 0.5 * (",
                    paste0(translogmatrix, collapse = " + "), ")", ")")

  # Stone price index
  if (priceindex == "S") {
    translog <- paste("(", paste(w, "*", logp, collapse = " + ") ,")" )
  }

  #### cobb douglas price aggregator ####
  cobbdoug <- paste("(", paste(paste("exp(", logp, "*", beta, ")"),
                               collapse = " * "), ")")

  #### combine everything to the final model ####
  eqs <- list()

  if (modeltype == "AI" & !ray) {

    #### eqs ####
    for (i in 1:neqs) {
      glogp <- paste(paste0(gamma[i,],"*", logp), collapse=" + ")
      eqs[i] <- paste(w[i], "~", alpha[i], "+", glogp, "+", beta[i], "*",
                      "(", logexp, "-", translog, ")" )
    }
  }

  if (modeltype == "QAI" & !ray) {

    #### eqs ####
    for (i in 1:neqs) {
      glogp <- paste(paste0(gamma[i,],"*", logp), collapse=" + ")
      eqs[i] <- paste(w[i], "~", alpha[i], "+", glogp, "+", beta[i], "*",
                      "(", logexp, "-", translog, ") +", lambda[i], "* ((",
                      logexp, "-", translog, ")^2 /", cobbdoug,")")
    }
  }

  if (modeltype == "AI" & ray) {

    translog_star <- paste("(", m0, "+", translog, ")")

    #### eqs ####
    for (i in 1:neqs) {
      glogp <- paste(paste0(gamma[i,],"*", logp), collapse=" + ")
      eqs[i] <- paste(w[i], "~", alpha[i], "+", glogp, "+", beta_star[i], "*",
                      "(", logexp, "-", translog_star, ")" )
    }
  }

  if (modeltype == "QAI" & ray) {

    translog_star <- paste("(", m0, "+", translog, ")")

    #### eqs ####
    for (i in 1:neqs) {
      glogp <- paste(paste0(gamma[i,],"*", logp), collapse=" + ")
      eqs[i] <- paste(w[i], "~", alpha[i], "+", glogp, "+", beta_star[i], "*",
                      "(", logexp, "-", translog_star, ") +", lambda[i], "* ((",
                      logexp, "-", translog_star, ")^2 /",
                      "(", cobbdoug, "*", cpz, ")", " )")
    }
  }

  if (modeltype == "eQAI" & !ray) {

    mue <- paste0("mue", 1:neqs)

    for (i in 1:neqs) {

      eqs[i] <- paste(mue[i], "~", 1, "+ 1 /", w[i], "*",
                      "(", beta[i], " +  2 *", lambda[i], "* ((",
                      logexp, "-", translog, ") /",
                      "(", cobbdoug, ")", " ))")
    }
  }

  if (modeltype == "eQAI" & ray) {

    mue <- paste0("mue", 1:neqs)
    translog_star <- paste("(", m0, "+", translog, ")")

    for (i in 1:neqs) {

      eqs[i] <- paste(mue[i], "~", 1, "+ 1 /", w[i], "*",
                      "(", beta_star[i], " +  2 *", lambda[i], "* ((",
                      logexp, "-", translog_star, ") /",
                      "(", cobbdoug, "*", cpz, ")", " ))")
    }
  }



  eqs <- unlist(eqs)

  model <- list()

  if (!(modeltype == "eQAI" )) {

    # drop one equation
    for (i in 1:(neqs-1))
      model[[i]] <- as.formula(eqs[i])

  } else {

    for (i in 1:(neqs))
      model[[i]] <- as.formula(eqs[i])
  }

  return(model)
}

#' Estimation of an Almost-Ideal Demand System
#'
#' Estimation of an Almost-Ideal Demand System using nlsur().
#'
#' @param w character vector of m budgetshares used for estimation.
#' @param p character vector of m prices.
#' @param x single character vector of total expenditure.
#' @param z character vector of k demographic variables.
#' @param a0 start value for translog price index.
#' @param data data.frame containing the variables.
#' @param scale logical if TRUE Rays (1983) scaling is used.
#' @param logp logical if prices are log prices.
#' @param logexp logical if expenditure is log expenditure.
#' @param ... additional options passed to nlsur
#'
#' @references Deaton, Angus S., Muellbauer, John: An Almost Ideal Demand
#'  System, The American Economic Review 70(3), American Economic Association,
#'  312-326, 1980
#' @references Ray, Ranjan: Measuring the costs of children: An alternative
#'  approach, Journal of Public Economics 22(1), 89-102, 1983
#'
#' @seealso qai and ai.model
#'
#' @export
ai <- function(w, p, x, z, a0 = 0, data, scale = FALSE,
               logp = TRUE, logexp = TRUE, ...) {

  if (missing(z))
    z <- character(0)

  vars <- c(w,p,x,z)
  ndat <- names(data)

  if(!all(vars %in% ndat)) {
    stop("Selected Variable not found in Dataset.")
  }

  if(!scale) {
    model <- ai.model(w = w, p = p, exp = x, alph0 = a0, modeltype = "AI",
                      logp = logp, logexp = logexp)
  } else {
    model <- ai.model(w = w, p = p, exp = x, demogr = z, alph0 = a0,
                      ray = TRUE, modeltype = "AI",
                      logp = logp, logexp = logexp)
  }

  res <- nlsur(eqns = model, data = data, type = 3, ...)

  res
}

#' Estimation of an Quadratic Almost-Ideal Demand System
#'
#' Estimation of an Quadratic Almost-Ideal Demand System using nlsur().
#'
#' @param w character vector of m budgetshares used for estimation.
#' @param p character vector of m prices.
#' @param x single character vector of total expenditure.
#' @param z character vector of k demographic variables.
#' @param a0 start value for translog price index.
#' @param data data.frame containing the variables.
#' @param scale logical if TRUE Rays (1983) scaling is used.
#' @param logp logical if prices are log prices.
#' @param logexp logical if expenditure is log expenditure.
#' @param ... additional options passed to nlsur
#'
#' @references Banks, James, Blundell, Richard, Lewbel, Arthur: Quadratic Engel
#'  Curves and Consumer Demand, The Review of Economics and Statistics 79(4),
#'  The MIT Press, 527-539, 1997
#' @references Ray, Ranjan: Measuring the costs of children: An alternative
#'  approach, Journal of Public Economics 22(1), 89-102, 1983
#'
#' @seealso ai and ai.model
#'
#' @export
qai <- function(w, p, x, z, a0 = 0, data, scale = FALSE,
                logp = TRUE, logexp = TRUE, ...) {

  if (missing(z))
    z <- character(0)

  vars <- c(w,p,x,z)
  ndat <- names(data)

  if(!all(vars %in% ndat)) {
    stop("Selected Variable not found in Dataset.")
  }

  if(!scale) {
    model <- ai.model(w = w, p = p, exp = x, alph0 = a0, modeltype = "QAI",
                      logp = logp, logexp = logexp)
  } else {
    model <- ai.model(w = w, p = p, exp = x, demogr = z, alph0 = a0,
                      ray = TRUE, modeltype = "QAI",
                      logp = logp, logexp = logexp)
  }

  res <- nlsur(eqns = model, data = data, type = 3, ...)

  # required for eQAI
  attr(res, "w") <- w
  attr(res, "p") <- p
  attr(res, "x") <- x
  attr(res, "z") <- z
  attr(res, "logp")   <- logp
  attr(res, "logexp") <- logexp
  attr(res, "scale")  <- scale

  res
}

#' Estimation of elasticies of the Quadratic Almost-Ideal Demand System
#'
#' Estimates the expenditure elasticity for goods
#'
#' @param object qai result
#' @param data data vector used for estimation
#' @param usemean evaluate at mean
#'
#' @references Poi, Brian P.: Easy demand-system estimation with quaids, The
#'  Stata Journal 12(3), 433-446, 2012
#'
#' @seealso ai and ai.model
#'
#' @export
eQAI <- function(object, data, usemean = FALSE) {

  # extract var names
  w <- attr(object, "w")
  p <- attr(object, "p")
  x <- attr(object, "x")
  z <- attr(object, "z")
  logp   <- attr(object, "logp")
  logexp <- attr(object, "logexp")
  scale  <- attr(object, "scale")

  if(!scale) {
    eqs <- ai.model(w = w, p = p, exp = x, modeltype = "eQAI",
                    logp = logp, logexp = logexp)
  } else {
    eqs <- ai.model(w = w, p = p, exp = x, demogr = z,
                    ray = TRUE, modeltype = "eQAI",
                    logp = logp, logexp = logexp)
  }

  if (!usemean) {

    eqns_lhs <- lapply(X = eqs, FUN = function(x)x[[2L]])
    eqns_rhs <- lapply(X = eqs, FUN = function(x)x[[3L]])
    vnam     <- sapply(X = eqns_lhs, FUN = as.character)

    data2 <- data.frame(data, as.list(coef(object)))

    # create fit: predict result
    fit <- lapply(X = eqns_rhs, FUN = eval, envir = data2)

    # replace is.infinite() with NA
    fit <- lapply(X = fit, FUN = function(x) replace(x, is.infinite(x), NA) )

    fit <- data.frame(fit)
    names(fit) <- vnam

  } else {

    fit <- NULL

    # log means are required
    logmean <- function(x) exp ( mean( log( x ) ) )

    # calculate means
    wm <- sapply(w, FUN=function(x) mean(data[[x]], use.na = FALSE))
    pm <- sapply(p, FUN=function(x) mean(data[[x]], use.na = FALSE))
    xm <- sapply(x, FUN=function(x) mean(data[[x]], use.na = FALSE))

    if (!logp) pm <- sapply(p, FUN=function(x) logmean(data[[x]]))
    if (!logexp) xm <- sapply(x, FUN=function(x) logmean(data[[x]]))

    ms <- data.frame(t(c(wm, pm, xm)))

    if(scale) {
      zm <- sapply(z, FUN=function(x) mean(data[[x]], use.na = FALSE))
      ms <- data.frame(t(c(wm, pm, xm, zm)))
    }

    for (eq in eqs) {

      vars <- all.vars(eq[[3]])
      vars <- vars[!(vars %in% names(coef(object)))]

      # replace values in equation with fixed estimates
      for (i in vars) {
        # no underscore ahead and behind i
        eq <- gsub(pattern = paste0("(?<!_)", i, "(?!_)"),
                   replacement = ms[names(ms) == i],
                   x = eq, perl = TRUE)
      }

      fits <- nlcom(object = object, form = eq[3])
      fit <- rbind(fit, fits)
    }

    rownames(fit) <- sapply(eqs, FUN = function(x) x[[2]])

  }

  fit
}
