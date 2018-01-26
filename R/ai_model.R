#' Q/AI model function
#'
#' Modelfunction to create Deatons and Muellbauers (1980) famouse
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

  eqs <- unlist(eqs)

  model <- list()

  # drop one equation
  for (i in 1:(neqs-1))
    # for (i in 1:(neqs))
    model[[i]] <- as.formula(eqs[i])

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

  res
}
