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
#' @details While Ray and Stata use log(m0) for the demographic variables, there
#' is no guarantee, that m0 is positive. The model relies much on the correct
#' starting values. Therefore log(abs(m0)) is used.
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

  if (missing(demogr) && ray)
    stop("Demogr. scaling selected, but no demographic variables specified.")

  if (missing(demogr))
    demogr <- character(0)


  #### Create AI equation system ####
  # w_i = alpha + \sum gamma log(p) + beta log(exp/P)

  neqs <- length(w)

  #### alpha ####
  # restricted
  alpha <- NULL
  alpha <- paste0("a", sprintf("%02d", seq_len(neqs)))
  if (priceindex == "translog" || priceindex == "S")
    alpha[neqs] <- paste0("(1-", paste(alpha[-neqs], collapse = "-"), ")")

  #### gammas ####
  # gammas are restricted.
  # gammas create a matrix.
  # row sums add to zero. gammas are positive semidefinit.
  # g11 g12 -g11-g12
  # g12 g22 -g12-g22
  gamma <- NULL
  gamma <- matrix(NA, neqs, neqs)
  for (i in seq_len(neqs)) {
    for (j in seq_len(neqs)) {
      gamma[i, j] <- paste0("g", sprintf("%02d", i), sprintf("%02d", j))
    }
  }

  # positive semidefinit
  gamma[lower.tri(gamma)] <- t(gamma)[lower.tri(gamma)]

  # make gammas add up to zero
  for (i in seq_len(neqs))
    gamma[i, neqs] <- paste0("(-", paste(gamma[i, -neqs], collapse = "-"), ")")
  for (i in seq_len(neqs))
    gamma[neqs, i] <- paste0("(-", paste(gamma[-neqs, i], collapse = "-"), ")")

  #### logp ####
  if (!isTRUE(logp)) {
    logp <- paste0("log(", p, ")")
  } else {
    logp <- p
  }

  #### beta ####
  # restricted to zero
  beta <- NULL
  beta <- paste0("b", sprintf("%02d", seq_len(neqs)))
  beta[neqs] <- paste0("(-", paste(beta[-neqs], collapse = "-"), ")")

  if (ray) {

    # rho is the scaling parameter
    rho <- paste0("rho_", demogr)

    m0 <- paste("( 1 + ", paste(rho, "*", demogr, collapse = " + "), ")")

    # log(m0) + translog
    # m0 >= 0 otherwise log(-1) == NaN
    m0 <- paste("log(abs(", m0, "))")

    # sum ( eta_i ) = 0.
    eta <- matrix(NA, length(demogr), neqs)
    for (i in seq_along(demogr)) {
      for (j in seq_len(neqs)) {
        eta[i, j] <- paste0("eta_", demogr[i], "_", sprintf("%02d", j))
      }
    }
    for (i in seq_len(length(demogr)))
      eta[i, neqs] <- paste0("- (", paste(eta[i, -neqs], collapse = " + "), ")")

    # eta(z)
    etaz <- matrix(NA, length(demogr), neqs)
    # eta_i * z
    for (i in seq_len(neqs))
      etaz[, i] <- paste0(eta[, i], "*", demogr)

    # c(p,z)
    cpz <- paste("(",
                 paste(
                   paste("exp(", logp, "*", t(etaz), ")"),
                   collapse = " * "),
                 ")")


    # beta* is (beta_i + eta_i * z)
    beta_star <- rep(NA, neqs)

    for (j in seq_len(neqs)) {
      beta_star[j] <-  paste(beta[j], " + ",
                             paste(etaz[, j], collapse = " + "),
                             collapse = "+")
      beta_star[j] <- paste("(", beta_star[j], ")")
    }

  }

  #### log(a) ####
  alpha0 <- alph0
  alogp <- paste(paste(alpha, "*", logp), collapse = " + ")

  #### lambda ####
  # restricted to sum lambda = 0
  lambda <- NULL
  lambda <- paste0("l", sprintf("%02d", seq_len(neqs)))
  lambda[neqs] <- paste0("(-", paste(lambda[-neqs], collapse = "-"), ")")

  #### log(exp) ####
  if (!isTRUE(logexp)) {
    logexp <- paste0("log(", exp, ")")
  } else {
    logexp <- exp
  }

  #### translog price index ####
  pmatrix <- matrix(NA, neqs, neqs)

  for (i in seq_len(neqs)) {
    pmatrix[i, ] <- paste0(logp[i], "*", logp)
  }

  translogmatrix <- matrix(NA, neqs, neqs)
  for (i in seq_len(neqs))
    translogmatrix[i, ] <- paste0(gamma[i, ], "*", pmatrix[i, ])

  translog <- paste("(", alpha0, "+", alogp, "+ 0.5 * (",
                    paste0(translogmatrix, collapse = " + "), ")", ")")

  # Stone price index
  if (priceindex == "S") {
    translog <- paste("(", paste(w, "*", logp, collapse = " + "), ")")
  }

  #### cobb douglas price aggregator ####
  cobbdoug <- paste("(", paste(paste("exp(", logp, "*", beta, ")"),
                               collapse = " * "), ")")

  #### combine everything to the final model ####
  eqs <- list()

  if (modeltype == "AI" && !ray) {

    #### eqs ####
    for (i in seq_len(neqs)) {
      glogp <- paste(paste0(gamma[i, ], "*", logp), collapse = " + ")
      eqs[i] <- paste(w[i], "~", alpha[i], "+", glogp, "+", beta[i], "*",
                      "(", logexp, "-", translog, ")")
    }
  }

  if (modeltype == "eAI" && !ray) {

    mu <- paste0("mu_", seq_len(neqs))

    for (i in seq_len(neqs)) {

      # FixMe: this is what Stata indicates. Is it correct?
      eqs[i] <- paste(mu[i], "~", 1, "+ 1 /", w[i], "*",
                      "(", beta[i], ")")
    }
  }

  if (modeltype == "uceAI" && !ray) {

    ue <- sort(rep(paste0("ue_", seq_len(neqs)), neqs))
    ue <- paste0(ue, "_", seq_len(neqs))

    count <- 0

    for (i in seq_len(neqs)) {
      for (j in seq_len(neqs)) {

        count <- count + 1

        glogp <- paste(paste0(gamma[j, ], "*", logp), collapse = " + ")

        delta <- 0; if (i == j) delta <- -1

        eqs[count] <- paste(ue[count], "~", delta,
                            "+ 1 /", w[i], "* (", gamma[i, j], "-",
                            beta[i],
                            "* (", alpha[j], "+", glogp, "))")

      }
    }
  }

  if (modeltype == "ceAI" && !ray) {

    ce <- sort(rep(paste0("ce_", seq_len(neqs)), neqs))
    ce <- paste0(ce, "_", seq_len(neqs))

    count <- 0

    for (i in seq_len(neqs)) {
      for (j in seq_len(neqs)) {

        count <- count + 1

        glogp <- paste(paste0(gamma[j, ], "*", logp), collapse = " + ")

        delta <- 0; if (i == j) delta <- -1

        eqs[count] <- paste(ce[count], "~ (", delta,
                            "+ 1 /", w[i], "* (", gamma[i, j], "-",
                            beta[i],
                            "* (", alpha[j], "+", glogp, "))) + ((",
                            1, "+ 1 /", w[i], "*", beta[i], ") *", w[j], ")")

      }
    }
  }

  if (modeltype == "QAI" && !ray) {

    #### eqs ####
    for (i in seq_len(neqs)) {
      glogp <- paste(paste0(gamma[i, ], "*", logp), collapse = " + ")
      eqs[i] <- paste(w[i], "~", alpha[i], "+", glogp, "+", beta[i], "*",
                      "(", logexp, "-", translog, ") +", lambda[i], "* ((",
                      logexp, "-", translog, ")^2 /", cobbdoug, ")")
    }
  }

  if (modeltype == "eQAI" && !ray) {

    mu <- paste0("mu_", seq_len(neqs))

    for (i in seq_len(neqs)) {

      eqs[i] <- paste(mu[i], "~", 1, "+ 1 /", w[i], "*",
                      "(", beta[i], " +  2 *", lambda[i], "* ((",
                      logexp, "-", translog, ") /",
                      "(", cobbdoug, ")", " ))")
    }
  }

  if (modeltype == "uceQAI" && !ray) {

    ue <- sort(rep(paste0("ue_", seq_len(neqs)), neqs))
    ue <- paste0(ue, "_", seq_len(neqs))

    count <- 0

    for (i in seq_len(neqs)) {
      for (j in seq_len(neqs)) {

        count <- count + 1

        glogp <- paste(paste0(gamma[j, ], "*", logp), collapse = " + ")

        delta <- 0; if (i == j) delta <- -1

        eqs[count] <- paste(ue[count], "~", delta,
                            "+ 1 /", w[i], "* (", gamma[i, j], "- (",
                            beta[i], "+ 2 *", lambda[i], "* ((",
                            logexp, "-", translog, ") /",
                            "(", cobbdoug, ")", " ))",
                            "* (", alpha[j], "+", glogp, ")",
                            "- (", beta[j], "* ", lambda[i],
                            "* ((", logexp, "-", translog, ")^2 /",
                            "(", cobbdoug, ")", " )))")

      }
    }
  }

  if (modeltype == "ceQAI" && !ray) {

    ce <- sort(rep(paste0("ce_", seq_len(neqs)), neqs))
    ce <- paste0(ce, "_", seq_len(neqs))

    count <- 0

    for (i in seq_len(neqs)) {
      for (j in seq_len(neqs)) {

        count <- count + 1

        glogp <- paste(paste0(gamma[j, ], "*", logp), collapse = " + ")

        delta <- 0; if (i == j) delta <- -1

        eqs[count] <- paste(ce[count], "~ (", delta,
                            "+ 1 /", w[i], "* (", gamma[i, j], "- (",
                            beta[i], "+ 2 *", lambda[i], "* ((",
                            logexp, "-", translog, ") /",
                            "(", cobbdoug, ")", " ))",
                            "* (", alpha[j], "+", glogp, ")",
                            "- (", beta[j], "* ", lambda[i],
                            "* ((", logexp, "-", translog, ")^2 /",
                            "(", cobbdoug, ")", " )))) + ((",
                            1, "+ 1 /", w[i], "*",
                            "(", beta[i], " +  2 *", lambda[i], "* ((",
                            logexp, "-", translog, ") /",
                            "(", cobbdoug, ")", " ))) *", w[j], ")")

      }
    }
  }

  if (modeltype == "AI" && ray) {

    translog_star <- paste("(", m0, "+", translog, ")")

    #### eqs ####
    for (i in seq_len(neqs)) {
      glogp <- paste(paste0(gamma[i, ], "*", logp), collapse = " + ")
      eqs[i] <- paste(w[i], "~", alpha[i], "+", glogp, "+", beta_star[i], "*",
                      "(", logexp, "-", translog_star, ")")
    }
  }

  if (modeltype == "eAI" && ray) {

    mu <- paste0("mu_", seq_len(neqs))

    for (i in seq_len(neqs)) {

      # FixMe: this is what Stata indicates. Is it correct?
      eqs[i] <- paste(mu[i], "~", 1, "+ 1 /", w[i], "*",
                      "(", beta_star[i], ")")
    }
  }

  if (modeltype == "uceAI" && ray) {

    ue <- sort(rep(paste0("ue_", seq_len(neqs)), neqs))
    ue <- paste0(ue, "_", seq_len(neqs))

    count <- 0

    for (i in seq_len(neqs)) {
      for (j in seq_len(neqs)) {

        count <- count + 1

        glogp <- paste(paste0(gamma[j, ], "*", logp), collapse = " + ")

        delta <- 0; if (i == j) delta <- -1

        eqs[count] <- paste(ue[count], "~", delta,
                            "+ 1 /", w[i], "* (", gamma[i, j], "-",
                            beta_star[i],
                            "* (", alpha[j], "+", glogp, "))")

      }
    }
  }

  if (modeltype == "ceAI" && ray) {

    ce <- sort(rep(paste0("ce_", seq_len(neqs)), neqs))
    ce <- paste0(ce, "_", seq_len(neqs))

    count <- 0

    for (i in seq_len(neqs)) {
      for (j in seq_len(neqs)) {

        count <- count + 1

        glogp <- paste(paste0(gamma[j, ], "*", logp), collapse = " + ")

        delta <- 0; if (i == j) delta <- -1

        eqs[count] <- paste(ce[count], "~ (", delta,
                            "+ 1 /", w[i], "* (", gamma[i, j], "-",
                            beta_star[i],
                            "* (", alpha[j], "+", glogp, "))) + ((",
                            1, "+ 1 /", w[i], "*", beta_star[i],
                            ") *", w[j], ")")

      }
    }
  }


  if (modeltype == "QAI" && ray) {

    translog_star <- paste("(", m0, "+", translog, ")")

    #### eqs ####
    for (i in seq_len(neqs)) {
      glogp <- paste(paste0(gamma[i, ], "*", logp), collapse = " + ")
      eqs[i] <- paste(w[i], "~", alpha[i], "+", glogp, "+", beta_star[i], "*",
                      "(", logexp, "-", translog_star, ") +", lambda[i], "* ((",
                      logexp, "-", translog_star, ")^2 /",
                      "(", cobbdoug, "*", cpz, ")", " )")
    }
  }

  if (modeltype == "eQAI" && ray) {

    mu <- paste0("mu_", seq_len(neqs))
    translog_star <- paste("(", m0, "+", translog, ")")

    for (i in seq_len(neqs)) {

      eqs[i] <- paste(mu[i], "~", 1, "+ 1 /", w[i], "*",
                      "(", beta_star[i], " +  2 *", lambda[i], "* ((",
                      logexp, "-", translog_star, ") /",
                      "(", cobbdoug, "*", cpz, ")", " ))")
    }
  }

  if (modeltype == "uceQAI" && ray) {

    ue <- sort(rep(paste0("ue_", seq_len(neqs)), neqs))
    ue <- paste0(ue, "_", seq_len(neqs))


    translog_star <- paste("(", m0, "+", translog, ")")

    count <- 0

    for (i in seq_len(neqs)) {
      for (j in seq_len(neqs)) {

        count <- count + 1

        glogp <- paste(paste0(gamma[j, ], "*", logp), collapse = " + ")

        delta <- 0; if (i == j) delta <- -1

        eqs[count] <- paste(ue[count], "~", delta,
                            "+ 1 /", w[i], "* (", gamma[i, j], "- (",
                            beta_star[i], "+ 2 *", lambda[i], "* ((",
                            logexp, "-", translog_star, ") /",
                            "(", cobbdoug, "*", cpz, ")", " ))",
                            "* (", alpha[j], "+", glogp, ")",
                            "- (", beta_star[j], "* ", lambda[i],
                            "* ((", logexp, "-", translog_star, ")^2 /",
                            "(", cobbdoug, "*", cpz, ")", " )))")

      }
    }
  }

  if (modeltype == "ceQAI" && ray) {

    ce <- sort(rep(paste0("ce_", seq_len(neqs)), neqs))
    ce <- paste0(ce, "_", seq_len(neqs))


    translog_star <- paste("(", m0, "+", translog, ")")

    count <- 0

    for (i in seq_len(neqs)) {
      for (j in seq_len(neqs)) {

        count <- count + 1

        glogp <- paste(paste0(gamma[j, ], "*", logp), collapse = " + ")

        delta <- 0; if (i == j) delta <- -1

        eqs[count] <- paste(ce[count], "~ (", delta,
                            "+ 1 /", w[i], "* (", gamma[i, j], "- (",
                            beta_star[i], "+ 2 *", lambda[i], "* ((",
                            logexp, "-", translog_star, ") /",
                            "(", cobbdoug, "*", cpz, ")", " ))",
                            "* (", alpha[j], "+", glogp, ")",
                            "- (", beta_star[j], "* ", lambda[i],
                            "* ((", logexp, "-", translog_star, ")^2 /",
                            "(", cobbdoug, "*", cpz, ")", " )))) + (",
                            1, "+ 1 /", w[i], "*",
                            "(", beta_star[i], " +  2 *", lambda[i], "* ((",
                            logexp, "-", translog_star, ") /",
                            "(", cobbdoug, "*", cpz, ")", " ))) * ", w[j]
        )

      }
    }
  }



  eqs <- unlist(eqs)

  model <- list()

  if (!(modeltype == "eAI" || modeltype == "uceAI" || modeltype == "ceAI")
      && !(modeltype == "eQAI" || modeltype == "uceQAI" || modeltype == "ceQAI")) {

    # drop one equation
    for (i in seq_len(neqs - 1))
      model[[i]] <- as.formula(eqs[i])

  } else {

    for (i in seq_along(eqs))
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

  vars <- c(w, p, x)
  if (missing(z)) {
    z <- substitute()
    scale <- FALSE
  } else {
    vars <- c(vars, z)
  }
  ndat <- names(data)

  if (!all(vars %in% ndat)) {
    stop("Selected Variable not found in Dataset.")
  }

  if (!scale) {
    model <- ai.model(w = w, p = p, exp = x, alph0 = a0, modeltype = "AI",
                      logp = logp, logexp = logexp)
  } else {
    model <- ai.model(w = w, p = p, exp = x, demogr = z, alph0 = a0,
                      ray = TRUE, modeltype = "AI",
                      logp = logp, logexp = logexp)
  }

  res <- nlsur(eqns = model, data = data, type = 3, ...)

  # required for eQAI
  attr(res, "w") <- w
  attr(res, "p") <- p
  attr(res, "x") <- x
  if (!missing(z)) attr(res, "z") <- z
  attr(res, "a0")     <- a0
  attr(res, "logp")   <- logp
  attr(res, "logexp") <- logexp
  attr(res, "scale")  <- scale
  attr(res, "modeltyp") <- "AI"

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

  vars <- c(w, p, x)
  if (missing(z)) {
    z <- substitute()
    scale <- FALSE
  } else {
    vars <- c(vars, z)
  }

  ndat <- names(data)

  if (!all(vars %in% ndat)) {
    stop("Selected Variable not found in Dataset.")
  }

  if (!scale) {
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
  if (!missing(z)) attr(res, "z") <- z
  attr(res, "a0")     <- a0
  attr(res, "logp")   <- logp
  attr(res, "logexp") <- logexp
  attr(res, "scale")  <- scale
  attr(res, "modeltyp") <- "QAI"

  res
}

#' Estimation of elasticities of the (Quadratic) Almost-Ideal Demand System
#'
#' Estimates the income/expenditure elasticity, the uncompensated price
#'  elasticity and the compensated price elasticity
#'
#' @param object qai result
#' @param data data vector used for estimation
#' @param type 1 = expenditure; 2 = uncompensated; 3 = compensated
#' @param usemean evaluate at mean
#'
#' @details Formula for the expenditure (income) elasticity
#'
#' \deqn{
#' \mu_i = 1 + \frac{1}{w_i} \left[ \beta_i + \frac{2\lambda_i}{b(\mathbf{p})} *
#'  ln \left\{\frac{m}{a(\mathbf{p})} \right\}\right]
#' }
#'
#' Formula for the uncompensated price elasticity
#'
#' \deqn{
#' \epsilon_{ij} = \delta_{ij} + \frac{1}{w_i} \left( \gamma_{ij} - \beta_i +
#'  \frac{2\lambda_i}{b(\mathbf{p})} \right) \left[\ln \left\{
#'  \frac{m}{a(\mathbf{p})}\right\} \right] \times \\
#'  \left(\alpha_j + \sum_k \gamma_{jk} \ln p_k \right) -
#'  \frac{\beta_j \lambda_i}{b(\mathbf{p})} \left[
#'   \ln \left\{ \frac{m}{a(\mathbf{p})} \right\}\right]
#' }
#'
#' Compensated price elasticities (Slutsky equation)
#' \deqn{
#'  \epsilon_{ij}^{C} = \epsilon_{ij} + \mu_i w_j
#' }
#'
#' @examples
#' \dontrun{
#' library(nlsur)
#' library(readstata13)
#'
#' dd <- read.dta13("http://www.stata-press.com/data/r15/food.dta")
#'
#' w <- paste0("w", 1:4); p <- paste0("p", 1:4); x <- "expfd"
#'
#' est <- ai(w = w, p = p, x = x, data = dd, a0 = 10, scale = FALSE,
#'           logp = F, logexp = F)
#'
#' mu <- elasticities(est, data = dd, type = 1, usemean = FALSE)
#'
#' ue <- elasticities(est, data = dd, type = 2, usemean = FALSE)
#'
#' ce <- elasticities(est, data = dd, type = 3, usemean = FALSE)}
#'
#' @references Banks, James, Blundell, Richard, Lewbel, Arthur: Quadratic Engel
#'  Curves and Consumer Demand, The Review of Economics and Statistics 79(4),
#'  The MIT Press, 527-539, 1997
#' @references Poi, Brian P.: Easy demand-system estimation with quaids, The
#'  Stata Journal 12(3), 433-446, 2012
#'
#' @seealso ai and qai
#'
#' @export
elasticities <- function(object, data, type = 1, usemean = FALSE) {


  modeltyp  <- attr(object, "modeltyp")

  if (modeltyp == "AI") {
    if (type == 1)
      mtyp <- "eAI"

    if (type == 2)
      mtyp <- "uceAI"

    if (type == 3)
      mtyp <- "ceAI"
  }

  if (modeltyp == "QAI") {
    if (type == 1)
      mtyp <- "eQAI"

    if (type == 2)
      mtyp <- "uceQAI"

    if (type == 3)
      mtyp <- "ceQAI"
  }


  # extract var names
  w <- attr(object, "w")
  p <- attr(object, "p")
  x <- attr(object, "x")
  z <- attr(object, "z")
  a0 <- attr(object, "a0")
  logp   <- attr(object, "logp")
  logexp <- attr(object, "logexp")
  scale  <- attr(object, "scale")

  if (!scale) {
    eqs <- ai.model(w = w, p = p, exp = x, modeltype = mtyp,
                    logp = logp, logexp = logexp, alph0 = a0)
  } else {
    eqs <- ai.model(w = w, p = p, exp = x, demogr = z,
                    ray = TRUE, modeltype = mtyp,
                    logp = logp, logexp = logexp, alph0 = a0)
  }

  if (!usemean) {

    eqns_lhs <- lapply(X = eqs, FUN = function(x) x[[2L]])
    eqns_rhs <- lapply(X = eqs, FUN = function(x) x[[3L]])
    vnam     <- sapply(X = eqns_lhs, FUN = as.character)

    data2 <- data.frame(data, as.list(coef(object)))

    # create fit: predict result
    fit <- lapply(X = eqns_rhs, FUN = eval, envir = data2)

    # replace is.infinite() with NA
    fit <- lapply(X = fit, FUN = function(x) replace(x, is.infinite(x), NA))

    fit <- data.frame(fit)
    names(fit) <- vnam

  } else {

    fit <- NULL

    # log means are required
    logmean <- function(x) exp(mean(log(x)))

    # calculate means
    wm <- sapply(w, FUN = function(x) mean(data[[x]], use.na = FALSE))
    pm <- sapply(p, FUN = function(x) mean(data[[x]], use.na = FALSE))
    xm <- sapply(x, FUN = function(x) mean(data[[x]], use.na = FALSE))

    if (!logp) pm <- sapply(p, FUN = function(x) logmean(data[[x]]))
    if (!logexp) xm <- sapply(x, FUN = function(x) logmean(data[[x]]))

    ms <- data.frame(t(c(wm, pm, xm)))

    if (scale) {
      zm <- sapply(z, FUN = function(x) mean(data[[x]], use.na = FALSE))
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

    fit <- as.data.frame(fit)
  }

  fit
}
