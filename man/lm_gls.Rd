% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlsur.R
\name{lm_gls}
\alias{lm_gls}
\title{Calculate WLS using sparse matrix and qr}
\usage{
lm_gls(X, Y, W, neqs, tol = 1e-07, covb = FALSE)
}
\arguments{
\item{X}{n x m X matrix}

\item{Y}{n x k matrix}

\item{W}{n x n}

\item{neqs}{k}

\item{tol}{tolerance for qr}

\item{covb}{if true covb is calculated else theta}
}
\description{
calculate WLS using eigen similar to the approach in
MASS::lm.gls
}
