% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ssr_est}
\alias{ssr_est}
\title{Estimate residual sum of squares}
\usage{
ssr_est(r, s, w)
}
\arguments{
\item{r}{residuals}

\item{s}{weighting matrix}

\item{w}{vector of weights}
}
\description{
calculate SSR where
\eqn{SSR(\beta) = u'D'Du.}
}
