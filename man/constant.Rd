% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{constant}
\alias{constant}
\title{Check if formula contains constant}
\usage{
constant(x)
}
\arguments{
\item{x}{formula}
}
\description{
Check if formula contains constant
}
\details{
Primitive function to check a formula for a constant part.
Function checks first and last term on rhs for constant variables
at front and back position.
}
\examples{
\dontrun{
constant(y ~ x + a * z) # x
constant(y ~ x * b + 1) # 1
constant(y  ~ 0 + x) # NULL
constant(y  ~ x) # x
constant(y ~ x1 * b1 + b0 + x2 * b2) # wont find b0
constant( y  ~ (x*b +k) + a*y + b*z ) # wont find k
constant( y  ~ (k+ x*b) + a*y + b*z ) # k
constant( y  ~  a*y + b*z + (k + x*b) ) # wont find k
constant( y  ~  a*y + b*z + (x*b + k) ) # k
}
}
