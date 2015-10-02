#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP calc_ssr (arma::mat  r,  arma::mat s, Rcpp::List eqs) {

  arma::mat ssr(1,1);
  int neqs = eqs.size();
  int n    = r.n_rows;

  for (int j = 0; j < neqs; ++j) {
    for (int i = 0; i < n; ++i){
      ssr += pow( r.row(i) * s.col(j), 2);
    }
  }

  return wrap(ssr);
}

// [[Rcpp::export]]
SEXP calc_reg (arma::mat x, arma::mat r, arma::mat qS, int sizetheta,
               Rcpp::List eqns) {

  arma::mat XDX(sizetheta, sizetheta);
  arma::mat XDy(sizetheta, 1);

  x = x.t();

  int neqs = eqns.size();
  int n = r.n_rows;

  for (int i = 0; i < n; ++i) {
    arma::mat XI = x.col(i);
    XI.reshape(neqs, sizetheta, 1);
    arma::mat yi = r.row(i);
    yi.reshape(yi.n_cols, yi.n_rows);

    XDX += XI.t() * qS * XI;
    XDy += XI.t() * qS * yi;
  }

  XDX = 0.5 * ( XDX + XDX.t() );

  // return wrap(solve(XDX, XDy));
  return wrap(XDy);
}
