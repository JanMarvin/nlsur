#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP calc_ssr (arma::mat  r, arma::mat s, Rcpp::List eqs) {

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
arma::mat arma_reshape(arma::mat mm, int sizetheta) {

  // arma::mat mm = m.row(0);
  arma::vec v = vectorise(mm);

  int newsize = v.n_elem / sizetheta;
  mm.set_size(newsize, sizetheta);
  mm = mm.t();

  int k = 0;
  for (int j = 0; j < mm.n_cols; ++j) {
    for (int i = 0; i < mm.n_rows; ++i) {
      mm(i,j) = v(k);
      k += +1;
    }
  }

  return mm.t();
}

// [[Rcpp::export]]
SEXP calc_reg (arma::mat x, arma::mat r, arma::mat qS,
               int sizetheta, int neqs) {

  arma::mat XDX(sizetheta, sizetheta, fill::zeros);
  arma::mat XDy(sizetheta, 1, fill::zeros);

  int n = r.n_rows;
  int xicol = x.n_cols / neqs;

  Rprintf("Anzahl an Spalten: %d \n", x.n_cols);
  Rprintf("Anzahl an Gleichungen: %d \n", neqs);
  Rprintf("Neue Matrix hat dim: %d %d \n", neqs, xicol);
  Rprintf("n: %d \n", n);

  for (int i = 0; i < n; ++i) {
    arma::mat xi = x.row(i);
    arma::mat XI = arma_reshape(xi, xicol);

    arma::mat yi = r.row(i);
    arma::mat YI = arma_reshape(yi, 1);

    XDX += XI.t() * qS * XI;
    XDy += XI.t() * qS * YI;
  }
  XDX = 0.5 * ( XDX + XDX.t() );

  // return wrap(solve(XDX, XDy));
  return wrap(XDX);
  // return wrap(XDy);
}
