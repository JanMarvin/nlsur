#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP calc_ssr (arma::Mat<double> r, arma::Mat<double> s, Rcpp::List eqs) {

  arma::mat ssr(1,1, fill::zeros);
  int neqs = eqs.size();
  int n    = r.n_rows;

  s = s.t();

  for (int j = 0; j < neqs; ++j) {
    for (int i = 0; i < n; ++i){
      ssr += pow( r.row(i) * s.col(j), 2);
    }
  }

  return wrap(ssr);
}

// [[Rcpp::export]]
arma::Mat<double> arma_reshape(arma::Mat<double> mm, int sizetheta) {

  arma::vec v = vectorise(mm);

  int newsize = v.n_elem / sizetheta;
  mm.set_size(newsize, sizetheta);
  mm = mm.t();

  int k = 0;
  for (uint j = 0; j < mm.n_cols; ++j) {
    for (uint i = 0; i < mm.n_rows; ++i) {
      mm(i,j) = v(k);
      k += +1;
    }
  }

  return mm.t();
}

// [[Rcpp::export]]
SEXP calc_reg (arma::Mat<double> x, arma::Mat<double> r, arma::Mat<double> qS,
               int sizetheta, int neqs) {

  arma::Mat<double> XDX(sizetheta, sizetheta, fill::zeros);
  arma::Mat<double> XDy(sizetheta, 1, fill::zeros);

  // arma::Mat<double> qS = pinv(S);

  int n = r.n_rows;
  int xicol = x.n_cols / neqs;

  for (int i = 0; i < n; ++i) {
    arma::Mat<double> xi = x.row(i);
    arma::Mat<double> XI = arma_reshape(xi, xicol);

    arma::Mat<double> YI = r.row(i).t();

    XDX += XI.t() * qS * XI;
    XDy += XI.t() * qS * YI;
  }

  XDX = 0.5 * ( XDX + XDX.t() );

  return wrap(solve(XDX, XDy).t());
}
