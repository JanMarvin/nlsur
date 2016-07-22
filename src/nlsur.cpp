#include <RcppArmadillo.h>


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//' calc_ssr
//' @param r r
//' @param s s
//' @param w w
//' @export
// [[Rcpp::export]]
SEXP calc_ssr (arma::Mat<double> r, arma::Mat<double> s, arma::Col<double> w) {

  arma::mat ssr(1,1, fill::zeros);
  int n    = r.n_rows;
  int k    = r.n_cols;

  // n / sum(w)
  // only w contains information about the size of n
  double scale = w.n_elem / sum(w);

  s = s.t();

  for (int j = 0; j < k; ++j) {
    for (int i = 0; i < n; ++i){
      ssr += w(i) * pow( r.row(i) * s.col(j), 2);
    }
  }

  return wrap(ssr * scale);
}

//' arma_reshape
//' @param mm mm
//' @param sizetheta sizetheta
//' @export
// [[Rcpp::export]]
arma::Mat<double> arma_reshape(arma::Mat<double> mm, int sizetheta) {

  mm = mm.t();

  int newsize = mm.n_elem / sizetheta;
  mm.set_size(newsize, sizetheta);

  return mm.t();
}

//' calc_reg
//' @param x x
//' @param r r
//' @param qS qS
//' @param w w
//' @param sizetheta sizetheta
//' @param fullreg fullreg
//' @param tol tol
//' @export
// [[Rcpp::export]]
SEXP calc_reg (arma::Mat<double> x, arma::Mat<double> r, arma::Mat<double> qS,
               arma::Col<double> w, int sizetheta, bool fullreg, double tol) {

  arma::Mat<double> XDX(sizetheta, sizetheta, fill::zeros);
  arma::Mat<double> XDy(sizetheta, 1, fill::zeros);

  Function Rf_qr("qr");
  Function Rf_qrcoef("qr.coef");

  int n = r.n_rows;
  int k = r.n_cols;

  for (int i = 0; i < n; ++i) {

    arma::Mat<double> XI = arma_reshape(x.row(i), k);
    XDX += w(i) * XI.t() * qS * XI;

    if (fullreg) {
      arma::Mat<double> YI = r.row(i).t();
      XDy += w(i) * XI.t() * qS * YI;
    }

  }
  XDX = 0.5 * ( XDX + XDX.t() );

  if (fullreg) /* weighted regression */
    return Rf_qrcoef(Rf_qr(XDX, _["tol"] = tol), XDy);
  else         /* covb */
    return wrap(XDX);
}


//' wt_mean
//' @param x x
//' @param w w
//' @export
// [[Rcpp::export]]
SEXP wt_mean(arma::Col<double>& x, arma::Col<double>& w) {

  return( wrap(sum(w % x) / sum(w)) );

}


//' calc_robust
//' @param x x
//' @param u u
//' @param qS qS
//' @param w w
//' @param sizetheta sizetheta
//' @export
// [[Rcpp::export]]
SEXP calc_robust (arma::Mat<double> x, arma::Mat<double> u, arma::Mat<double> qS,
                  arma::Col<double> w, int sizetheta) {

  arma::Mat<double> XDuuDX(sizetheta, sizetheta, fill::zeros);

  int n = u.n_rows;
  int k = u.n_cols;

  for (int i = 0; i < n; ++i) {

    arma::Mat<double> XI = arma_reshape(x.row(i), k);

    arma::Mat<double> UI = u.row(i).t();

    XDuuDX += w(i) * XI.t() * qS * UI * UI.t() * qS * XI;
  }

  XDuuDX = 0.5 * ( XDuuDX + XDuuDX.t() );


  return wrap(XDuuDX);
}
