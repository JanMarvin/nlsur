#include <RcppArmadillo.h>


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//' Estimate residual sum of squares
//' @description calculate SSR where
//' \eqn{SSR(\beta) = u'D'Du.}
//' @param r residuals
//' @param s weighting matrix
//' @param w vector of weights
//' @export
// [[Rcpp::export]]
SEXP ssr_est(arma::Mat<double> r, arma::Mat<double> s, arma::Mat<double> w) {

  arma::mat ssr(1,1, fill::zeros);
  int n = r.n_rows, k = r.n_cols;

  // n / sum(w) : only w contains information about the size of n
  arma::rowvec scale = w.n_rows / arma::sum(w, 0);

  // s is transposed to avoid incompatible matrix dimensions
  for (int j = 0; j < k; ++j) {
    for (int i = 0; i < n; ++i) {
      ssr += (w(i,j) * pow( r.row(i) * s.row(j).t(), 2)) * scale(j);
    }

    Rcpp::checkUserInterrupt();
  }

  return wrap(ssr);
}

//' Reshape matrix for blockwise WLS estimation
//' @param mm a matrix
//' @param sizetheta integer of length(theta) to shrink mm into
//' @description reshape mm for blockwise multiplication in wls_est
//' @examples
//' mm <- matrix(c(11,21,31,41,
//'  12,22,32,42,
//'  13,23,33,43,
//'  14,24,34,44),
//'  ncol = 4)
//'
//' mm_a <- arma_reshape(mm, 2)
//'
//' mm_m <- matrix(t(mm), nrow = 2, byrow = TRUE)
//' @export
// [[Rcpp::export]]
arma::Mat<double> arma_reshape(arma::Mat<double> mm, int sizetheta) {

  mm = mm.t();

  int newsize = mm.n_elem / sizetheta;
  mm.set_size(newsize, sizetheta);

  return mm.t();
}

//' Blockwise WLS estimation
//' @description
//' Blockwise WLS estimation. Usually for \eqn{(X'X)^{-1} W^{-1} X'Y} X and Y
//' X, W and Y are of similar dimensions. In nlsur W is a cov-matrix of size
//' \eqn{k \times k} and usually way smaller than X. To avoid blowing all
//' matrices up for the estimation, a blockwise approach is used. X is shrunken
//' to match size k. W is D'D so XDX is calculated. XDy is only calculated if
//' wanted for a full WLS. For the cov-matrix only XDX is required.
//' @param x matrix of derivatives
//' @param r residual matrix
//' @param qS weighting matrix of sizetheta x sizetheta
//' @param w vector of weights
//' @param sizetheta integer defining the amount of coefficients
//' @param fullreg bool defining if WLS or Cov is calculated
//' @param tol tolerance used for qr()
//' @details as reference see:
//' http://www.navipedia.net/index.php/Block-Wise_Weighted_Least_Square
//' @export
// [[Rcpp::export]]
SEXP wls_est(arma::Mat<double> x, arma::Mat<double> r, arma::Mat<double> qS,
             arma::Mat<double> w, int sizetheta, bool fullreg, double tol) {

  arma::Mat<double> XDX(sizetheta, sizetheta, fill::zeros);
  arma::Mat<double> XDy(sizetheta, 1, fill::zeros);

  Function Rf_qr("qr");
  Function Rf_qrcoef("qr.coef");

  int n = r.n_rows, k = r.n_cols;

  for (int i = 0; i < n; ++i) {

    arma::Mat<double> wts = w.row(i) % qS.each_row();
    // Rcpp::Rcout<< wts << std::endl;

    arma::Mat<double> XI = arma_reshape(x.row(i), k);
    XDX +=   (XI.t() * wts * XI);

    if (fullreg) {
      arma::Mat<double> YI = r.row(i).t();
      XDy += (XI.t() * wts * YI);
    }

    Rcpp::checkUserInterrupt();
  }

  // force symetry on the matrix
  XDX = 0.5 * (XDX + XDX.t());

  if (fullreg) /* weighted regression */
    return Rf_qrcoef(Rf_qr(XDX, _["tol"] = tol), XDy);
  else         /* covb */
    return wrap(XDX);
}

//' Calculate a weighted mean
//' @param x matrix of derivatives
//' @param w vector of weights
//' @export
// [[Rcpp::export]]
SEXP wt_mean(arma::Col<double>& x, arma::Col<double>& w) {

  return( wrap(sum(w % x) / sum(w)) );

}

//' Calculate a robust covariance matrix
//'
//' @param x matrix of derivatives
//' @param u u
//' @param qS weighting matrix
//' @param w vector of weights
//' @param sizetheta sizetheta
//' @description As discussed in Wooldridge (2002, 160)
//' @export
// [[Rcpp::export]]
SEXP cov_robust(arma::Mat<double> x, arma::Mat<double> u, arma::Mat<double> qS,
                arma::Mat<double> w, int sizetheta) {

  arma::Mat<double> XDuuDX(sizetheta, sizetheta, fill::zeros);

  int n = u.n_rows, k = u.n_cols;

  for (int i = 0; i < n; ++i) {
    arma::Mat<double> wts = w.row(i) % qS.each_row();

    arma::Mat<double> XI = arma_reshape(x.row(i), k);

    arma::Mat<double> UI = u.row(i).t();

    // prev was w(i) * XI.t()
    // to avoid dual weighting, only weight the first matrix. Correct?
    XDuuDX += XI.t() * wts * UI * UI.t() * qS * XI;
  }

  // force symetry on the matrix
  XDuuDX = 0.5 * ( XDuuDX + XDuuDX.t() );

  return wrap(XDuuDX);
}
