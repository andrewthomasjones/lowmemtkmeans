// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cluster_BIC
double cluster_BIC(arma::mat& data, arma::mat& centres);
RcppExport SEXP tkmeans_cluster_BIC(SEXP dataSEXP, SEXP centresSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type centres(centresSEXP);
    __result = Rcpp::wrap(cluster_BIC(data, centres));
    return __result;
END_RCPP
}
// tkmeans
arma::mat tkmeans(arma::mat& M, int k, double alpha, int nstart, int iter, double tol, bool verbose);
RcppExport SEXP tkmeans_tkmeans(SEXP MSEXP, SEXP kSEXP, SEXP alphaSEXP, SEXP nstartSEXP, SEXP iterSEXP, SEXP tolSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type nstart(nstartSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    __result = Rcpp::wrap(tkmeans(M, k, alpha, nstart, iter, tol, verbose));
    return __result;
END_RCPP
}
// scale_mat_inplace
arma::mat scale_mat_inplace(arma::mat& M);
RcppExport SEXP tkmeans_scale_mat_inplace(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat& >::type M(MSEXP);
    __result = Rcpp::wrap(scale_mat_inplace(M));
    return __result;
END_RCPP
}
// nearest_cluster
arma::uvec nearest_cluster(arma::mat& data, arma::mat& centres);
RcppExport SEXP tkmeans_nearest_cluster(SEXP dataSEXP, SEXP centresSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type centres(centresSEXP);
    __result = Rcpp::wrap(nearest_cluster(data, centres));
    return __result;
END_RCPP
}
