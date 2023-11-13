// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rdirichlet_cpp
arma::vec rdirichlet_cpp(arma::vec alpha);
RcppExport SEXP _baysc_rdirichlet_cpp(SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(rdirichlet_cpp(alpha));
    return rcpp_result_gen;
END_RCPP
}
// rcat_cpp
int rcat_cpp(arma::rowvec probs);
RcppExport SEXP _baysc_rcat_cpp(SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcat_cpp(probs));
    return rcpp_result_gen;
END_RCPP
}
// mvrnorm_cpp
arma::mat mvrnorm_cpp(const int& n, const arma::vec& mu, const arma::mat& sigma);
RcppExport SEXP _baysc_mvrnorm_cpp(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnorm_cpp(n, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// rtruncnorm_cpp
double rtruncnorm_cpp(const int& n, const double& a, const double& b, const double& mean, const double& sd);
RcppExport SEXP _baysc_rtruncnorm_cpp(SEXP nSEXP, SEXP aSEXP, SEXP bSEXP, SEXP meanSEXP, SEXP sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const double& >::type sd(sdSEXP);
    rcpp_result_gen = Rcpp::wrap(rtruncnorm_cpp(n, a, b, mean, sd));
    return rcpp_result_gen;
END_RCPP
}
// logSumExp_cpp
double logSumExp_cpp(const arma::rowvec& x);
RcppExport SEXP _baysc_logSumExp_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logSumExp_cpp(x));
    return rcpp_result_gen;
END_RCPP
}
// update_pi
void update_pi(arma::vec& pi, const arma::vec& w_all, const arma::vec& c_all, const int& K, const arma::vec& alpha);
RcppExport SEXP _baysc_update_pi(SEXP piSEXP, SEXP w_allSEXP, SEXP c_allSEXP, SEXP KSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w_all(w_allSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type c_all(c_allSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type alpha(alphaSEXP);
    update_pi(pi, w_all, c_all, K, alpha);
    return R_NilValue;
END_RCPP
}
// update_c
void update_c(arma::vec& c_all, const int& n, const int& K, const int& J, const arma::cube& theta, const arma::mat& x_mat, const arma::vec& pi, const arma::vec& z_all, const arma::mat& V, const arma::mat& xi, const arma::vec& y_all);
RcppExport SEXP _baysc_update_c(SEXP c_allSEXP, SEXP nSEXP, SEXP KSEXP, SEXP JSEXP, SEXP thetaSEXP, SEXP x_matSEXP, SEXP piSEXP, SEXP z_allSEXP, SEXP VSEXP, SEXP xiSEXP, SEXP y_allSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type c_all(c_allSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int& >::type J(JSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x_mat(x_matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type z_all(z_allSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y_all(y_allSEXP);
    update_c(c_all, n, K, J, theta, x_mat, pi, z_all, V, xi, y_all);
    return R_NilValue;
END_RCPP
}
// update_c_wolca
void update_c_wolca(arma::vec& c_all, const int& n, const int& K, const int& J, const arma::cube& theta, const arma::mat& x_mat, const arma::vec& pi);
RcppExport SEXP _baysc_update_c_wolca(SEXP c_allSEXP, SEXP nSEXP, SEXP KSEXP, SEXP JSEXP, SEXP thetaSEXP, SEXP x_matSEXP, SEXP piSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type c_all(c_allSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int& >::type J(JSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x_mat(x_matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pi(piSEXP);
    update_c_wolca(c_all, n, K, J, theta, x_mat, pi);
    return R_NilValue;
END_RCPP
}
// update_theta
void update_theta(arma::cube& theta, const int& J, const int& K, const int& R, const arma::mat& eta, const arma::vec& w_all, const arma::vec& c_all, arma::mat& x_mat);
RcppExport SEXP _baysc_update_theta(SEXP thetaSEXP, SEXP JSEXP, SEXP KSEXP, SEXP RSEXP, SEXP etaSEXP, SEXP w_allSEXP, SEXP c_allSEXP, SEXP x_matSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const int& >::type J(JSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w_all(w_allSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type c_all(c_allSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x_mat(x_matSEXP);
    update_theta(theta, J, K, R, eta, w_all, c_all, x_mat);
    return R_NilValue;
END_RCPP
}
// update_xi
arma::mat update_xi(arma::mat& xi, const int& n, const int& K, const arma::vec& w_all, const arma::vec& c_all, const arma::vec& z_all, const arma::mat& V, const List& mu0, const List& Sig0);
RcppExport SEXP _baysc_update_xi(SEXP xiSEXP, SEXP nSEXP, SEXP KSEXP, SEXP w_allSEXP, SEXP c_allSEXP, SEXP z_allSEXP, SEXP VSEXP, SEXP mu0SEXP, SEXP Sig0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w_all(w_allSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type c_all(c_allSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type z_all(z_allSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const List& >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< const List& >::type Sig0(Sig0SEXP);
    rcpp_result_gen = Rcpp::wrap(update_xi(xi, n, K, w_all, c_all, z_all, V, mu0, Sig0));
    return rcpp_result_gen;
END_RCPP
}
// update_z
arma::vec update_z(arma::vec& z_all, const int& n, const arma::mat& V, const arma::mat& xi, const arma::vec& c_all, const arma::vec& y_all);
RcppExport SEXP _baysc_update_z(SEXP z_allSEXP, SEXP nSEXP, SEXP VSEXP, SEXP xiSEXP, SEXP c_allSEXP, SEXP y_allSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type z_all(z_allSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type c_all(c_allSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y_all(y_allSEXP);
    rcpp_result_gen = Rcpp::wrap(update_z(z_all, n, V, xi, c_all, y_all));
    return rcpp_result_gen;
END_RCPP
}
// update_loglik
void update_loglik(arma::vec& loglik, const int& n, const int& J, const arma::vec& c_all, const arma::cube& theta, const arma::mat& x_mat, const arma::vec& pi, const arma::vec& z_all, const arma::mat& V, const arma::mat& xi, const arma::vec& y_all);
RcppExport SEXP _baysc_update_loglik(SEXP loglikSEXP, SEXP nSEXP, SEXP JSEXP, SEXP c_allSEXP, SEXP thetaSEXP, SEXP x_matSEXP, SEXP piSEXP, SEXP z_allSEXP, SEXP VSEXP, SEXP xiSEXP, SEXP y_allSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type loglik(loglikSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int& >::type J(JSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type c_all(c_allSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x_mat(x_matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type z_all(z_allSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y_all(y_allSEXP);
    update_loglik(loglik, n, J, c_all, theta, x_mat, pi, z_all, V, xi, y_all);
    return R_NilValue;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_stan_fit4SWOLCA_main_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_baysc_rdirichlet_cpp", (DL_FUNC) &_baysc_rdirichlet_cpp, 1},
    {"_baysc_rcat_cpp", (DL_FUNC) &_baysc_rcat_cpp, 1},
    {"_baysc_mvrnorm_cpp", (DL_FUNC) &_baysc_mvrnorm_cpp, 3},
    {"_baysc_rtruncnorm_cpp", (DL_FUNC) &_baysc_rtruncnorm_cpp, 5},
    {"_baysc_logSumExp_cpp", (DL_FUNC) &_baysc_logSumExp_cpp, 1},
    {"_baysc_update_pi", (DL_FUNC) &_baysc_update_pi, 5},
    {"_baysc_update_c", (DL_FUNC) &_baysc_update_c, 11},
    {"_baysc_update_c_wolca", (DL_FUNC) &_baysc_update_c_wolca, 7},
    {"_baysc_update_theta", (DL_FUNC) &_baysc_update_theta, 8},
    {"_baysc_update_xi", (DL_FUNC) &_baysc_update_xi, 9},
    {"_baysc_update_z", (DL_FUNC) &_baysc_update_z, 6},
    {"_baysc_update_loglik", (DL_FUNC) &_baysc_update_loglik, 11},
    {"_rcpp_module_boot_stan_fit4SWOLCA_main_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4SWOLCA_main_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_baysc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
