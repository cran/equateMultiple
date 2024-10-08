// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// objectivefzRcpp
double objectivefzRcpp(NumericVector par, int T, NumericVector ab, NumericVector wt, NumericMatrix aj1T, NumericMatrix bj1T, NumericMatrix cj1T, int nummet, int itmp, double D, int base);
RcppExport SEXP _equateMultiple_objectivefzRcpp(SEXP parSEXP, SEXP TSEXP, SEXP abSEXP, SEXP wtSEXP, SEXP aj1TSEXP, SEXP bj1TSEXP, SEXP cj1TSEXP, SEXP nummetSEXP, SEXP itmpSEXP, SEXP DSEXP, SEXP baseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ab(abSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type aj1T(aj1TSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bj1T(bj1TSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type cj1T(cj1TSEXP);
    Rcpp::traits::input_parameter< int >::type nummet(nummetSEXP);
    Rcpp::traits::input_parameter< int >::type itmp(itmpSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type base(baseSEXP);
    rcpp_result_gen = Rcpp::wrap(objectivefzRcpp(par, T, ab, wt, aj1T, bj1T, cj1T, nummet, itmp, D, base));
    return rcpp_result_gen;
END_RCPP
}
// hessRcpp
NumericMatrix hessRcpp(NumericVector par, int T, NumericVector ab, NumericVector wt, NumericMatrix aj1T, NumericMatrix bj1T, NumericMatrix cj1T, int nummet, int itmp, double D, int base);
RcppExport SEXP _equateMultiple_hessRcpp(SEXP parSEXP, SEXP TSEXP, SEXP abSEXP, SEXP wtSEXP, SEXP aj1TSEXP, SEXP bj1TSEXP, SEXP cj1TSEXP, SEXP nummetSEXP, SEXP itmpSEXP, SEXP DSEXP, SEXP baseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ab(abSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type aj1T(aj1TSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bj1T(bj1TSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type cj1T(cj1TSEXP);
    Rcpp::traits::input_parameter< int >::type nummet(nummetSEXP);
    Rcpp::traits::input_parameter< int >::type itmp(itmpSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type base(baseSEXP);
    rcpp_result_gen = Rcpp::wrap(hessRcpp(par, T, ab, wt, aj1T, bj1T, cj1T, nummet, itmp, D, base));
    return rcpp_result_gen;
END_RCPP
}
// gradRcpp
NumericVector gradRcpp(NumericVector par, int T, NumericVector ab, NumericVector wt, NumericMatrix aj1T, NumericMatrix bj1T, NumericMatrix cj1T, int nummet, int itmp, double D, int base);
RcppExport SEXP _equateMultiple_gradRcpp(SEXP parSEXP, SEXP TSEXP, SEXP abSEXP, SEXP wtSEXP, SEXP aj1TSEXP, SEXP bj1TSEXP, SEXP cj1TSEXP, SEXP nummetSEXP, SEXP itmpSEXP, SEXP DSEXP, SEXP baseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ab(abSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type aj1T(aj1TSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bj1T(bj1TSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type cj1T(cj1TSEXP);
    Rcpp::traits::input_parameter< int >::type nummet(nummetSEXP);
    Rcpp::traits::input_parameter< int >::type itmp(itmpSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type base(baseSEXP);
    rcpp_result_gen = Rcpp::wrap(gradRcpp(par, T, ab, wt, aj1T, bj1T, cj1T, nummet, itmp, D, base));
    return rcpp_result_gen;
END_RCPP
}
// partialABgammaRcpp
List partialABgammaRcpp(NumericVector par, int T, NumericVector ab, NumericVector wt, NumericMatrix aj1T, NumericMatrix bj1T, NumericMatrix cj1T, int nummet, int itmp, double D, int base, int nb, IntegerVector posnomi);
RcppExport SEXP _equateMultiple_partialABgammaRcpp(SEXP parSEXP, SEXP TSEXP, SEXP abSEXP, SEXP wtSEXP, SEXP aj1TSEXP, SEXP bj1TSEXP, SEXP cj1TSEXP, SEXP nummetSEXP, SEXP itmpSEXP, SEXP DSEXP, SEXP baseSEXP, SEXP nbSEXP, SEXP posnomiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ab(abSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type aj1T(aj1TSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type bj1T(bj1TSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type cj1T(cj1TSEXP);
    Rcpp::traits::input_parameter< int >::type nummet(nummetSEXP);
    Rcpp::traits::input_parameter< int >::type itmp(itmpSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type base(baseSEXP);
    Rcpp::traits::input_parameter< int >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type posnomi(posnomiSEXP);
    rcpp_result_gen = Rcpp::wrap(partialABgammaRcpp(par, T, ab, wt, aj1T, bj1T, cj1T, nummet, itmp, D, base, nb, posnomi));
    return rcpp_result_gen;
END_RCPP
}
// ipfRcpp
List ipfRcpp(NumericMatrix aj1T, int base, double eps);
RcppExport SEXP _equateMultiple_ipfRcpp(SEXP aj1TSEXP, SEXP baseSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type aj1T(aj1TSEXP);
    Rcpp::traits::input_parameter< int >::type base(baseSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(ipfRcpp(aj1T, base, eps));
    return rcpp_result_gen;
END_RCPP
}
// profLikRcpp
double profLikRcpp(arma::vec par, arma::vec coef, arma::uvec t, List X_list, List itmvar, int numforms, arma::uvec notbase, arma::uvec DffcltNum, arma::uvec DscrmnNum, arma::uvec pos);
RcppExport SEXP _equateMultiple_profLikRcpp(SEXP parSEXP, SEXP coefSEXP, SEXP tSEXP, SEXP X_listSEXP, SEXP itmvarSEXP, SEXP numformsSEXP, SEXP notbaseSEXP, SEXP DffcltNumSEXP, SEXP DscrmnNumSEXP, SEXP posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type coef(coefSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type t(tSEXP);
    Rcpp::traits::input_parameter< List >::type X_list(X_listSEXP);
    Rcpp::traits::input_parameter< List >::type itmvar(itmvarSEXP);
    Rcpp::traits::input_parameter< int >::type numforms(numformsSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type notbase(notbaseSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type DffcltNum(DffcltNumSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type DscrmnNum(DscrmnNumSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type pos(posSEXP);
    rcpp_result_gen = Rcpp::wrap(profLikRcpp(par, coef, t, X_list, itmvar, numforms, notbase, DffcltNum, DscrmnNum, pos));
    return rcpp_result_gen;
END_RCPP
}
// profLikRcpp_1PL
double profLikRcpp_1PL(arma::vec par, arma::vec coef, arma::uvec t, List X_list, List itmvar, int numforms, arma::uvec notbase, arma::uvec pos);
RcppExport SEXP _equateMultiple_profLikRcpp_1PL(SEXP parSEXP, SEXP coefSEXP, SEXP tSEXP, SEXP X_listSEXP, SEXP itmvarSEXP, SEXP numformsSEXP, SEXP notbaseSEXP, SEXP posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type coef(coefSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type t(tSEXP);
    Rcpp::traits::input_parameter< List >::type X_list(X_listSEXP);
    Rcpp::traits::input_parameter< List >::type itmvar(itmvarSEXP);
    Rcpp::traits::input_parameter< int >::type numforms(numformsSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type notbase(notbaseSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type pos(posSEXP);
    rcpp_result_gen = Rcpp::wrap(profLikRcpp_1PL(par, coef, t, X_list, itmvar, numforms, notbase, pos));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_equateMultiple_objectivefzRcpp", (DL_FUNC) &_equateMultiple_objectivefzRcpp, 11},
    {"_equateMultiple_hessRcpp", (DL_FUNC) &_equateMultiple_hessRcpp, 11},
    {"_equateMultiple_gradRcpp", (DL_FUNC) &_equateMultiple_gradRcpp, 11},
    {"_equateMultiple_partialABgammaRcpp", (DL_FUNC) &_equateMultiple_partialABgammaRcpp, 13},
    {"_equateMultiple_ipfRcpp", (DL_FUNC) &_equateMultiple_ipfRcpp, 3},
    {"_equateMultiple_profLikRcpp", (DL_FUNC) &_equateMultiple_profLikRcpp, 10},
    {"_equateMultiple_profLikRcpp_1PL", (DL_FUNC) &_equateMultiple_profLikRcpp_1PL, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_equateMultiple(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
