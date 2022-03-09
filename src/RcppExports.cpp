// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcpp_sum
double rcpp_sum(NumericVector v);
RcppExport SEXP _chronoG_rcpp_sum(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_sum(v));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_rbinom
NumericVector rcpp_rbinom(int n, int size, double prob);
RcppExport SEXP _chronoG_rcpp_rbinom(SEXP nSEXP, SEXP sizeSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_rbinom(n, size, prob));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_rpois
NumericVector rcpp_rpois(int n, double lambda);
RcppExport SEXP _chronoG_rcpp_rpois(SEXP nSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_rpois(n, lambda));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_mutate
IntegerVector rcpp_mutate(IntegerVector gt, double mu);
RcppExport SEXP _chronoG_rcpp_mutate(SEXP gtSEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type gt(gtSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_mutate(gt, mu));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_mutate_length_matrix
IntegerMatrix rcpp_mutate_length_matrix(IntegerMatrix gt, double mu, int gens);
RcppExport SEXP _chronoG_rcpp_mutate_length_matrix(SEXP gtSEXP, SEXP muSEXP, SEXP gensSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type gt(gtSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< int >::type gens(gensSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_mutate_length_matrix(gt, mu, gens));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_chronoG_rcpp_sum", (DL_FUNC) &_chronoG_rcpp_sum, 1},
    {"_chronoG_rcpp_rbinom", (DL_FUNC) &_chronoG_rcpp_rbinom, 3},
    {"_chronoG_rcpp_rpois", (DL_FUNC) &_chronoG_rcpp_rpois, 2},
    {"_chronoG_rcpp_mutate", (DL_FUNC) &_chronoG_rcpp_mutate, 2},
    {"_chronoG_rcpp_mutate_length_matrix", (DL_FUNC) &_chronoG_rcpp_mutate_length_matrix, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_chronoG(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}