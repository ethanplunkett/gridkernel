// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// calckernelc
double calckernelc(Rcpp::NumericMatrix x, Rcpp::NumericMatrix k, int row, int col);
RcppExport SEXP _gridkernel_calckernelc(SEXP xSEXP, SEXP kSEXP, SEXP rowSEXP, SEXP colSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type row(rowSEXP);
    Rcpp::traits::input_parameter< int >::type col(colSEXP);
    rcpp_result_gen = Rcpp::wrap(calckernelc(x, k, row, col));
    return rcpp_result_gen;
END_RCPP
}
// kernelsmoothc
Rcpp::NumericMatrix kernelsmoothc(Rcpp::NumericMatrix x, Rcpp::NumericMatrix k);
RcppExport SEXP _gridkernel_kernelsmoothc(SEXP xSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(kernelsmoothc(x, k));
    return rcpp_result_gen;
END_RCPP
}
// upscalec
Rcpp::NumericMatrix upscalec(Rcpp::NumericMatrix x, int factor, double maxPNA);
RcppExport SEXP _gridkernel_upscalec(SEXP xSEXP, SEXP factorSEXP, SEXP maxPNASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type factor(factorSEXP);
    Rcpp::traits::input_parameter< double >::type maxPNA(maxPNASEXP);
    rcpp_result_gen = Rcpp::wrap(upscalec(x, factor, maxPNA));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gridkernel_calckernelc", (DL_FUNC) &_gridkernel_calckernelc, 4},
    {"_gridkernel_kernelsmoothc", (DL_FUNC) &_gridkernel_kernelsmoothc, 2},
    {"_gridkernel_upscalec", (DL_FUNC) &_gridkernel_upscalec, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_gridkernel(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
