// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cdebug
void cdebug(bool dodebug);
RcppExport SEXP _durmod_cdebug(SEXP dodebugSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< bool >::type dodebug(dodebugSEXP);
    cdebug(dodebug);
    return R_NilValue;
END_RCPP
}
// cloglik
NumericVector cloglik(List dataset, List pset, List control, const bool gdiff, const bool dogradient, const bool dofisher);
RcppExport SEXP _durmod_cloglik(SEXP datasetSEXP, SEXP psetSEXP, SEXP controlSEXP, SEXP gdiffSEXP, SEXP dogradientSEXP, SEXP dofisherSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type dataset(datasetSEXP);
    Rcpp::traits::input_parameter< List >::type pset(psetSEXP);
    Rcpp::traits::input_parameter< List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< const bool >::type gdiff(gdiffSEXP);
    Rcpp::traits::input_parameter< const bool >::type dogradient(dogradientSEXP);
    Rcpp::traits::input_parameter< const bool >::type dofisher(dofisherSEXP);
    rcpp_result_gen = Rcpp::wrap(cloglik(dataset, pset, control, gdiff, dogradient, dofisher));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_durmod_cdebug", (DL_FUNC) &_durmod_cdebug, 1},
    {"_durmod_cloglik", (DL_FUNC) &_durmod_cloglik, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_durmod(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
