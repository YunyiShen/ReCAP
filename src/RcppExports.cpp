// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// getLeslie
arma::mat getLeslie(const arma::mat& Surv, const arma::mat& Fec, const double& SRB);
RcppExport SEXP _ReCAP_getLeslie(SEXP SurvSEXP, SEXP FecSEXP, SEXP SRBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Surv(SurvSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Fec(FecSEXP);
    Rcpp::traits::input_parameter< const double& >::type SRB(SRBSEXP);
    rcpp_result_gen = Rcpp::wrap(getLeslie(Surv, Fec, SRB));
    return rcpp_result_gen;
END_RCPP
}
// ProjectAllCpp
List ProjectAllCpp(const arma::mat& Surv, const arma::mat& Harvpar, const arma::mat& Fec, const arma::mat& SRB, const List& aK0, const bool& global, const bool& null, const arma::mat& bl, const int& period, const IntegerVector& nage);
RcppExport SEXP _ReCAP_ProjectAllCpp(SEXP SurvSEXP, SEXP HarvparSEXP, SEXP FecSEXP, SEXP SRBSEXP, SEXP aK0SEXP, SEXP globalSEXP, SEXP nullSEXP, SEXP blSEXP, SEXP periodSEXP, SEXP nageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Surv(SurvSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Harvpar(HarvparSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Fec(FecSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type SRB(SRBSEXP);
    Rcpp::traits::input_parameter< const List& >::type aK0(aK0SEXP);
    Rcpp::traits::input_parameter< const bool& >::type global(globalSEXP);
    Rcpp::traits::input_parameter< const bool& >::type null(nullSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type bl(blSEXP);
    Rcpp::traits::input_parameter< const int& >::type period(periodSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type nage(nageSEXP);
    rcpp_result_gen = Rcpp::wrap(ProjectAllCpp(Surv, Harvpar, Fec, SRB, aK0, global, null, bl, period, nage));
    return rcpp_result_gen;
END_RCPP
}
// getAerialCountPost
arma::mat getAerialCountPost(const List& Proj, const arma::mat& A);
RcppExport SEXP _ReCAP_getAerialCountPost(SEXP ProjSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type Proj(ProjSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(getAerialCountPost(Proj, A));
    return rcpp_result_gen;
END_RCPP
}
// getAerialCountPre
arma::mat getAerialCountPre(const List& Proj, const arma::mat& A);
RcppExport SEXP _ReCAP_getAerialCountPre(SEXP ProjSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type Proj(ProjSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(getAerialCountPre(Proj, A));
    return rcpp_result_gen;
END_RCPP
}
// get_obs_LambdasA
arma::mat get_obs_LambdasA(const arma::mat& Living_total);
RcppExport SEXP _ReCAP_get_obs_LambdasA(SEXP Living_totalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Living_total(Living_totalSEXP);
    rcpp_result_gen = Rcpp::wrap(get_obs_LambdasA(Living_total));
    return rcpp_result_gen;
END_RCPP
}
// get_hypo_Lambdas_helper
arma::mat get_hypo_Lambdas_helper(const arma::mat& Harv_n, const arma::mat& living, const arma::mat& H_n, const arma::mat& Surv_np1, const arma::mat& Fec_np1, const double& SRB_np1, const arma::mat& H_np1);
RcppExport SEXP _ReCAP_get_hypo_Lambdas_helper(SEXP Harv_nSEXP, SEXP livingSEXP, SEXP H_nSEXP, SEXP Surv_np1SEXP, SEXP Fec_np1SEXP, SEXP SRB_np1SEXP, SEXP H_np1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Harv_n(Harv_nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type living(livingSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H_n(H_nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Surv_np1(Surv_np1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Fec_np1(Fec_np1SEXP);
    Rcpp::traits::input_parameter< const double& >::type SRB_np1(SRB_np1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H_np1(H_np1SEXP);
    rcpp_result_gen = Rcpp::wrap(get_hypo_Lambdas_helper(Harv_n, living, H_n, Surv_np1, Fec_np1, SRB_np1, H_np1));
    return rcpp_result_gen;
END_RCPP
}
// get_hypo_Lambdas
arma::mat get_hypo_Lambdas(const arma::mat& Harvest, const arma::mat& Living, const arma::mat& Harvpar, const arma::mat& Surv, const arma::mat& Fec, const arma::mat& SRB);
RcppExport SEXP _ReCAP_get_hypo_Lambdas(SEXP HarvestSEXP, SEXP LivingSEXP, SEXP HarvparSEXP, SEXP SurvSEXP, SEXP FecSEXP, SEXP SRBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Harvest(HarvestSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Living(LivingSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Harvpar(HarvparSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Surv(SurvSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Fec(FecSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type SRB(SRBSEXP);
    rcpp_result_gen = Rcpp::wrap(get_hypo_Lambdas(Harvest, Living, Harvpar, Surv, Fec, SRB));
    return rcpp_result_gen;
END_RCPP
}
// eyes
arma::mat eyes(const int& n);
RcppExport SEXP _ReCAP_eyes(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(eyes(n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ReCAP_getLeslie", (DL_FUNC) &_ReCAP_getLeslie, 3},
    {"_ReCAP_ProjectAllCpp", (DL_FUNC) &_ReCAP_ProjectAllCpp, 10},
    {"_ReCAP_getAerialCountPost", (DL_FUNC) &_ReCAP_getAerialCountPost, 2},
    {"_ReCAP_getAerialCountPre", (DL_FUNC) &_ReCAP_getAerialCountPre, 2},
    {"_ReCAP_get_obs_LambdasA", (DL_FUNC) &_ReCAP_get_obs_LambdasA, 1},
    {"_ReCAP_get_hypo_Lambdas_helper", (DL_FUNC) &_ReCAP_get_hypo_Lambdas_helper, 7},
    {"_ReCAP_get_hypo_Lambdas", (DL_FUNC) &_ReCAP_get_hypo_Lambdas, 6},
    {"_ReCAP_eyes", (DL_FUNC) &_ReCAP_eyes, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_ReCAP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
