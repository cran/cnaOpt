// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>
#include "typedefs.h"

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// C_order_desc
Rcpp::IntegerVector C_order_desc(const Rcpp::NumericVector& x);
RcppExport SEXP _cnaOpt_C_order_desc(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(C_order_desc(x));
    return rcpp_result_gen;
END_RCPP
}
// C_ccoKeep
LogicalVector C_ccoKeep(const NumericVector con, const NumericVector cov);
RcppExport SEXP _cnaOpt_C_ccoKeep(SEXP conSEXP, SEXP covSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type con(conSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type cov(covSEXP);
    rcpp_result_gen = Rcpp::wrap(C_ccoKeep(con, cov));
    return rcpp_result_gen;
END_RCPP
}
// C_msubset
List C_msubset(const dblList x, const LogicalVector s);
RcppExport SEXP _cnaOpt_C_msubset(SEXP xSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const dblList >::type x(xSEXP);
    Rcpp::traits::input_parameter< const LogicalVector >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(C_msubset(x, s));
    return rcpp_result_gen;
END_RCPP
}
// C_getOptim
List C_getOptim(const dblList x);
RcppExport SEXP _cnaOpt_C_getOptim(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const dblList >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(C_getOptim(x));
    return rcpp_result_gen;
END_RCPP
}
// C_append
NumericVector C_append(const NumericVector x, const NumericVector y);
RcppExport SEXP _cnaOpt_C_append(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(C_append(x, y));
    return rcpp_result_gen;
END_RCPP
}
// C_mappend
List C_mappend(const dblList x, const dblList y);
RcppExport SEXP _cnaOpt_C_mappend(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const dblList >::type x(xSEXP);
    Rcpp::traits::input_parameter< const dblList >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(C_mappend(x, y));
    return rcpp_result_gen;
END_RCPP
}
// repTrueFalse
LogicalVector repTrueFalse(int len, int ntrue);
RcppExport SEXP _cnaOpt_repTrueFalse(SEXP lenSEXP, SEXP ntrueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    Rcpp::traits::input_parameter< int >::type ntrue(ntrueSEXP);
    rcpp_result_gen = Rcpp::wrap(repTrueFalse(len, ntrue));
    return rcpp_result_gen;
END_RCPP
}
// C_iterate2
List C_iterate2(dblList dx, dblList dminxy, double Sx_base, double Sy, int blksize, bool verbose);
RcppExport SEXP _cnaOpt_C_iterate2(SEXP dxSEXP, SEXP dminxySEXP, SEXP Sx_baseSEXP, SEXP SySEXP, SEXP blksizeSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< dblList >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< dblList >::type dminxy(dminxySEXP);
    Rcpp::traits::input_parameter< double >::type Sx_base(Sx_baseSEXP);
    Rcpp::traits::input_parameter< double >::type Sy(SySEXP);
    Rcpp::traits::input_parameter< int >::type blksize(blksizeSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(C_iterate2(dx, dminxy, Sx_base, Sy, blksize, verbose));
    return rcpp_result_gen;
END_RCPP
}
// resize
List resize(const List x, const int n);
RcppExport SEXP _cnaOpt_resize(SEXP xSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(resize(x, n));
    return rcpp_result_gen;
END_RCPP
}
// C_mhs_iteration
List C_mhs_iteration(const int m, const LogicalMatrix m_x, const LogicalMatrix m_sol);
RcppExport SEXP _cnaOpt_C_mhs_iteration(SEXP mSEXP, SEXP m_xSEXP, SEXP m_solSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const LogicalMatrix >::type m_x(m_xSEXP);
    Rcpp::traits::input_parameter< const LogicalMatrix >::type m_sol(m_solSEXP);
    rcpp_result_gen = Rcpp::wrap(C_mhs_iteration(m, m_x, m_sol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cnaOpt_C_order_desc", (DL_FUNC) &_cnaOpt_C_order_desc, 1},
    {"_cnaOpt_C_ccoKeep", (DL_FUNC) &_cnaOpt_C_ccoKeep, 2},
    {"_cnaOpt_C_msubset", (DL_FUNC) &_cnaOpt_C_msubset, 2},
    {"_cnaOpt_C_getOptim", (DL_FUNC) &_cnaOpt_C_getOptim, 1},
    {"_cnaOpt_C_append", (DL_FUNC) &_cnaOpt_C_append, 2},
    {"_cnaOpt_C_mappend", (DL_FUNC) &_cnaOpt_C_mappend, 2},
    {"_cnaOpt_repTrueFalse", (DL_FUNC) &_cnaOpt_repTrueFalse, 2},
    {"_cnaOpt_C_iterate2", (DL_FUNC) &_cnaOpt_C_iterate2, 6},
    {"_cnaOpt_resize", (DL_FUNC) &_cnaOpt_resize, 2},
    {"_cnaOpt_C_mhs_iteration", (DL_FUNC) &_cnaOpt_C_mhs_iteration, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_cnaOpt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
