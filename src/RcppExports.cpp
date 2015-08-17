// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// NLL
double NLL(const arma::mat S, const arma::mat P);
RcppExport SEXP rags2ridges_NLL(SEXP SSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type P(PSEXP);
    __result = Rcpp::wrap(NLL(S, P));
    return __result;
END_RCPP
}
// PNLL
double PNLL(const arma::mat S, const arma::mat P, const arma::mat T, const double lambda);
RcppExport SEXP rags2ridges_PNLL(SEXP SSEXP, SEXP PSEXP, SEXP TSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type T(TSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    __result = Rcpp::wrap(PNLL(S, P, T, lambda));
    return __result;
END_RCPP
}
// NLL_fused
double NLL_fused(const Rcpp::List Slist, const Rcpp::List Plist, const arma::vec ns);
RcppExport SEXP rags2ridges_NLL_fused(SEXP SlistSEXP, SEXP PlistSEXP, SEXP nsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Rcpp::List >::type Slist(SlistSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type Plist(PlistSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type ns(nsSEXP);
    __result = Rcpp::wrap(NLL_fused(Slist, Plist, ns));
    return __result;
END_RCPP
}
// PNLL_fused
double PNLL_fused(const Rcpp::List Slist, const Rcpp::List Plist, const arma::vec ns, const Rcpp::List Tlist, const arma::mat lambda);
RcppExport SEXP rags2ridges_PNLL_fused(SEXP SlistSEXP, SEXP PlistSEXP, SEXP nsSEXP, SEXP TlistSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Rcpp::List >::type Slist(SlistSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type Plist(PlistSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type Tlist(TlistSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type lambda(lambdaSEXP);
    __result = Rcpp::wrap(PNLL_fused(Slist, Plist, ns, Tlist, lambda));
    return __result;
END_RCPP
}
// armaPooledS
arma::mat armaPooledS(const Rcpp::List& Slist, const Rcpp::NumericVector ns, const int mle);
RcppExport SEXP rags2ridges_armaPooledS(SEXP SlistSEXP, SEXP nsSEXP, SEXP mleSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Slist(SlistSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< const int >::type mle(mleSEXP);
    __result = Rcpp::wrap(armaPooledS(Slist, ns, mle));
    return __result;
END_RCPP
}
// armaPooledP
arma::mat armaPooledP(const Rcpp::List& Plist, const Rcpp::NumericVector ns, const int mle);
RcppExport SEXP rags2ridges_armaPooledP(SEXP PlistSEXP, SEXP nsSEXP, SEXP mleSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Plist(PlistSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< const int >::type mle(mleSEXP);
    __result = Rcpp::wrap(armaPooledP(Plist, ns, mle));
    return __result;
END_RCPP
}
// armaRidgePAnyTarget
arma::mat armaRidgePAnyTarget(const arma::mat& S, const arma::mat& target, const double lambda, int invert);
RcppExport SEXP rags2ridges_armaRidgePAnyTarget(SEXP SSEXP, SEXP targetSEXP, SEXP lambdaSEXP, SEXP invertSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type invert(invertSEXP);
    __result = Rcpp::wrap(armaRidgePAnyTarget(S, target, lambda, invert));
    return __result;
END_RCPP
}
// armaRidgePScalarTarget
arma::mat armaRidgePScalarTarget(const arma::mat& S, const double alpha, const double lambda, int invert);
RcppExport SEXP rags2ridges_armaRidgePScalarTarget(SEXP SSEXP, SEXP alphaSEXP, SEXP lambdaSEXP, SEXP invertSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type invert(invertSEXP);
    __result = Rcpp::wrap(armaRidgePScalarTarget(S, alpha, lambda, invert));
    return __result;
END_RCPP
}
// armaRidgeP
arma::mat armaRidgeP(const arma::mat& S, const arma::mat& target, const double lambda, int invert);
RcppExport SEXP rags2ridges_armaRidgeP(SEXP SSEXP, SEXP targetSEXP, SEXP lambdaSEXP, SEXP invertSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type invert(invertSEXP);
    __result = Rcpp::wrap(armaRidgeP(S, target, lambda, invert));
    return __result;
END_RCPP
}
// armaFusedUpdateI
arma::mat armaFusedUpdateI(int g0, const Rcpp::List& Plist, const Rcpp::List& Slist, const Rcpp::List& Tlist, const arma::vec& ns, const arma::mat& lambda);
RcppExport SEXP rags2ridges_armaFusedUpdateI(SEXP g0SEXP, SEXP PlistSEXP, SEXP SlistSEXP, SEXP TlistSEXP, SEXP nsSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type g0(g0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Plist(PlistSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Slist(SlistSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Tlist(TlistSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lambda(lambdaSEXP);
    __result = Rcpp::wrap(armaFusedUpdateI(g0, Plist, Slist, Tlist, ns, lambda));
    return __result;
END_RCPP
}
// armaFusedUpdateII
arma::mat armaFusedUpdateII(int g0, const Rcpp::List& Plist, const Rcpp::List& Slist, const Rcpp::List& Tlist, const arma::vec ns, const arma::mat lambda);
RcppExport SEXP rags2ridges_armaFusedUpdateII(SEXP g0SEXP, SEXP PlistSEXP, SEXP SlistSEXP, SEXP TlistSEXP, SEXP nsSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type g0(g0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Plist(PlistSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Slist(SlistSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Tlist(TlistSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type lambda(lambdaSEXP);
    __result = Rcpp::wrap(armaFusedUpdateII(g0, Plist, Slist, Tlist, ns, lambda));
    return __result;
END_RCPP
}
// armaFusedUpdateIII
arma::mat armaFusedUpdateIII(int g0, const Rcpp::List& Plist, const Rcpp::List& Slist, const Rcpp::List& Tlist, const arma::vec& ns, const arma::mat& lambda);
RcppExport SEXP rags2ridges_armaFusedUpdateIII(SEXP g0SEXP, SEXP PlistSEXP, SEXP SlistSEXP, SEXP TlistSEXP, SEXP nsSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type g0(g0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Plist(PlistSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Slist(SlistSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Tlist(TlistSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lambda(lambdaSEXP);
    __result = Rcpp::wrap(armaFusedUpdateIII(g0, Plist, Slist, Tlist, ns, lambda));
    return __result;
END_RCPP
}
// armaRidgeP_fused
Rcpp::List armaRidgeP_fused(const Rcpp::List& Slist, const arma::vec& ns, const Rcpp::List& Tlist, const arma::mat& lambda, const Rcpp::List& Plist, const int maxit, const double eps, const bool verbose);
RcppExport SEXP rags2ridges_armaRidgeP_fused(SEXP SlistSEXP, SEXP nsSEXP, SEXP TlistSEXP, SEXP lambdaSEXP, SEXP PlistSEXP, SEXP maxitSEXP, SEXP epsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Slist(SlistSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Tlist(TlistSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Plist(PlistSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    __result = Rcpp::wrap(armaRidgeP_fused(Slist, ns, Tlist, lambda, Plist, maxit, eps, verbose));
    return __result;
END_RCPP
}
// rmvnormal
arma::mat rmvnormal(const int n, arma::rowvec mu, arma::mat sigma);
RcppExport SEXP rags2ridges_rmvnormal(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    __result = Rcpp::wrap(rmvnormal(n, mu, sigma));
    return __result;
END_RCPP
}
// armaRWishart
arma::cube armaRWishart(const int n, const arma::mat& sigma, const double nu);
RcppExport SEXP rags2ridges_armaRWishart(SEXP nSEXP, SEXP sigmaSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type nu(nuSEXP);
    __result = Rcpp::wrap(armaRWishart(n, sigma, nu));
    return __result;
END_RCPP
}
// armaRInvWishart
arma::cube armaRInvWishart(const int n, const arma::mat& psi, const double nu);
RcppExport SEXP rags2ridges_armaRInvWishart(SEXP nSEXP, SEXP psiSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const double >::type nu(nuSEXP);
    __result = Rcpp::wrap(armaRInvWishart(n, psi, nu));
    return __result;
END_RCPP
}
