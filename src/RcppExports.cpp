// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// NLL
double NLL(const arma::mat S, const arma::mat P);
RcppExport SEXP rags2ridges_NLL(SEXP SSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(NLL(S, P));
    return rcpp_result_gen;
END_RCPP
}
// PNLL
double PNLL(const arma::mat S, const arma::mat P, const arma::mat T, const double lambda);
RcppExport SEXP rags2ridges_PNLL(SEXP SSEXP, SEXP PSEXP, SEXP TSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type T(TSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(PNLL(S, P, T, lambda));
    return rcpp_result_gen;
END_RCPP
}
// NLL_fused
double NLL_fused(const Rcpp::List Slist, const Rcpp::List Plist, const arma::vec ns);
RcppExport SEXP rags2ridges_NLL_fused(SEXP SlistSEXP, SEXP PlistSEXP, SEXP nsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type Slist(SlistSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type Plist(PlistSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type ns(nsSEXP);
    rcpp_result_gen = Rcpp::wrap(NLL_fused(Slist, Plist, ns));
    return rcpp_result_gen;
END_RCPP
}
// PNLL_fused
double PNLL_fused(const Rcpp::List Slist, const Rcpp::List Plist, const arma::vec ns, const Rcpp::List Tlist, const arma::mat lambda);
RcppExport SEXP rags2ridges_PNLL_fused(SEXP SlistSEXP, SEXP PlistSEXP, SEXP nsSEXP, SEXP TlistSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type Slist(SlistSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type Plist(PlistSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type Tlist(TlistSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(PNLL_fused(Slist, Plist, ns, Tlist, lambda));
    return rcpp_result_gen;
END_RCPP
}
// armaPooledS
arma::mat armaPooledS(const Rcpp::List& Slist, const Rcpp::NumericVector ns, const int mle);
RcppExport SEXP rags2ridges_armaPooledS(SEXP SlistSEXP, SEXP nsSEXP, SEXP mleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Slist(SlistSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< const int >::type mle(mleSEXP);
    rcpp_result_gen = Rcpp::wrap(armaPooledS(Slist, ns, mle));
    return rcpp_result_gen;
END_RCPP
}
// armaPooledP
arma::mat armaPooledP(const Rcpp::List& Plist, const Rcpp::NumericVector ns, const int mle);
RcppExport SEXP rags2ridges_armaPooledP(SEXP PlistSEXP, SEXP nsSEXP, SEXP mleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Plist(PlistSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< const int >::type mle(mleSEXP);
    rcpp_result_gen = Rcpp::wrap(armaPooledP(Plist, ns, mle));
    return rcpp_result_gen;
END_RCPP
}
// armaEigShrink
arma::vec armaEigShrink(const arma::vec dVec, const double lambda, const double cons);
RcppExport SEXP rags2ridges_armaEigShrink(SEXP dVecSEXP, SEXP lambdaSEXP, SEXP consSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type dVec(dVecSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type cons(consSEXP);
    rcpp_result_gen = Rcpp::wrap(armaEigShrink(dVec, lambda, cons));
    return rcpp_result_gen;
END_RCPP
}
// armaEigShrinkAnyTarget
arma::vec armaEigShrinkAnyTarget(const arma::mat& S, const arma::mat& target, const double lambda);
RcppExport SEXP rags2ridges_armaEigShrinkAnyTarget(SEXP SSEXP, SEXP targetSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(armaEigShrinkAnyTarget(S, target, lambda));
    return rcpp_result_gen;
END_RCPP
}
// armaEigShrinkArchI
arma::vec armaEigShrinkArchI(const arma::vec dVec, const double lambda, const double cons);
RcppExport SEXP rags2ridges_armaEigShrinkArchI(SEXP dVecSEXP, SEXP lambdaSEXP, SEXP consSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type dVec(dVecSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type cons(consSEXP);
    rcpp_result_gen = Rcpp::wrap(armaEigShrinkArchI(dVec, lambda, cons));
    return rcpp_result_gen;
END_RCPP
}
// armaRidgePAnyTarget
arma::mat armaRidgePAnyTarget(const arma::mat& S, const arma::mat& target, const double lambda, int invert);
RcppExport SEXP rags2ridges_armaRidgePAnyTarget(SEXP SSEXP, SEXP targetSEXP, SEXP lambdaSEXP, SEXP invertSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type invert(invertSEXP);
    rcpp_result_gen = Rcpp::wrap(armaRidgePAnyTarget(S, target, lambda, invert));
    return rcpp_result_gen;
END_RCPP
}
// armaRidgePScalarTarget
arma::mat armaRidgePScalarTarget(const arma::mat& S, const double alpha, const double lambda, int invert);
RcppExport SEXP rags2ridges_armaRidgePScalarTarget(SEXP SSEXP, SEXP alphaSEXP, SEXP lambdaSEXP, SEXP invertSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type invert(invertSEXP);
    rcpp_result_gen = Rcpp::wrap(armaRidgePScalarTarget(S, alpha, lambda, invert));
    return rcpp_result_gen;
END_RCPP
}
// armaRidgeP
arma::mat armaRidgeP(const arma::mat& S, const arma::mat& target, const double lambda, int invert);
RcppExport SEXP rags2ridges_armaRidgeP(SEXP SSEXP, SEXP targetSEXP, SEXP lambdaSEXP, SEXP invertSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type invert(invertSEXP);
    rcpp_result_gen = Rcpp::wrap(armaRidgeP(S, target, lambda, invert));
    return rcpp_result_gen;
END_RCPP
}
// armaFusedUpdateI
arma::mat armaFusedUpdateI(int g0, const Rcpp::List& Plist, const Rcpp::List& Slist, const Rcpp::List& Tlist, const arma::vec& ns, const arma::mat& lambda);
RcppExport SEXP rags2ridges_armaFusedUpdateI(SEXP g0SEXP, SEXP PlistSEXP, SEXP SlistSEXP, SEXP TlistSEXP, SEXP nsSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type g0(g0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Plist(PlistSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Slist(SlistSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Tlist(TlistSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(armaFusedUpdateI(g0, Plist, Slist, Tlist, ns, lambda));
    return rcpp_result_gen;
END_RCPP
}
// armaFusedUpdateII
arma::mat armaFusedUpdateII(int g0, const Rcpp::List& Plist, const Rcpp::List& Slist, const Rcpp::List& Tlist, const arma::vec ns, const arma::mat lambda);
RcppExport SEXP rags2ridges_armaFusedUpdateII(SEXP g0SEXP, SEXP PlistSEXP, SEXP SlistSEXP, SEXP TlistSEXP, SEXP nsSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type g0(g0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Plist(PlistSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Slist(SlistSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Tlist(TlistSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(armaFusedUpdateII(g0, Plist, Slist, Tlist, ns, lambda));
    return rcpp_result_gen;
END_RCPP
}
// armaFusedUpdateIII
arma::mat armaFusedUpdateIII(int g0, const Rcpp::List& Plist, const Rcpp::List& Slist, const Rcpp::List& Tlist, const arma::vec& ns, const arma::mat& lambda);
RcppExport SEXP rags2ridges_armaFusedUpdateIII(SEXP g0SEXP, SEXP PlistSEXP, SEXP SlistSEXP, SEXP TlistSEXP, SEXP nsSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type g0(g0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Plist(PlistSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Slist(SlistSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Tlist(TlistSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(armaFusedUpdateIII(g0, Plist, Slist, Tlist, ns, lambda));
    return rcpp_result_gen;
END_RCPP
}
// armaRidgeP_fused
Rcpp::List armaRidgeP_fused(const Rcpp::List& Slist, const arma::vec& ns, const Rcpp::List& Tlist, const arma::mat& lambda, const Rcpp::List& Plist, const int maxit, const double eps, const bool relative, const bool verbose);
RcppExport SEXP rags2ridges_armaRidgeP_fused(SEXP SlistSEXP, SEXP nsSEXP, SEXP TlistSEXP, SEXP lambdaSEXP, SEXP PlistSEXP, SEXP maxitSEXP, SEXP epsSEXP, SEXP relativeSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Slist(SlistSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Tlist(TlistSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Plist(PlistSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const bool >::type relative(relativeSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(armaRidgeP_fused(Slist, ns, Tlist, lambda, Plist, maxit, eps, relative, verbose));
    return rcpp_result_gen;
END_RCPP
}
// rmvnormal
arma::mat rmvnormal(const int n, arma::rowvec mu, arma::mat sigma);
RcppExport SEXP rags2ridges_rmvnormal(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(rmvnormal(n, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// armaRWishart
arma::cube armaRWishart(const int n, const arma::mat& sigma, const double nu);
RcppExport SEXP rags2ridges_armaRWishart(SEXP nSEXP, SEXP sigmaSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(armaRWishart(n, sigma, nu));
    return rcpp_result_gen;
END_RCPP
}
// armaRInvWishart
arma::cube armaRInvWishart(const int n, const arma::mat& psi, const double nu);
RcppExport SEXP rags2ridges_armaRInvWishart(SEXP nSEXP, SEXP psiSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(armaRInvWishart(n, psi, nu));
    return rcpp_result_gen;
END_RCPP
}
// armaRidgePchordalInitWorkhorse
arma::mat armaRidgePchordalInitWorkhorse(arma::mat S, const double lambda, arma::mat target, std::string type, Rcpp::List Cliques, Rcpp::List Separators);
RcppExport SEXP rags2ridges_armaRidgePchordalInitWorkhorse(SEXP SSEXP, SEXP lambdaSEXP, SEXP targetSEXP, SEXP typeSEXP, SEXP CliquesSEXP, SEXP SeparatorsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type target(targetSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Cliques(CliquesSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Separators(SeparatorsSEXP);
    rcpp_result_gen = Rcpp::wrap(armaRidgePchordalInitWorkhorse(S, lambda, target, type, Cliques, Separators));
    return rcpp_result_gen;
END_RCPP
}
// armaPenLLreparPforNLM
Rcpp::NumericVector armaPenLLreparPforNLM(const arma::vec x, const arma::mat E1, const arma::mat E2, const arma::mat S, const double lambda, const arma::mat target, const arma::uvec nonzerosR, const arma::uvec nonzerosC);
RcppExport SEXP rags2ridges_armaPenLLreparPforNLM(SEXP xSEXP, SEXP E1SEXP, SEXP E2SEXP, SEXP SSEXP, SEXP lambdaSEXP, SEXP targetSEXP, SEXP nonzerosRSEXP, SEXP nonzerosCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type E1(E1SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type E2(E2SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type nonzerosR(nonzerosRSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type nonzerosC(nonzerosCSEXP);
    rcpp_result_gen = Rcpp::wrap(armaPenLLreparPforNLM(x, E1, E2, S, lambda, target, nonzerosR, nonzerosC));
    return rcpp_result_gen;
END_RCPP
}
// armaPenLLreparP
const double armaPenLLreparP(const arma::vec x, const arma::mat E1, const arma::mat E2, const arma::mat S, const double lambda, const arma::mat target, const arma::uvec nonzerosR, const arma::uvec nonzerosC);
RcppExport SEXP rags2ridges_armaPenLLreparP(SEXP xSEXP, SEXP E1SEXP, SEXP E2SEXP, SEXP SSEXP, SEXP lambdaSEXP, SEXP targetSEXP, SEXP nonzerosRSEXP, SEXP nonzerosCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type E1(E1SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type E2(E2SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type nonzerosR(nonzerosRSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type nonzerosC(nonzerosCSEXP);
    rcpp_result_gen = Rcpp::wrap(armaPenLLreparP(x, E1, E2, S, lambda, target, nonzerosR, nonzerosC));
    return rcpp_result_gen;
END_RCPP
}
// armaPenLLreparPgrad
arma::vec armaPenLLreparPgrad(const arma::vec x, const arma::mat E1, const arma::mat E2, const arma::mat S, const double lambda, const arma::mat target, const arma::uvec nonzerosR, const arma::uvec nonzerosC);
RcppExport SEXP rags2ridges_armaPenLLreparPgrad(SEXP xSEXP, SEXP E1SEXP, SEXP E2SEXP, SEXP SSEXP, SEXP lambdaSEXP, SEXP targetSEXP, SEXP nonzerosRSEXP, SEXP nonzerosCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type E1(E1SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type E2(E2SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type nonzerosR(nonzerosRSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type nonzerosC(nonzerosCSEXP);
    rcpp_result_gen = Rcpp::wrap(armaPenLLreparPgrad(x, E1, E2, S, lambda, target, nonzerosR, nonzerosC));
    return rcpp_result_gen;
END_RCPP
}
// armaPenLLreparGradArchI
arma::vec armaPenLLreparGradArchI(const arma::vec x, const arma::mat E1, const arma::mat E2, const arma::mat S, const double lambda, const arma::mat target, const arma::uvec nonzerosR, const arma::uvec nonzerosC);
RcppExport SEXP rags2ridges_armaPenLLreparGradArchI(SEXP xSEXP, SEXP E1SEXP, SEXP E2SEXP, SEXP SSEXP, SEXP lambdaSEXP, SEXP targetSEXP, SEXP nonzerosRSEXP, SEXP nonzerosCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type E1(E1SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type E2(E2SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type nonzerosR(nonzerosRSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type nonzerosC(nonzerosCSEXP);
    rcpp_result_gen = Rcpp::wrap(armaPenLLreparGradArchI(x, E1, E2, S, lambda, target, nonzerosR, nonzerosC));
    return rcpp_result_gen;
END_RCPP
}
// armaPenLLreparGradArchII
arma::vec armaPenLLreparGradArchII(const arma::vec x, const arma::mat E1, const arma::mat E2, const arma::mat S, const double lambda, const arma::mat target, const arma::uvec nonzerosR, const arma::uvec nonzerosC);
RcppExport SEXP rags2ridges_armaPenLLreparGradArchII(SEXP xSEXP, SEXP E1SEXP, SEXP E2SEXP, SEXP SSEXP, SEXP lambdaSEXP, SEXP targetSEXP, SEXP nonzerosRSEXP, SEXP nonzerosCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type E1(E1SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type E2(E2SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type nonzerosR(nonzerosRSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type nonzerosC(nonzerosCSEXP);
    rcpp_result_gen = Rcpp::wrap(armaPenLLreparGradArchII(x, E1, E2, S, lambda, target, nonzerosR, nonzerosC));
    return rcpp_result_gen;
END_RCPP
}
