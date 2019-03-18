// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_rags2ridges_RCPPEXPORTS_H_GEN_
#define RCPP_rags2ridges_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace rags2ridges {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("rags2ridges", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("rags2ridges", "_rags2ridges_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in rags2ridges");
            }
        }
    }

    inline double NLL(const arma::mat S, const arma::mat P) {
        typedef SEXP(*Ptr_NLL)(SEXP,SEXP);
        static Ptr_NLL p_NLL = NULL;
        if (p_NLL == NULL) {
            validateSignature("double(*NLL)(const arma::mat,const arma::mat)");
            p_NLL = (Ptr_NLL)R_GetCCallable("rags2ridges", "_rags2ridges_NLL");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_NLL(Shield<SEXP>(Rcpp::wrap(S)), Shield<SEXP>(Rcpp::wrap(P)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double PNLL(const arma::mat S, const arma::mat P, const arma::mat T, const double lambda) {
        typedef SEXP(*Ptr_PNLL)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_PNLL p_PNLL = NULL;
        if (p_PNLL == NULL) {
            validateSignature("double(*PNLL)(const arma::mat,const arma::mat,const arma::mat,const double)");
            p_PNLL = (Ptr_PNLL)R_GetCCallable("rags2ridges", "_rags2ridges_PNLL");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_PNLL(Shield<SEXP>(Rcpp::wrap(S)), Shield<SEXP>(Rcpp::wrap(P)), Shield<SEXP>(Rcpp::wrap(T)), Shield<SEXP>(Rcpp::wrap(lambda)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double NLL_fused(const Rcpp::List Slist, const Rcpp::List Plist, const arma::vec ns) {
        typedef SEXP(*Ptr_NLL_fused)(SEXP,SEXP,SEXP);
        static Ptr_NLL_fused p_NLL_fused = NULL;
        if (p_NLL_fused == NULL) {
            validateSignature("double(*NLL_fused)(const Rcpp::List,const Rcpp::List,const arma::vec)");
            p_NLL_fused = (Ptr_NLL_fused)R_GetCCallable("rags2ridges", "_rags2ridges_NLL_fused");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_NLL_fused(Shield<SEXP>(Rcpp::wrap(Slist)), Shield<SEXP>(Rcpp::wrap(Plist)), Shield<SEXP>(Rcpp::wrap(ns)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double PNLL_fused(const Rcpp::List Slist, const Rcpp::List Plist, const arma::vec ns, const Rcpp::List Tlist, const arma::mat lambda) {
        typedef SEXP(*Ptr_PNLL_fused)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_PNLL_fused p_PNLL_fused = NULL;
        if (p_PNLL_fused == NULL) {
            validateSignature("double(*PNLL_fused)(const Rcpp::List,const Rcpp::List,const arma::vec,const Rcpp::List,const arma::mat)");
            p_PNLL_fused = (Ptr_PNLL_fused)R_GetCCallable("rags2ridges", "_rags2ridges_PNLL_fused");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_PNLL_fused(Shield<SEXP>(Rcpp::wrap(Slist)), Shield<SEXP>(Rcpp::wrap(Plist)), Shield<SEXP>(Rcpp::wrap(ns)), Shield<SEXP>(Rcpp::wrap(Tlist)), Shield<SEXP>(Rcpp::wrap(lambda)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline arma::mat _armaPooledS(const Rcpp::List& Slist, const Rcpp::NumericVector ns, const int mle = 0) {
        typedef SEXP(*Ptr__armaPooledS)(SEXP,SEXP,SEXP);
        static Ptr__armaPooledS p__armaPooledS = NULL;
        if (p__armaPooledS == NULL) {
            validateSignature("arma::mat(*_armaPooledS)(const Rcpp::List&,const Rcpp::NumericVector,const int)");
            p__armaPooledS = (Ptr__armaPooledS)R_GetCCallable("rags2ridges", "_rags2ridges__armaPooledS");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p__armaPooledS(Shield<SEXP>(Rcpp::wrap(Slist)), Shield<SEXP>(Rcpp::wrap(ns)), Shield<SEXP>(Rcpp::wrap(mle)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::mat _armaPooledP(const Rcpp::List& Plist, const Rcpp::NumericVector ns, const int mle = 0) {
        typedef SEXP(*Ptr__armaPooledP)(SEXP,SEXP,SEXP);
        static Ptr__armaPooledP p__armaPooledP = NULL;
        if (p__armaPooledP == NULL) {
            validateSignature("arma::mat(*_armaPooledP)(const Rcpp::List&,const Rcpp::NumericVector,const int)");
            p__armaPooledP = (Ptr__armaPooledP)R_GetCCallable("rags2ridges", "_rags2ridges__armaPooledP");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p__armaPooledP(Shield<SEXP>(Rcpp::wrap(Plist)), Shield<SEXP>(Rcpp::wrap(ns)), Shield<SEXP>(Rcpp::wrap(mle)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::vec _armaEigShrink(const arma::vec dVec, const double lambda, const double cons = 0) {
        typedef SEXP(*Ptr__armaEigShrink)(SEXP,SEXP,SEXP);
        static Ptr__armaEigShrink p__armaEigShrink = NULL;
        if (p__armaEigShrink == NULL) {
            validateSignature("arma::vec(*_armaEigShrink)(const arma::vec,const double,const double)");
            p__armaEigShrink = (Ptr__armaEigShrink)R_GetCCallable("rags2ridges", "_rags2ridges__armaEigShrink");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p__armaEigShrink(Shield<SEXP>(Rcpp::wrap(dVec)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(cons)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::vec >(rcpp_result_gen);
    }

    inline arma::vec _armaEigShrinkAnyTarget(const arma::mat& S, const arma::mat& target, const double lambda) {
        typedef SEXP(*Ptr__armaEigShrinkAnyTarget)(SEXP,SEXP,SEXP);
        static Ptr__armaEigShrinkAnyTarget p__armaEigShrinkAnyTarget = NULL;
        if (p__armaEigShrinkAnyTarget == NULL) {
            validateSignature("arma::vec(*_armaEigShrinkAnyTarget)(const arma::mat&,const arma::mat&,const double)");
            p__armaEigShrinkAnyTarget = (Ptr__armaEigShrinkAnyTarget)R_GetCCallable("rags2ridges", "_rags2ridges__armaEigShrinkAnyTarget");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p__armaEigShrinkAnyTarget(Shield<SEXP>(Rcpp::wrap(S)), Shield<SEXP>(Rcpp::wrap(target)), Shield<SEXP>(Rcpp::wrap(lambda)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::vec >(rcpp_result_gen);
    }

    inline arma::vec _armaEigShrinkArchI(const arma::vec dVec, const double lambda, const double cons) {
        typedef SEXP(*Ptr__armaEigShrinkArchI)(SEXP,SEXP,SEXP);
        static Ptr__armaEigShrinkArchI p__armaEigShrinkArchI = NULL;
        if (p__armaEigShrinkArchI == NULL) {
            validateSignature("arma::vec(*_armaEigShrinkArchI)(const arma::vec,const double,const double)");
            p__armaEigShrinkArchI = (Ptr__armaEigShrinkArchI)R_GetCCallable("rags2ridges", "_rags2ridges__armaEigShrinkArchI");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p__armaEigShrinkArchI(Shield<SEXP>(Rcpp::wrap(dVec)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(cons)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::vec >(rcpp_result_gen);
    }

    inline arma::mat _armaRidgePAnyTarget(const arma::mat& S, const arma::mat& target, const double lambda, int invert = 2) {
        typedef SEXP(*Ptr__armaRidgePAnyTarget)(SEXP,SEXP,SEXP,SEXP);
        static Ptr__armaRidgePAnyTarget p__armaRidgePAnyTarget = NULL;
        if (p__armaRidgePAnyTarget == NULL) {
            validateSignature("arma::mat(*_armaRidgePAnyTarget)(const arma::mat&,const arma::mat&,const double,int)");
            p__armaRidgePAnyTarget = (Ptr__armaRidgePAnyTarget)R_GetCCallable("rags2ridges", "_rags2ridges__armaRidgePAnyTarget");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p__armaRidgePAnyTarget(Shield<SEXP>(Rcpp::wrap(S)), Shield<SEXP>(Rcpp::wrap(target)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(invert)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::mat _armaRidgePScalarTarget(const arma::mat& S, const double alpha, const double lambda, int invert = 2) {
        typedef SEXP(*Ptr__armaRidgePScalarTarget)(SEXP,SEXP,SEXP,SEXP);
        static Ptr__armaRidgePScalarTarget p__armaRidgePScalarTarget = NULL;
        if (p__armaRidgePScalarTarget == NULL) {
            validateSignature("arma::mat(*_armaRidgePScalarTarget)(const arma::mat&,const double,const double,int)");
            p__armaRidgePScalarTarget = (Ptr__armaRidgePScalarTarget)R_GetCCallable("rags2ridges", "_rags2ridges__armaRidgePScalarTarget");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p__armaRidgePScalarTarget(Shield<SEXP>(Rcpp::wrap(S)), Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(invert)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::mat _armaRidgeP(const arma::mat& S, const arma::mat& target, const double lambda, int invert = 2) {
        typedef SEXP(*Ptr__armaRidgeP)(SEXP,SEXP,SEXP,SEXP);
        static Ptr__armaRidgeP p__armaRidgeP = NULL;
        if (p__armaRidgeP == NULL) {
            validateSignature("arma::mat(*_armaRidgeP)(const arma::mat&,const arma::mat&,const double,int)");
            p__armaRidgeP = (Ptr__armaRidgeP)R_GetCCallable("rags2ridges", "_rags2ridges__armaRidgeP");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p__armaRidgeP(Shield<SEXP>(Rcpp::wrap(S)), Shield<SEXP>(Rcpp::wrap(target)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(invert)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::mat _armaFusedUpdateI(int g0, const Rcpp::List& Plist, const Rcpp::List& Slist, const Rcpp::List& Tlist, const arma::vec& ns, const arma::mat& lambda) {
        typedef SEXP(*Ptr__armaFusedUpdateI)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr__armaFusedUpdateI p__armaFusedUpdateI = NULL;
        if (p__armaFusedUpdateI == NULL) {
            validateSignature("arma::mat(*_armaFusedUpdateI)(int,const Rcpp::List&,const Rcpp::List&,const Rcpp::List&,const arma::vec&,const arma::mat&)");
            p__armaFusedUpdateI = (Ptr__armaFusedUpdateI)R_GetCCallable("rags2ridges", "_rags2ridges__armaFusedUpdateI");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p__armaFusedUpdateI(Shield<SEXP>(Rcpp::wrap(g0)), Shield<SEXP>(Rcpp::wrap(Plist)), Shield<SEXP>(Rcpp::wrap(Slist)), Shield<SEXP>(Rcpp::wrap(Tlist)), Shield<SEXP>(Rcpp::wrap(ns)), Shield<SEXP>(Rcpp::wrap(lambda)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::mat _armaFusedUpdateII(int g0, const Rcpp::List& Plist, const Rcpp::List& Slist, const Rcpp::List& Tlist, const arma::vec ns, const arma::mat lambda) {
        typedef SEXP(*Ptr__armaFusedUpdateII)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr__armaFusedUpdateII p__armaFusedUpdateII = NULL;
        if (p__armaFusedUpdateII == NULL) {
            validateSignature("arma::mat(*_armaFusedUpdateII)(int,const Rcpp::List&,const Rcpp::List&,const Rcpp::List&,const arma::vec,const arma::mat)");
            p__armaFusedUpdateII = (Ptr__armaFusedUpdateII)R_GetCCallable("rags2ridges", "_rags2ridges__armaFusedUpdateII");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p__armaFusedUpdateII(Shield<SEXP>(Rcpp::wrap(g0)), Shield<SEXP>(Rcpp::wrap(Plist)), Shield<SEXP>(Rcpp::wrap(Slist)), Shield<SEXP>(Rcpp::wrap(Tlist)), Shield<SEXP>(Rcpp::wrap(ns)), Shield<SEXP>(Rcpp::wrap(lambda)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::mat _armaFusedUpdateIII(int g0, const Rcpp::List& Plist, const Rcpp::List& Slist, const Rcpp::List& Tlist, const arma::vec& ns, const arma::mat& lambda) {
        typedef SEXP(*Ptr__armaFusedUpdateIII)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr__armaFusedUpdateIII p__armaFusedUpdateIII = NULL;
        if (p__armaFusedUpdateIII == NULL) {
            validateSignature("arma::mat(*_armaFusedUpdateIII)(int,const Rcpp::List&,const Rcpp::List&,const Rcpp::List&,const arma::vec&,const arma::mat&)");
            p__armaFusedUpdateIII = (Ptr__armaFusedUpdateIII)R_GetCCallable("rags2ridges", "_rags2ridges__armaFusedUpdateIII");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p__armaFusedUpdateIII(Shield<SEXP>(Rcpp::wrap(g0)), Shield<SEXP>(Rcpp::wrap(Plist)), Shield<SEXP>(Rcpp::wrap(Slist)), Shield<SEXP>(Rcpp::wrap(Tlist)), Shield<SEXP>(Rcpp::wrap(ns)), Shield<SEXP>(Rcpp::wrap(lambda)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline Rcpp::List _armaRidgeP_fused(const Rcpp::List& Slist, const arma::vec& ns, const Rcpp::List& Tlist, const arma::mat& lambda, const Rcpp::List& Plist, const int maxit = 100, const double eps = 1e-7, const bool relative = true, const bool verbose = false) {
        typedef SEXP(*Ptr__armaRidgeP_fused)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr__armaRidgeP_fused p__armaRidgeP_fused = NULL;
        if (p__armaRidgeP_fused == NULL) {
            validateSignature("Rcpp::List(*_armaRidgeP_fused)(const Rcpp::List&,const arma::vec&,const Rcpp::List&,const arma::mat&,const Rcpp::List&,const int,const double,const bool,const bool)");
            p__armaRidgeP_fused = (Ptr__armaRidgeP_fused)R_GetCCallable("rags2ridges", "_rags2ridges__armaRidgeP_fused");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p__armaRidgeP_fused(Shield<SEXP>(Rcpp::wrap(Slist)), Shield<SEXP>(Rcpp::wrap(ns)), Shield<SEXP>(Rcpp::wrap(Tlist)), Shield<SEXP>(Rcpp::wrap(lambda)), Shield<SEXP>(Rcpp::wrap(Plist)), Shield<SEXP>(Rcpp::wrap(maxit)), Shield<SEXP>(Rcpp::wrap(eps)), Shield<SEXP>(Rcpp::wrap(relative)), Shield<SEXP>(Rcpp::wrap(verbose)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline arma::mat rmvnormal(const int n, arma::rowvec mu, arma::mat sigma) {
        typedef SEXP(*Ptr_rmvnormal)(SEXP,SEXP,SEXP);
        static Ptr_rmvnormal p_rmvnormal = NULL;
        if (p_rmvnormal == NULL) {
            validateSignature("arma::mat(*rmvnormal)(const int,arma::rowvec,arma::mat)");
            p_rmvnormal = (Ptr_rmvnormal)R_GetCCallable("rags2ridges", "_rags2ridges_rmvnormal");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_rmvnormal(Shield<SEXP>(Rcpp::wrap(n)), Shield<SEXP>(Rcpp::wrap(mu)), Shield<SEXP>(Rcpp::wrap(sigma)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::cube _armaRWishart(const int n, const arma::mat& sigma, const double nu) {
        typedef SEXP(*Ptr__armaRWishart)(SEXP,SEXP,SEXP);
        static Ptr__armaRWishart p__armaRWishart = NULL;
        if (p__armaRWishart == NULL) {
            validateSignature("arma::cube(*_armaRWishart)(const int,const arma::mat&,const double)");
            p__armaRWishart = (Ptr__armaRWishart)R_GetCCallable("rags2ridges", "_rags2ridges__armaRWishart");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p__armaRWishart(Shield<SEXP>(Rcpp::wrap(n)), Shield<SEXP>(Rcpp::wrap(sigma)), Shield<SEXP>(Rcpp::wrap(nu)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::cube >(rcpp_result_gen);
    }

    inline arma::cube _armaRInvWishart(const int n, const arma::mat& psi, const double nu) {
        typedef SEXP(*Ptr__armaRInvWishart)(SEXP,SEXP,SEXP);
        static Ptr__armaRInvWishart p__armaRInvWishart = NULL;
        if (p__armaRInvWishart == NULL) {
            validateSignature("arma::cube(*_armaRInvWishart)(const int,const arma::mat&,const double)");
            p__armaRInvWishart = (Ptr__armaRInvWishart)R_GetCCallable("rags2ridges", "_rags2ridges__armaRInvWishart");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p__armaRInvWishart(Shield<SEXP>(Rcpp::wrap(n)), Shield<SEXP>(Rcpp::wrap(psi)), Shield<SEXP>(Rcpp::wrap(nu)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::cube >(rcpp_result_gen);
    }

}

#endif // RCPP_rags2ridges_RCPPEXPORTS_H_GEN_
