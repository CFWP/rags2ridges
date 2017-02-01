////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* ----------------------------------------------------------------------------------------------------------------------------------------------------------
SUPPORT FUNCTIONS FOR RIDGE PRECISION ESTIMATION WITH CHORDAL SUPPORT
INTENDED TO BE INCLUDED IN RAGS2RIDGES-PACKAGE WHEN DIRECT IMPORTATION FROM THAT SOURCE IS POSSIBLE
---------------------------------------------------------------------------------------------------------------------------------------------------------- */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export(".armaRidgePchordalInit")]]
arma::mat armaRidgePchordalInitWorkhorse(arma::mat S, const double lambda, arma::mat target, std::string type, Rcpp::List Cliques, Rcpp::List Separators){
	////////////////////////////////////////////////////////////////////////////
	// construct (guess of) ridge MLE of precision
	// various types (Alt, ArchI, ArchII) of the ridge estimor are implemented
	////////////////////////////////////////////////////////////////////////////

	// start with null matrix
	arma::mat P = arma::zeros<arma::mat>(S.n_rows, S.n_rows);

	// first archetypal ridge estimator
	if (type == "ArchI"){ 
		arma::mat targetInv = arma::inv_sympd(target); 
		// account for clique contributions
		for (int k = 0; k < Cliques.size(); ++k){
			arma::ivec slh = Cliques[k];
			slh = slh - 1;
			arma::uvec cliqNodes = arma::conv_to<arma::uvec>::from(slh);
			P.submat(cliqNodes, cliqNodes) =  P.submat(cliqNodes, cliqNodes) + 
				arma::inv_sympd((1-lambda) * S.submat(cliqNodes, cliqNodes) + lambda * targetInv.submat(cliqNodes, cliqNodes));
		}
		// account for separator contributions
	 	for (int k = 0; k < Separators.size(); ++k){
			arma::ivec slh = Separators[k];
			if (slh.n_elem > 0){
				slh = slh - 1;
				arma::uvec sepNodes = arma::conv_to<arma::uvec>::from(slh);
        	 		P.submat(sepNodes, sepNodes) =  arma::symmatl(P.submat(sepNodes, sepNodes)) - 
					arma::inv_sympd((1-lambda) * S.submat(sepNodes, sepNodes) + lambda * targetInv.submat(sepNodes, sepNodes));
			}
		}
	}

	// second archetypal ridge estimator
	if (type == "ArchII"){ 
		// account for clique contributions
		for (int k = 0; k < Cliques.size(); ++k){
			arma::ivec slh = Cliques[k];
			slh = slh - 1;
			arma::uvec cliqNodes = arma::conv_to<arma::uvec>::from(slh);
			P.submat(cliqNodes, cliqNodes) = arma::symmatl(P.submat(cliqNodes, cliqNodes)) + 
				arma::inv_sympd(S.submat(cliqNodes, cliqNodes) + lambda * arma::eye(cliqNodes.n_elem, cliqNodes.n_elem));
		}
		// account for separator contributions
	 	for (int k = 0; k < Separators.size(); ++k){
			arma::ivec slh = Separators[k];
			if (slh.n_elem > 0){
				slh = slh - 1;
				arma::uvec sepNodes = arma::conv_to<arma::uvec>::from(slh);
				P.submat(sepNodes, sepNodes) = arma::symmatl(P.submat(sepNodes, sepNodes)) - 
					arma::inv_sympd(S.submat(sepNodes, sepNodes) + lambda * arma::eye(sepNodes.n_elem, sepNodes.n_elem));
			}
		}
	}

	// alternative ridge estimator
	if (type == "Alt"){ 
		// account for clique contributions
		for (int k = 0; k < Cliques.size(); ++k){
			arma::ivec slh = Cliques[k];
			slh = slh - 1;
			arma::uvec cliqNodes = arma::conv_to<arma::uvec>::from(slh);
			if (cliqNodes.n_elem > 1){ P.submat(cliqNodes, cliqNodes) = P.submat(cliqNodes, cliqNodes) + 
				armaRidgeP(S.submat(cliqNodes, cliqNodes), target.submat(cliqNodes, cliqNodes), lambda); 
			}
			if (cliqNodes.n_elem == 1){  
				arma::mat slh2 = S.submat(cliqNodes, cliqNodes) - lambda * target.submat(cliqNodes, cliqNodes);
				double slh3 = slh2(0,0);
				P.submat(cliqNodes, cliqNodes) = P.submat(cliqNodes, cliqNodes) + 1/(sqrt(lambda + 0.25f * slh3 * slh3) + slh3 / 2);
			}
		}
		// account for separator contributions
	 	for (int k = 0; k < Separators.size(); ++k){
			arma::ivec slh = Separators[k];
			if (slh.n_elem > 0){
				slh = slh - 1;
				arma::uvec sepNodes = arma::conv_to<arma::uvec>::from(slh);
				if (sepNodes.n_elem > 1){ P.submat(sepNodes, sepNodes) = P.submat(sepNodes, sepNodes) - 
					armaRidgeP(S.submat(sepNodes, sepNodes), target.submat(sepNodes, sepNodes), lambda); 
				}
				if (sepNodes.n_elem == 1){  
					arma::mat slh2 = S.submat(sepNodes, sepNodes) - lambda * target.submat(sepNodes, sepNodes);
					double slh3 = slh2(0,0);
				 	P.submat(sepNodes, sepNodes) = P.submat(sepNodes, sepNodes) - 1/(sqrt(lambda + 0.25f * slh3 * slh3) + slh3 / 2);
				}
			}
		}
	}
	
	// return initial ridge guess for chordal precision matrix
	return(P);
}

// [[Rcpp::export(".armaPenLLreparPforNLM")]]
Rcpp::NumericVector armaPenLLreparPforNLM(const arma::vec x, const arma::mat E1, const arma::mat E2, const arma::mat S, const double lambda, const arma::mat target, const arma::uvec nonzerosR, const arma::uvec nonzerosC){
	/////////////////////////////////////////////////////////////////////////////////////////
	// ridge penalized log-likelihood (and its gradient) for the reparametrized precision matrix.
	// for the reparametrization refer Dahl et al. (2005).
	// passed on to the 'nlm' optimization function
	/////////////////////////////////////////////////////////////////////////////////////////

	// construct precision matrix from alternative parametrization
	const arma::mat P = E1 * arma::diagmat(x) * arma::trans(E2) + E2 * arma::diagmat(x) * arma::trans(E1);

	// return (minus) penalized log-likelihood
	Rcpp::NumericVector penLL(1);
	penLL = -log(det(P)) + sum(arma::diagvec(P * S)) + 0.5 * lambda * sum(arma::diagvec((P-target) * (P-target)));

	// return gradient of (minus) penalized log-likelihood: nonzero elements only
	int p = S.n_rows;
	const arma::mat gradMat = 2 * (S - arma::inv_sympd(P) + lambda * (P - target));
	arma::uvec elemID = (nonzerosC - 1) * p + nonzerosR - 1;
	arma::vec grad = gradMat.elem(elemID);
	penLL.attr("gradient") = grad;

	return penLL;
}

// [[Rcpp::export(".armaPenLLreparP")]]
const double armaPenLLreparP(const arma::vec x, const arma::mat E1, const arma::mat E2, const arma::mat S, const double lambda, const arma::mat target, const arma::uvec nonzerosR, const arma::uvec nonzerosC){
	/////////////////////////////////////////////////////////////////////////////////////////
	// ridge penalized log-likelihood for the reparametrized precision matrix.
	// for the reparametrization refer Dahl et al. (2005)
	// passed on to the 'optim' optimization function
	/////////////////////////////////////////////////////////////////////////////////////////

	// construct precision matrix from alternative parametrization
	const arma::mat P = E1 * arma::diagmat(x) * arma::trans(E2) + E2 * arma::diagmat(x) * arma::trans(E1);

	// return (minus) penalized log-likelihood
    // double logDetP; double detPsign;
    // arma::log_det(logDetP, detPsign, P);	
	// const double penLL = -logDetP + sum(arma::diagvec(P * S)) + 0.5 * lambda * sum(arma::diagvec((P-target) * (P-target)));
	const double penLL = -log(det(P)) + sum(arma::diagvec(P * S)) + 0.5 * lambda * sum(arma::diagvec((P-target) * (P-target)));	
	return penLL;
}

// [[Rcpp::export(".armaPenLLreparPgrad")]]
arma::vec armaPenLLreparPgrad(const arma::vec x, const arma::mat E1, const arma::mat E2, const arma::mat S, const double lambda, const arma::mat target, const arma::uvec nonzerosR, const arma::uvec nonzerosC){
	/////////////////////////////////////////////////////////////////////////////////////////
	// gradient of ridge penalized log-likelihood for the reparametrized precision matrix.
	// for the reparametrization refer Dahl et al. (2005)
	// passed on to the 'optim' optimization function
	/////////////////////////////////////////////////////////////////////////////////////////

	// construct precision matrix from alternative parametrization
	const arma::mat P = E1 * arma::diagmat(x) * arma::trans(E2) + E2 * arma::diagmat(x) * arma::trans(E1);

	// return gradient of (minus) penalized log-likelihood: nonzero elements only
	int p = S.n_rows;
	const arma::mat gradMat = 2 * (arma::symmatl(S) - arma::symmatl(arma::inv_sympd(P)) + lambda * (arma::symmatl(P) - arma::symmatl(target)));
	arma::uvec elemID = (nonzerosC - 1) * p + nonzerosR - 1;
	return gradMat.elem(elemID);
}

// [[Rcpp::export(".armaPenLLreparGradArchI")]]
arma::vec armaPenLLreparGradArchI(const arma::vec x, const arma::mat E1, const arma::mat E2, const arma::mat S, const double lambda, const arma::mat target, const arma::uvec nonzerosR, const arma::uvec nonzerosC){
	/////////////////////////////////////////////////////////////////////////////////////////
	// gradient of 'ArchI'-ridge penalized log-likelihood for the reparametrized precision matrix.
	// for the reparametrization refer Dahl et al. (2005)
	// passed on to the 'optim' optimization function
	/////////////////////////////////////////////////////////////////////////////////////////

	// construct precision matrix from alternative parametrization
	const arma::mat P = E1 * arma::diagmat(x) * arma::trans(E2) + E2 * arma::diagmat(x) * arma::trans(E1);
 
	// return gradient of (minus) penalized log-likelihood: nonzero elements only
	int p = S.n_rows;
	const arma::mat gradMat = 2 * ((1-lambda) * arma::symmatl(S) + lambda * arma::symmatl(arma::inv_sympd(target)) - arma::symmatl(arma::inv_sympd(P)));
	arma::uvec elemID = (nonzerosC - 1) * p + nonzerosR - 1;
	return gradMat.elem(elemID);
}

// [[Rcpp::export(".armaPenLLreparGradArchII")]]
arma::vec armaPenLLreparGradArchII(const arma::vec x, const arma::mat E1, const arma::mat E2, const arma::mat S, const double lambda, const arma::mat target, const arma::uvec nonzerosR, const arma::uvec nonzerosC){
	/////////////////////////////////////////////////////////////////////////////////////////
	// gradient of 'ArchII-ridge penalized log-likelihood for the reparametrized precision matrix.
	// for the reparametrization refer Dahl et al. (2005)
	// passed on to the 'optim' optimization function
	/////////////////////////////////////////////////////////////////////////////////////////

	// construct precision matrix from alternative parametrization
	const arma::mat P = E1 * arma::diagmat(x) * arma::trans(E2) + E2 * arma::diagmat(x) * arma::trans(E1);

	// return gradient of (minus) penalized log-likelihood: nonzero elements only
	int p = S.n_rows;
	const arma::mat gradMat = 2 * (arma::symmatl(S) + lambda * arma::eye(p,p) - arma::inv_sympd(P));
	arma::uvec elemID = (nonzerosC - 1) * p + nonzerosR - 1;
	return gradMat.elem(elemID);
}

