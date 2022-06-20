#' functions to update references

#' function to check convergence and print error code if there is any
#' @param opt.res a list optimization results
check.converge <- function(opt.res){
	converge.code <- sapply(opt.res, "[[", "convergence")
	all.messages <- sapply(opt.res, "[[", "message")
	
	if(!all(converge.code==0)){
		print(all.messages[converge.code!=0])		
	}
	else{
		cat("Rcgmin seems to have converged \n")
	}
	if(!all(!grepl("Too many gradient evaluations",all.messages))) {
		cat("Please increase maxit in opt.control \n")
	}
	
}


#' function that transforms phi_t with a log fold change gamma_t: softmax(phi_tg*exp(gamma_tg)) 
#' @param phi_t a numeric vector of length G (sums up to one)
#' @param gamma_t a numeric vector of length G
#' @return a numeric vector of length G (sums up to one)
transform.phi_t <- function (phi_t, gamma_t){
		
	stablizing.constant <- max(gamma_t) 		
  	gamma.stab <- gamma_t - stablizing.constant
  	psi_t <- phi_t * exp(gamma.stab)
  	psi_t <- psi_t/sum(psi_t)
  	
  	return(psi_t)
}

#' function that transforms phi matrix with a log fold change matrix gamma
#' @param phi a numeric matrix of dimension T*G (each row sums up to one)
#' @param gamma a numeric matrix of dimension T*G
#' @return a numeric matrix of dimension T*G (each row sums up to one)
transform.phi <- function (phi, gamma){
		
  	psi <- matrix(NA, 
  				  nrow=nrow(phi), ncol=ncol(phi),
  				  dimnames=dimnames(phi))
  	
  	for(t in 1:nrow(psi)) psi[t,] <- transform.phi_t (phi[t,], gamma[t,])
  	
  	return(psi)
}


#' function that generates psi_mal (a N*G matrix) using MLE estimators
#' @param Z_ng_mal a numeric matrix of dimension N*G 
#' @param pseudo.min minimum value used to adjust zero entries
#' @return a numeric matrix of dimension N*G
get.MLE.psi_mal <- function(Z_ng_mal,
							pseudo.min){
					   	  	
	mle.psi_mal <- Z_ng_mal / rowSums(Z_ng_mal)
	
	mle.psi_mal.pseudo.min <- norm.to.one(ref = mle.psi_mal,
										  pseudo.min = pseudo.min)
	
	return(mle.psi_mal.pseudo.min)
}

#' function that updates the reference matrix based on initial Gibbs sampling
#' @param Z a Z_ngk array with cell states merged
#' @param phi_prime a matrix of dimension T*G, with T = number of cell types
#' @param map a list to store the correspondence between cell states and cell types.
#"		of the format list(cell.type1=c(cell.stateA, cell.stateB), ...)
#' @param key a charater string to denote the word that corresponds to the malignant cell type, belongs to names(map)
#'		  set to NULL if there is no malignant cells in the problem.
#' @param optimizer a character string to denote which algorithm to use
#'		  "paper" = the one used by the BayesPrism paper. "EB": the new algorithm that approximates using empirical Bayes-like approach
#' @param sigma a numeric value required if algo="paper"
#' @param opt.control a list containing the parameters to control optimization
#' @return an S$ object of the class "reference". More specifically returns the refPhi class if key=NULL, else return the refTumor class. 
updateReference <- function(Z,
							phi_prime,
							map,
							key,
							optimizer = c("MAP","MLE"),
							opt.control){
	
	cat("Update the reference matrix ... \n")
		
	sigma <- opt.control$sigma
	opt.control$sigma <- NULL
		
	optimizer <- opt.control$optimizer
	opt.control$optimizer <- NULL
	
	if(is.na(key)){
		#if no reference for maligant cells
		Z_gt <- colSums(Z, dims=1) 
		
		if(optimizer == "MAP"){
			psi <- optimize.psi (phi = phi_prime@phi,
					   			   Z_gt = Z_gt,
					   			   prior.num = -1 / (2* sigma ^2),
					   			   opt.control = opt.control)$psi
		}
		if(optimizer == "MLE"){
			psi <- optimize.psi.oneGamma (phi = phi_prime@phi,
					   			  Z_gt = Z_gt,
					   			  opt.control = opt.control)$psi
		}
		return(new("refPhi", phi = psi))
	}
	else{
		# if reference for maligant cells is present
		
		# get MLE for psi tumor in each bulk
		
		Z_ng_mal <- Z[,,key]
		if(is.null(dim(Z_ng_mal))) 
			Z_ng_mal <- matrix(Z_ng_mal, nrow=1, dimnames=dimnames(Z)[1:2])
		
		psi_mal <- get.MLE.psi_mal(Z_ng_mal = Z_ng_mal, 
								   pseudo.min = phi_prime@pseudo.min)
		
		#get MAP for psi 
		cellType.env <-  names(map)[names(map)!= key]
		Z_gt_env <- colSums(Z[,,cellType.env,drop=F], dims=1) 
		phi_env <- phi_prime@phi[cellType.env,,drop=F]
		
		if(optimizer == "MAP"){
			psi_env <- optimize.psi (phi = phi_env,
					   			     Z_gt = Z_gt_env,
					   			     prior.num = -1 / (2* sigma ^2),
					   			     opt.control = opt.control)$psi
		}
		if(optimizer == "MLE"){
			psi_env <- optimize.psi.oneGamma (phi = phi_env,
					   			  Z_gt = Z_gt_env,
					   			  opt.control = opt.control)$psi
		}
		
		return(new("refTumor", 
					psi_mal = psi_mal,
					psi_env = psi_env,
					key = key))
	}
}
