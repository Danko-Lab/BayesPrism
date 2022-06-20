# optimization function using a single gamma for all cell types 


transform.phi_transpose <- function (phi_transpose, gamma){
		
  	psi <- matrix(NA, 
  				  nrow=ncol(phi_transpose), ncol=nrow(phi_transpose),
  				  dimnames=list(dimnames(phi_transpose)[[2]],dimnames(phi_transpose)[[1]]))
  	
  	for(t in 1:nrow(psi)) psi[t,] <- transform.phi_t (phi_transpose[,t], gamma)
  	
  	return(psi)
}



#' function to compute log poserior over gamma_t conditional on mean parameter
log.mle.gamma <- function (gamma, 
						   phi_transpose,
						   phi.log_transpose, #G*K
						   Z_tg,
						   Z_t){
	  
	x <- phi.log_transpose + gamma #K*G
	psi.log <- t(x) - apply(x,2,logsumexp) #K*G
	log.likelihood <- sum(Z_tg * psi.log, na.rm=TRUE)
	  	
	return(-log.likelihood)
}




#' function to compute gradient of log poserior over gamma_t conditional on mean parameter
log.mle.gamma.grad <- function(gamma, 
					   		   phi_transpose,
					   		   phi.log_transpose,
					   		   Z_tg,
					   		   Z_t){
  
  psi <- transform.phi_transpose (phi_transpose = phi_transpose, gamma = gamma)

  log.likelihood.grad <- colSums( Z_tg - (Z_t * psi) )
  
  return(-log.likelihood.grad)
}





#' function to optimize over gamma using an empiricalBayes-like approach
optimize.psi.oneGamma <- function(phi,
					   			  Z_gt,
					   			  opt.control,
					   			  optimizer="Rcgmin"){
						  
	opt.control$n.cores <- NULL
						   		
	#obtain a single MLE estimators for gamma
	
	if(optimizer=="Rcgmin")
		opt.res <- Rcgminu(par= rep(0,ncol(phi)),
	  				   fn= log.mle.gamma,
	  				   gr= log.mle.gamma.grad,
	  				   control= opt.control, 
	  				   phi_transpose = t(phi),
	  				   phi.log_transpose = t(log(phi)),
	  				   Z_tg = t(Z_gt), 
	  				   Z_t = colSums(Z_gt))
	if(optimizer=="BFGS")
		opt.res <- optim(par= rep(0,ncol(phi)),
	  				   fn= log.mle.gamma,
	  				   gr= log.mle.gamma.grad,
	  				   control= opt.control,
	  				   method="BFGS", 
	  				   phi_transpose = t(phi),
	  				   phi.log_transpose = t(log(phi)),
	  				   Z_tg = t(Z_gt), 
	  				   Z_t = colSums(Z_gt))
	
	#check.converge(list(opt.res))

	opt.gamma <- opt.res$par
	value <- opt.res$value
	
	psi <- transform.phi(phi, do.call(rbind, lapply(1:nrow(phi), function(i) opt.gamma)))
	dimnames(psi) <- dimnames(phi)
	
	return(list(psi = psi, value = value, gamma= opt.gamma))				   	
}





