# optimization function used by the BayesPrism paper (if optimizer="paper")

#' function to compute log(sum(exp(a_1) + exp(a_2) + exp(a_3) + ...))
logsumexp <- function (x) {
	y = max(x)
	y + log(sum(exp(x - y)))
}


#' function to compute log poserior over gamma_t
log.posterior.gamma <- function (gamma_t,
							     phi_t,
							     phi_t.log,
							     Z_gt.t,
							     Z_t.t,
							     prior.num){
	
	x <- phi_t.log + gamma_t
	psi_t.log <- x - logsumexp(x)
	log.likelihood <- sum(Z_gt.t * psi_t.log, na.rm=TRUE)
	
	log.prior <- sum(prior.num * gamma_t^2) #elementalwise prod
  	log.posterior <- log.likelihood + log.prior
  	
	return(-log.posterior)
}


#' function to compute gradient of log poserior over gamma_t
log.posterior.gamma.grad <- function (gamma_t, 
							   		  phi_t,
							   		  phi_t.log,
							   		  Z_gt.t,
							   		  Z_t.t,
							   		  prior.num){
   	  
  psi_t <- transform.phi_t (phi_t= phi_t, gamma_t)

  log.likelihood.grad <- Z_gt.t - (Z_t.t * psi_t)  

  log.prior.grad <- 2* prior.num * gamma_t
  log.posterior.grad <-  log.likelihood.grad + log.prior.grad
  
  return(-log.posterior.grad)
}




#' function to optimize over gamma
optimize.psi<-function(phi,
					   Z_gt,
					   prior.num,
					   opt.control){
	
	#strip off dimnames to reduce memory use 
	
	phi.dimnames <- dimnames(phi)
	dimnames(phi) <- NULL
	dimnames(Z_gt) <- NULL
					   		
	Z_t <- colSums(Z_gt)
				  			
	cpu.fun <- function(t) {
		require(BayesPrism)
		Rcgminu(par= rep(0,ncol(phi)),
	  			fn= log.posterior.gamma,
	  			gr= log.posterior.gamma.grad,
	  			control= opt.control, 
	  			phi_t = phi[t,],
	  			phi_t.log = log(phi[t,]),
	  			Z_gt.t = Z_gt[,t], 
	  			Z_t.t = Z_t[t],
	  			prior.num = prior.num)
	}
	
	sfInit(parallel = TRUE, cpus = opt.control$n.cores, type = "SOCK" )
	opt.control$n.cores <- NULL
	sfExport("phi", "Z_gt", "Z_t", "prior.num", "opt.control")

	environment(cpu.fun) <- globalenv()
	opt.res <- sfLapply( 1:nrow(phi), cpu.fun)
	sfStop()
	gc()
	
	#check.converge(opt.res)
	
	value <- sum(unlist(lapply(opt.res,"[[", "value")))
	
	opt.gamma <- do.call(rbind,lapply(opt.res, "[[", "par"))
	
	#06-02-2020, bug fix : if too close to the true solution (unstable gradient estimates), sometimes Rcgminu will get stuck at an extreme value
	#set to 0 if extremes > 20 (biologically unlikely value), i.e. use MLE instead	
	opt.gamma[apply(abs(opt.gamma),1,max) > 20,] <- 0
	
	psi <- transform.phi(phi, opt.gamma)
	dimnames(psi) <- phi.dimnames
	
	return(list(psi = psi, value = value))				   	
}










