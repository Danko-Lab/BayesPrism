

#' E step function
#####################################################################################################

#' function to sample Z and theta conditional on each individual bulk
#' @param X_n a numeric vector of the count vector of Nth bulk sample
#' @param theta.ini_n the input theta vector by concatenating omega and theta_env of the Nth bulk sample
#' @param phi.hat a reference matrix by combining the rows of eta and psi_env
#' @param alpha the dirichlet hyper-parameter
#' @param gibbs.idx a numeric vector to denote the index of samples to be retained from MCMC chain 
#' @prarm conditional.idx rows of psi_env in phi.hat to be fixed during sampling
#' @param compute.elbo logical to denote if compute ELBO
sample.n <- function(X_n, theta.ini_n, phi.hat, alpha, gibbs.idx, conditional.idx, compute.elbo){
		
	G <- ncol(phi.hat)
	K.tot <- nrow(phi.hat)
	
	theta_n.i <- theta.ini_n
	
	theta_n.conditional <- theta.ini_n[conditional.idx]
	tau_n <- 1-sum(theta_n.conditional)
	free.idx <- 1: (length(theta.ini_n) - length(theta_n.conditional))
	
	Z_n.i <- array(NA, c(G, K.tot))
	theta_n.sum <- rep(0, length(theta.ini_n))
	Z_n.sum <- array(0,c(G, K.tot))
	
	multinom.coef <- 0
	z.logtheta <- 0
	theta.dirichlet.alpha <- 0
	
	for(i in 1:max(gibbs.idx)){
	
		#sample Z for patient i
		prob.mat <- phi.hat * theta_n.i #multiply by each column
		
		for (g in 1:G) Z_n.i[g,] <- rmultinom(n = 1, size = X_n[g], prob = prob.mat[,g])
						
		# sample theta for patient i
		Z_nk.i <- colSums(Z_n.i) #total count for each cell type
		theta_n.i <- rdirichlet(alpha = Z_nk.i[free.idx] + alpha)
		theta_n.i <- c(theta_n.i * tau_n, theta_n.conditional)
				
		if(i %in% gibbs.idx) {
			Z_n.sum <- Z_n.sum + Z_n.i
			theta_n.sum <- theta_n.sum + theta_n.i
			
			if(compute.elbo){
				multinom.coef <-  multinom.coef  - sum(lfactorial(Z_n.i))
				z.logtheta <- z.logtheta + sum(Z_nk.i[theta_n.i!=0] * log(theta_n.i[theta_n.i!=0]))
				theta.dirichlet.alpha <- theta.dirichlet.alpha + sum((alpha-1) * log(theta_n.i[theta_n.i!=0]))
			}
		}
		if((i %% 50) == 0) gc()
	}
	
	samples.size <- length(gibbs.idx)
	
	gibbs.constant <- multinom.coef + z.logtheta + theta.dirichlet.alpha
	
	return(list(Z_n = Z_n.sum / samples.size, 
				theta_n = theta_n.sum / samples.size, 
				gibbs.constant = gibbs.constant / samples.size))
}
	
	
#' function to run Gibbs sampling for embeding learning ()
#'
#' @param eta, a tumor embedding matrix of dimension K*G (K=# of tumor embeddings/factors)
#' @param psi_env, a matrix of non-malignant cell expression of profile of dimension (T-1) * G 
#' @param theta_env theta matrix of dimension N*(T-1) to denote  non-malignant cell fractions 
#' @param X bulk matrix of N*G
#' @param compute.elbo logical to denote if compute ELBO
#'
#' @import snowfall
#' @return a list of a jointPost object and a number to represent terms in posterior probability 
run.gibbs.refEta <- function(eta,
							 psi_env,
							 theta_env,
				  			 X,
				  			 gibbs.control, 
				  			 compute.elbo){
	
	K.tum <- nrow(eta)
	K.tot <- nrow(eta) + nrow(psi_env)
	
	conditional.idx <- (K.tum+1): K.tot	
	phi.hat <- rbind(eta, psi_env)

	
	#devide the tumor fraction of each sample by K 
	tau <- 1 - rowSums(theta_env)
	omega.ini.k <- tau / K.tum
	omega.ini <- do.call(rbind, lapply(omega.ini.k, function(omega.ini.k_n) rep(omega.ini.k_n, K.tum)))
	colnames(omega.ini) <- rownames(eta)
	
	#concatenate omega.ini and theta_env
	theta.ini <- cbind(omega.ini, theta_env)

	#get control parameters
	gibbs.control <- gibbs.control
	alpha <- gibbs.control$alpha
	
	gibbs.idx <- get.gibbs.idx(gibbs.control)
	seed <- gibbs.control$seed
	
	sfInit(parallel = TRUE, cpus = gibbs.control$n.cores, type = "SOCK" )
	sfExport("phi.hat", "X", "alpha", "conditional.idx", "theta.ini", "gibbs.idx", "seed", "compute.elbo")

	cpu.fun <- function(n) {
		if(!is.null(seed)) set.seed(seed)
		require(BayesPrism)
		sample.n(X_n = X[n,], 
				 theta.ini_n = theta.ini[n,], 
				 phi.hat = phi.hat, 
				 alpha = alpha, 
				 gibbs.idx = gibbs.idx, 
				 conditional.idx = conditional.idx, 
				 compute.elbo = compute.elbo)
	}
	environment(cpu.fun) <- globalenv()
	gibbs.list <- sfLapply( 1:nrow(X), cpu.fun)
	sfStop()
		
	jointPost <- newJointPost( bulkID = rownames(X),
						   	   geneID = colnames(X), 
						   	   cellType = c(rownames(eta), rownames(psi_env)),
						   	   gibbs.list = gibbs.list )
	
	return(jointPost)
}
	
	
#' function to compute ELBO
#' @param opt.value the log.posterior term from opt.psi
#' @param psi_env updated reference for non-malignant cells
#' @param jointPost a jointPost object from run.gibbs.refEta
#'
#' @return a numeric value of ELBO
compute.elbo <- function(opt.value,
						 psi_env,
						 jointPost){
		
	Z_gk_env <- colSums(jointPost@Z[,,rownames(psi_env)], dims=1)
	
	#add back the minus likelihood term of environment cells
	elbo.env <- - sum(log(psi_env) * t(Z_gk_env), na.rm=TRUE) #elemental wise prod

	#add back the constant from dirichlet multionomial. note the sign is inversed
	elbo <-  opt.value + elbo.env - jointPost@constant
	
	return(elbo)					  	
}





#' function to run Gibbs_EM 
#'
#' @param eta_prior, a matrix of dimension K*G (K=# of tumor embeddings/factors) to denote the prior means of eta
#' @param eta_post, a matrix of dimension K*G (K=# of tumor embeddings/factors) to denote the initial point of eta
#' @param psi_env, a matrix of non-malignant cell expression of profile of dimension (T-1) * G 
#' @param theta_env theta matrix of dimension N*(T-1) to denote  non-malignant cell fractions 
#' @param X bulk matrix of N*G
#' @param cycle number of EM cycles
#' @param gibbs.control a list containing parameters of the Gibbs sampler
#'		chain.length: length of MCMC chain. Default=1000;
#'		burn.in: length of burn in period. Default=500;
#'		thinning: retain every # of MCMC samples after the burn in period to reduce auto-correlation. Default=2;
#'		n.cores: number of cores to use. Default =1. Will be over-written by n.cores in the main argument.
#'		seed: seed number to use for repoducibility. Default = 123. set to NULL if use pseudo random.
#'		alpha: a numeric vector to represent the parameter of dirichlet distribution 
#' @param opt.control a list containing parameters for the optimization step:
#'		maxit: maximum number of cycles to interate. Default=100000
#'		sigma: hyper-parameter of the prior if optimizer="paper". Default=2.
#' @param compute.elbo logical to denote if compute ELBO
#'
#' @return a list of eta, omega and ELBO to represent the EM results

run.EM <- function(eta_prior,
				   eta_post,
				   psi_env,
				   theta_env,
				   X,
				   cycle,
				   gibbs.control,
				   opt.control,
				   compute.elbo){


	K.tum <- nrow(eta_prior)
	sigma <- opt.control$sigma
	opt.control$sigma <- NULL
	
	
	every.n <- 1 #default report ELBO for each EM cycle
	
	if(is.numeric(compute.elbo)) {
		#compute ELBO
		every.n <- compute.elbo #report ELBO every N cycles
		compute.elbo <- TRUE
	}
	
	
	#initiate starting point
	if(is.null(eta_post)) eta_post <- eta_prior
		
	elbo.vec <- Inf
	
	ptm <- proc.time()

	EM.cycle.idx<-1	
	
	while(EM.cycle.idx <= cycle){ #EM convergence criterion
		cat("Starting EM cycles # ", EM.cycle.idx, "\n")
		#print("running gibbs sampling...")
		
		if((EM.cycle.idx %% 5) == 0) gc()
		
		# E step using Gibbs sampling
		gibbs.res <- run.gibbs.refEta (eta = eta_post,
							 		   psi_env = psi_env,
							 		   theta_env = theta_env,
				  			 		   X = X, 
				  			 		   gibbs.control = gibbs.control,
				  			 		   compute.elbo = compute.elbo)
		
		#M step
		opt.res <- optimize.psi(phi = eta_prior,
					   			Z_gt = colSums(gibbs.res@Z[,,1:K.tum], dims=1),
					   			prior.num = -1 / (2* sigma ^2),
					   			opt.control = opt.control)
		
		#update eta	   
		eta_post <- opt.res$psi		
		
		if(compute.elbo & (EM.cycle.idx %% every.n) == 0){
			elbo <- compute.elbo (opt.value = opt.res$value,
						 	      psi_env = psi_env,
						 	      jointPost = gibbs.res)
		 elbo.vec <- c(elbo.vec, elbo)
		 cat("ELBO = ", elbo, "\n")
		}
						
		#estimate total EM time
		if(EM.cycle.idx==1){
		  total.time <- proc.time() - ptm
		  total.time <- as.numeric(total.time["elapsed"])
		  estimated.time <- total.time * cycle
		
		  current.time <- Sys.time()
	
		  cat("Current time: ", as.character(current.time), "\n")
		  cat("Estimated time to complete: ", my_seconds_to_period(estimated.time), "\n")
		  cat("Estimated finishing time: ", as.character(current.time + estimated.time), "\n")
		}
		
		EM.cycle.idx <- EM.cycle.idx+1	

	}

	cat("\n")
	omega <- gibbs.res@theta[,1:nrow(eta_post)]

	return(list(eta_prior = eta_prior,
				eta_post = eta_post,
				omega = omega, 
				elbo = elbo.vec))					
}


#' @param bp a BayesPrism object
#' @param eta_prior, a K*G matrix (each row is normalized to sum to one)
#' @param cycle number of EM cycles
#' @param gibbs.control a list containing parameters of the Gibbs sampler
#'		chain.length: length of MCMC chain. Default=1000;
#'		burn.in: length of burn in period. Default=500;
#'		thinning: retain every # of MCMC samples after the burn in period to reduce auto-correlation. Default=2;
#'		n.cores: number of cores to use. Default =1. Will be over-written by n.cores in the main argument.
#'		seed: seed number to use for repoducibility. Default = 123. set to NULL if use pseudo random.
#'		alpha: a numeric vector to represent the parameter of dirichlet distribution 
#' @param opt.control a list containing parameters for the optimization step:
#'		maxit: maximum number of cycles to interate. Default=100000
#'		sigma: hyper-parameter of the prior if optimizer="paper". Default=2.
#' @param EM.res a previous EM result containing eta and omega. eta from EM.res will be used to initate additional EM runs if supplied.
#' @param compute.elbo logical to denote if compute ELBO
#' @return a list of eta, omega and ELBO to represent the EM results

learn.embedding <- function(bp,
							eta_prior = NULL,
							cycle = 50,
							gibbs.control=list(),
							opt.control=list(),
							EM.res = NULL,
							compute.elbo = F){
	
	#determine if bp.res has updated gibbs result
	if(!bp@control_param$update.gibbs)
		stop("please update gibbs sampling first")
	#dtermine if run using a tumor mode
	if(!is(bp@reference.update, "refTumor"))
		stop("please run BayesPrism using tumor mode by specifying the key parameter")

	#use parameters from BayesPrism run if no change
	gibbs.control.bp <- bp@control_param$gibbs.control
	opt.control.bp <- bp@control_param$opt.control
	
	gibbs.control.bp[names(gibbs.control)] <- gibbs.control
	opt.control.bp[names(opt.control)] <- opt.control

	opt.control <- valid.opt.control(opt.control.bp)
	gibbs.control <- valid.gibbs.control(gibbs.control.bp)

	opt.control$Z_t.min <- NULL		  			
	opt.control$optimizer <- NULL
	

	
	if(!is.null(eta_prior)){
		#normalize input eta
		eta_prior <- norm.to.one (ref = eta_prior, pseudo.min = bp@prism@phi_cellState@pseudo.min)
		if(is.null(rownames(eta_prior))) rownames(eta_prior) <- paste("program", 1:nrow(eta_prior), sep="-")
	}
	
	#extract results from previous run
	if(!is.null(EM.res)){ 
		cat("Using eta_prior from previous EM cycles. eta_prior will be overwritten. \n")
		eta_prior <- EM.res$eta_prior
	}
	
	if(is.null(eta_prior))
		stop("Please provide either eta_prior or EM.res \n")

	if(ncol(eta_prior)!= ncol(bp@reference.update@psi_env))
		stop("Number of genes in eta_prior needs to match that in bp.")
	
	EM.res.new <- run.EM (eta_prior = eta_prior,
						  eta_post = EM.res$eta_post,
				   		  psi_env = bp@reference.update@psi_env,
				   		  theta_env = bp@posterior.theta_f@theta[,rownames(bp@reference.update@psi_env)], 
				   		  X = bp@prism@mixture,
				   		  cycle = cycle,
				   		  gibbs.control = gibbs.control,
				   		  opt.control = opt.control,
				   		  compute.elbo = compute.elbo)
	
	if(!is.null(EM.res) & compute.elbo) 
		EM.res.new$elbo[1] <- EM.res$elbo[length(EM.res$elbo)]
				
	return(EM.res.new)

}






#' function to validate nmf.control
#' @param control a named list of parameters required to control NMF  
valid.nmf.control <- function(control){
	ctrl <- list(nrun=200, seed=123, .opt=list())
	namc <- names(control)
	
	if (!all(namc %in% names(ctrl)))
        stop("unknown names in nmf.control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control

	if(ctrl$nrun <0) stop("nrun needs to be positive")
	if(ctrl$seed <0) stop("seed needs to be positive")
	
	
	return(ctrl)
}




#' @param bp a BayesPrism object
#' @param K number of tumor programs (embeddings)
#' @param cycle number of EM cycles
#' @param gibbs.control a list containing parameters of the Gibbs sampler
#'		chain.length: length of MCMC chain. Default=1000;
#'		burn.in: length of burn in period. Default=500;
#'		thinning: retain every # of MCMC samples after the burn in period to reduce auto-correlation. Default=2;
#'		n.cores: number of cores to use. Default =1. Will be over-written by n.cores in the main argument.
#'		seed: seed number to use for repoducibility. Default = 123. set to NULL if use pseudo random.
#'		alpha: a numeric vector to represent the parameter of dirichlet distribution 
#' @param opt.control a list containing parameters for the optimization step:
#'		maxit: maximum number of cycles to interate. Default=100000
#'		sigma: hyper-parameter of the prior if optimizer="paper". Default=2.
#' @param nmf.control a list containing parameters for NMF
#' @param compute.elbo logical to denote if compute ELBO
#' @param ... other arguments passed to nmf
#' @return a list of eta, omega and ELBO to represent the EM results
learn.embedding.nmf <- function(bp,
								  K,
								  cycle = 50,
								  gibbs.control=list(),
								  opt.control=list(),
								  nmf.control=list(),
								  EM.res = NULL,
								  compute.elbo = F,
								  ...){

	#determine if bp.res has updated gibbs result
	if(!bp@control_param$update.gibbs)
		stop("please update gibbs sampling first")
	#dtermine if run using a tumor mode
	if(!is(bp@reference.update, "refTumor"))
		stop("please run BayesPrism using tumor mode by specifying the key parameter")

	#use parameters from BayesPrism run if no change
	gibbs.control.bp <- bp@control_param$gibbs.control
	opt.control.bp <- bp@control_param$opt.control
	
	gibbs.control.bp[names(gibbs.control)] <- gibbs.control
	opt.control.bp[names(opt.control)] <- opt.control

	opt.control <- valid.opt.control(opt.control.bp)
	gibbs.control <- valid.gibbs.control(gibbs.control.bp)

	opt.control$Z_t.min <- NULL		  			
	opt.control$optimizer <- NULL

	nmf.control$".opt" <- paste("vp",min(gibbs.control$n.cores, opt.control$n.cores),sep="")
	nmf.control <- valid.nmf.control(nmf.control)

	#run NMF
	psi_mal <- bp@reference.update@psi_mal
	
	nmf.res <- NMF::nmf(x=t(psi_mal), 
				  		rank = K,
				  		nrun = nmf.control$nrun, 
				  		.opt = nmf.control$".opt", 
				  		seed = nmf.control$seed,
				  		...)
				
	nmf.eta <- t(NMF::basis(nmf.res))

	#run embedding learning
	learn.embedding(bp = bp,
					eta_prior = nmf.eta,
					cycle = cycle,
					gibbs.control = gibbs.control,
					opt.control = opt.control,
					EM.res = EM.res,
					compute.elbo = compute.elbo)
	
}








