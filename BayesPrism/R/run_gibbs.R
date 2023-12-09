#' function to return the index set of MCMC chain
#' @param gibbs.control a list containing parameters of the Gibbs sampler 
#' @return a numeric vector of the index of samples to be retained from MCMC chain 
get.gibbs.idx <- function(gibbs.control){
	chain.length <- gibbs.control$chain.length
	burn.in <- gibbs.control$burn.in
	thinning <- gibbs.control$thinning
	all.idx <- 1: chain.length
	burned.idx <-  all.idx[-(1: burn.in)]
	thinned.idx <- burned.idx[seq(1,length(burned.idx), by= thinning)]
	thinned.idx
}



#' function that generates one sample from Dirichlet distribution
#' param alpha a numeric vector to represent the parameter of dirichlet distribution
#' return a numeric vector (one sample from Dirichlet distribution)
rdirichlet <- function(alpha){
    l <- length(alpha)
    x <- rgamma(length(alpha), alpha)
    return(x/sum(x))
}


#' function to run Gibbs sampling for Z and theta on each bulk
#' @param X_n a numeric vector of reads count of the nth bulk sample
#' @param phi an array of dimension K*G to denote reference matrix
#' @param alpha a numeric value to denote the symmetric Dirichlet prior 
#' @param gibbs.idx a numeric vector to denote the index of samples to be retained from MCMC chain
#' @param compute.elbo a logical variable to denote if compute ELBO. Default=FALSE.
#' return a list containing the posterior mean of Z_n and theta_n
sample.Z.theta_n <- function(X_n, 
				        	 	 phi,
				        	 	alpha,
				        	 	gibbs.idx,
				        	 	compute.elbo=FALSE){
		
	G <- ncol(phi)
	K <- nrow(phi)
		
	theta_n.i <- rep(1/K, K)
	Z_n.i <- array(NA,c(G,K))
	
	Z_n.sum <- array(0,c(G,K))
	theta_n.sum <- rep(0, K)
	theta_n2.sum <- rep(0, K)
	
	#variable for computing ELBO
	multinom.coef <- 0

	for(i in 1:max(gibbs.idx)){
	
		#sample Z for patient n
		prob.mat <- phi * theta_n.i #multiply by each column

		for (g in 1:G) Z_n.i[g,] <- rmultinom(n = 1, 
										   size = X_n[g], 
										   prob = prob.mat[,g])				
		# sample theta for patient n
		Z_nk.i <- colSums(Z_n.i) #total count for each cell type
		theta_n.i <- rdirichlet(alpha = Z_nk.i + alpha)
				
		if(i %in% gibbs.idx) {
			#collect sample and compute posterior sum
			Z_n.sum <- Z_n.sum + Z_n.i
			theta_n.sum <- theta_n.sum + theta_n.i
			theta_n2.sum <- theta_n2.sum + theta_n.i^2
			
			if(compute.elbo){
				# multinom.coef.i <- sum((alpha-1) * log(theta_n.i), na.rm=TRUE) + 
								   # sum(Z_nk.i * log(theta_n.i), na.rm=TRUE) - 
								   # sum(lfactorial(Z_nk.i)) - sum(lfactorial(Z_n.i))
				multinom.coef.i <- sum(lfactorial(Z_nk.i)) - sum(lfactorial(Z_n.i))
				multinom.coef <- multinom.coef + multinom.coef.i
				
			}	
		}
		
		if((i %% 50) == 0) gc()
	}
	
	samples.size <- length(gibbs.idx)
	#gibbs.constant <- multinom.coef + z.logtheta + theta.dirichlet.alpha
	
	Z_n <- Z_n.sum / samples.size
	theta_n <- theta_n.sum / samples.size
	theta.cv_n <- sqrt(theta_n2.sum / samples.size - (theta_n^2)) / theta_n
	gibbs.constant <- multinom.coef / samples.size
		
	return(list(Z_n = Z_n, 
			    theta_n = theta_n,
			    theta.cv_n = theta.cv_n,
			    gibbs.constant = gibbs.constant))
}



#' function to run Gibbs sampling for only theta on each bulk (updated Gibbs)
#' @param X_n a numeric vector of reads count of the nth bulk sample
#' @param phi an array of dimension K*G to denote reference matrix of Nth sample
#' @param alpha a numeric value to denote the symmetric Dirichlet prior 
#' @param gibbs.idx a numeric vector to denote the index of samples to be retained from MCMC chain 
#' return a vector containing the posterior mean of theta_n
sample.theta_n <- function(X_n, 
				           phi,
				           alpha,
				           gibbs.idx){
		
	G <- ncol(phi)
	K <- nrow(phi)
		
	theta_n.i <- rep(1/K, K)
	Z_n.i <- array(NA,c(G,K))
	
	theta_n.sum <- rep(0, K)
	theta_n2.sum <- rep(0, K)

	for(i in 1:max(gibbs.idx)){
	
		#sample Z for patient n
		prob.mat <- phi * theta_n.i #multiply by each column
		
		for (g in 1:G) Z_n.i[g,] <- rmultinom(n = 1, 
										   size = X_n[g], 
										   prob = prob.mat[,g])				
		# sample theta for patient n
		theta_n.i <- rdirichlet(alpha = colSums(Z_n.i) + alpha)
				
		if(i %in% gibbs.idx) {
			theta_n.sum <- theta_n.sum + theta_n.i
			theta_n2.sum <- theta_n2.sum + theta_n.i^2
		}
		if((i %% 50) == 0) gc()
	}
	
	samples.size <- length(gibbs.idx)
	
	theta_n <- theta_n.sum / samples.size
	theta.cv_n <- sqrt(theta_n2.sum / samples.size - (theta_n^2)) / theta_n
	
	return(list(theta_n = theta_n,
			    theta.cv_n = theta.cv_n))

}


#' function to convert seconds to day/hour/min/sec. used to report the estimated the run time of gibbs sampler
#' adpated from https://stackoverflow.com/questions/27312292/convert-seconds-to-days-hoursminutesseconds
#' @param x numeric value in the unit of seconds
my_seconds_to_period <- function(x) {
	days = round(x %/% (60 * 60 * 24))
	hours = round((x - days*60*60*24) %/% (60 * 60))
	minutes = round((x - days*60*60*24 - hours*60*60) %/% 60) +1 
	days_str = ifelse(days == 0, "", paste0(days, "days "))
	hours_str = ifelse((hours == 0 & days == 0), "", paste0(hours, "hrs "))
	minutes_str = ifelse((minutes == 0 & days == 0 & hours == 0), "", paste0(minutes, "mins"))
	final_str = paste0(days_str, hours_str, minutes_str)
	return(final_str)
}

#‘ function used to estimate run time of Gibbs sampling (print to console)
#' @param gibbsSampler.obj a gibbsSampler object
#' @param final a logical variable denote whether report only updated theta_f (final =TRUE)
#' @param chain.length length of MCMC chain used to simulate. Default=50.
estimate.gibbs.time <- function(gibbsSampler.obj,
								final, 
								chain.length=50){
	
	ref <- gibbsSampler.obj@reference
	X <- gibbsSampler.obj@X
	gibbs.control <- gibbsSampler.obj@gibbs.control
	
	ptm <- proc.time()
	
	if(!final){
		#initial Gibbs
		sample.Z.theta_n (X_n = X[1,], 
								  phi = ref@phi, 
								  alpha = gibbs.control$alpha,
								  gibbs.idx = get.gibbs.idx(
				     		 	  	list(chain.length = chain.length, 
				     					 burn.in = chain.length*gibbs.control$burn.in/gibbs.control$chain.length,
				     					 thinning = gibbs.control$thinning)
				     		 	  )
						   )
	}
	else{
		#updated Gibbs
		if(is(ref,"refPhi")){
			sample.theta_n (X_n = X[1,], 
							phi = ref@phi, 
							alpha = gibbs.control$alpha,
							gibbs.idx = get.gibbs.idx(
				     			list(chain.length = chain.length, 
				     				burn.in = chain.length*gibbs.control$burn.in/gibbs.control$chain.length,
				     				thinning = gibbs.control$thinning)
				     		 	 )
				    			)
		}
		if(is(ref,"refTumor")){
		 	sample.theta_n (X_n = X[1,], 
							phi = rbind(ref@psi_mal[1,], ref@psi_env), 
							alpha = gibbs.control$alpha,
							gibbs.idx = get.gibbs.idx(
				     			list(chain.length = chain.length, 
				     				 burn.in = chain.length*gibbs.control$burn.in/gibbs.control$chain.length,
				     				 thinning = gibbs.control$thinning)
				     		 	)
				    			)	
		}	
	}
	
	total.time <- proc.time() - ptm
	total.time <- as.numeric(total.time["elapsed"])
	
	# seems to underestimate when transferring data comsumes large amount of time (spatial data)
	estimated.time <- gibbs.control $chain.length / chain.length * total.time * ceiling(nrow(X) / gibbs.control$n.cores) *2 		
	current.time <- Sys.time()
	
	cat("Current time: ", as.character(current.time), "\n")
	cat("Estimated time to complete: ", my_seconds_to_period(estimated.time), "\n")
	cat("Estimated finishing time: ", as.character(current.time + estimated.time), "\n")

}






#' function to run initial Gibbs sampling if reference is of the class refPhi
#'
#' @param gibbsSampler.obj, a gibbsSampler object
#' @import snowfall
#' @return a jointPost object with Z and theta entries 
run.gibbs.refPhi.ini <- function(gibbsSampler.obj,
							 	   compute.elbo){
	
	phi <- gibbsSampler.obj@reference@phi
	X <- gibbsSampler.obj@X
	gibbs.control <- gibbsSampler.obj@gibbs.control
	alpha <- gibbs.control$alpha
	
	gibbs.idx <- get.gibbs.idx(gibbs.control)
	seed <- gibbs.control$seed

	sample.Z.theta_n <- BayesPrism:::sample.Z.theta_n
		
	cat("Start run... \n")
	
	if(gibbs.control$n.cores>1){	
		#parallel using snowfall	
		sfInit(parallel = TRUE, cpus = gibbs.control$n.cores, type = "SOCK" )
					
		cpu.fun <- function(n) {
			if(!is.null(seed)) set.seed(seed)
			#load nth mixture from disk
			file.name.X_n <- paste(tmp.dir, "/mixture_",n,".rdata",sep="")
			load(file.name.X_n)
			
			sample.Z.theta_n (X_n = X_n, 
							  phi = phi, 
							  alpha = alpha, 
							  gibbs.idx = gibbs.idx, 
							  compute.elbo = compute.elbo)		
		}
		tmp.dir <- tempdir(check=TRUE)
		sfExport("phi", "alpha", "gibbs.idx", "seed", 
				 	"compute.elbo", "sample.Z.theta_n","tmp.dir")
		environment(cpu.fun) <- globalenv()
		gibbs.list <- sfLapply( 1:nrow(X), cpu.fun)
		sfStop()
	}
	else{
		#single thread
		cpu.fun <- function(n) {
				if(!is.null(seed)) set.seed(seed)
				cat(n," ")
				sample.Z.theta_n (X_n = X[n,], phi = phi, alpha = alpha, 
								  gibbs.idx = gibbs.idx, compute.elbo = compute.elbo)
		}
		gibbs.list <- lapply( 1:nrow(X), cpu.fun)
		cat("\n")
	}
	
	jointPost <- newJointPost(bulkID = rownames(X),
							  geneID = colnames(X), 
							  cellType = rownames(phi),
							  gibbs.list = gibbs.list )
	return(jointPost)
}






#' function to run final sampling if reference is of the class refPhi
#'
#' @param gibbsSampler.obj, a gibbsSampler object
#' @import snowfall
#' @return a theta_f matrix 
run.gibbs.refPhi.final <- function(gibbsSampler.obj,
							 	   compute.elbo){
	
	phi <- gibbsSampler.obj@reference@phi
	X <- gibbsSampler.obj@X
	gibbs.control <- gibbsSampler.obj@gibbs.control
	alpha <- gibbs.control$alpha
	
	gibbs.idx <- get.gibbs.idx(gibbs.control)
	seed <- gibbs.control$seed

	sample.theta_n <- BayesPrism:::sample.theta_n
		
	cat("Start run... \n")
	
	if(gibbs.control$n.cores>1){	
		#parallel using snowfall	
		sfInit(parallel = TRUE, cpus = gibbs.control$n.cores, type = "SOCK" )
		
		cpu.fun <- function(n) {
			if(!is.null(seed)) set.seed(seed)
			#load nth mixture from disk
			file.name.X_n <- paste(tmp.dir, "/mixture_",n,".rdata",sep="")
			load(file.name.X_n)
			
			sample.theta_n (X_n = X_n, 
								phi = phi, 
								alpha = alpha, 
								gibbs.idx = gibbs.idx)		
		}
		tmp.dir <- tempdir()
		sfExport("phi", "alpha", "gibbs.idx", "seed", 
				 "compute.elbo", "sample.theta_n", "tmp.dir")
		environment(cpu.fun) <- globalenv()
		gibbs.list <- sfLapply( 1:nrow(X), cpu.fun)
		sfStop() 
	}
	else{
		#single thread
		cpu.fun <- function(n) {
			if(!is.null(seed)) set.seed(seed)
			cat(n," ")
			sample.theta_n (X_n = X[n,], phi = phi, alpha = alpha, gibbs.idx = gibbs.idx)
		}
		gibbs.list <- lapply( 1:nrow(X), cpu.fun)
		cat("\n")
	}
	thetaPost <- newThetaPost (bulkID = rownames(X),
						 			   cellType = rownames(phi),
						 			   gibbs.list = gibbs.list)
	return(thetaPost)
}


#' function to run Gibbs sampling if reference is of the class refPhi
#'
#' @param gibbsSampler.obj, a gibbsSampler object
#' @import snowfall
#' @return a theta_f matrix 
run.gibbs.refTumor <- function(gibbsSampler.obj){
	
	psi_mal <- gibbsSampler.obj@reference@psi_mal
	psi_env <- gibbsSampler.obj@reference@psi_env
	key <- gibbsSampler.obj@reference@key
	
	X <- gibbsSampler.obj@X
	gibbs.control <- gibbsSampler.obj@gibbs.control
	alpha <- gibbs.control$alpha

	gibbs.idx <- get.gibbs.idx(gibbs.control)
	seed <- gibbs.control$seed
	
	sample.theta_n <- BayesPrism:::sample.theta_n
	

	tmp.dir <- tempdir()
	for(n in 1:nrow(psi_mal)) {
		psi_mal_n <- psi_mal[n,]
		file.name <- paste(tmp.dir, "/psi_mal_",n,".rdata",sep="")
		save(psi_mal_n, file= file.name)
	}
	
	cat("Start run... \n")
	
	cpu.fun <- function(n) {
		if(!is.null(seed)) set.seed(seed)
			
		#load reference of malignant cells and nth mixture from disk
		file.name.psi_mal_n <- paste(tmp.dir, "/psi_mal_",n,".rdata",sep="")
		load(file.name.psi_mal_n)
		file.name.X_n <- paste(tmp.dir, "/mixture_",n,".rdata",sep="")
		load(file.name.X_n)
			
		sample.theta_n (X_n = X_n, 
						phi = rbind(psi_mal_n, psi_env), 
						alpha = alpha,
						gibbs.idx = gibbs.idx)				    				
	}
	sfInit(parallel = TRUE, cpus = gibbs.control$n.cores, type = "SOCK" )
	tmp.dir <- tempdir()
	sfExport("psi_env", "X", "alpha", "gibbs.idx", "seed", "sample.theta_n", "tmp.dir")
	
	environment(cpu.fun) <- globalenv()
	gibbs.list <- sfLapply( 1:nrow(X), cpu.fun)
	sfStop()

	thetaPost <- newThetaPost (bulkID = rownames(X),
						 	   cellType = c(key, rownames(psi_env)),
						 	   gibbs.list = gibbs.list)
	
	unlink(tmp.dir, recursive = TRUE)
	return(thetaPost)
}



#' function to run Gibbs sampling 
#'
#' @param gibbsSampler.obj, a gibbsSampler object
#' @param final a logical variable denote whether report only updated theta_f (final=TRUE)
#' @param if.estimate a logical variable denote whether estimate the run time. Default=TRUE
#' @import snowfall
#' @return a gibbsSampler object with Z and theta entries, if final=FALSE,  or theta_f matrix if final=TRUE
run.gibbs <- function(gibbsSampler.obj, 
					  final,
					  if.estimate=TRUE,
					  compute.elbo=FALSE
					  ){
	
	if(final) cat("Run Gibbs sampling using updated reference ... \n")
	else cat("Run Gibbs sampling... \n")

	if(if.estimate)
		estimate.gibbs.time(gibbsSampler.obj = gibbsSampler.obj, 
							final = final)
	
	if(is(gibbsSampler.obj@reference,"refPhi")) {
		if(!final) return(run.gibbs.refPhi.ini(gibbsSampler.obj = gibbsSampler.obj, 
											   compute.elbo = compute.elbo))
		else return(run.gibbs.refPhi.final(gibbsSampler.obj = gibbsSampler.obj, 
										   compute.elbo = compute.elbo))
	}
	if(is(gibbsSampler.obj@reference,"refTumor")) {
		return(run.gibbs.refTumor(gibbsSampler.obj = gibbsSampler.obj))	
	}
}






