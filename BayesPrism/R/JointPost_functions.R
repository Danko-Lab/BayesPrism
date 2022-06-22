#all operations for the jointPost S4 class


#' constructor of the class jointPost
#'
#' @param mixtureID, a character vector to represent identifiyers of bulk samples
#' @param cellType, a character vector to represent identifiyers of cell types
#' @param geneID, a character vector to represent identifiyers of genes
#' @param gibbs.list, list containing the posterior means of gibbs sampling, with each element containing Z_n and theta_n 
#'
#' @return a jointPost S4 object
newJointPost <- function(bulkID,
						 geneID,
						 cellType,
						 gibbs.list){
		
	N <- length(bulkID)
	G <- length(geneID)
	K <- length(cellType)
	
	stopifnot(length(gibbs.list) == N)
						 	
	Z <- array(NA, 
			   dim = c(N,G,K),
			   dimnames=list(bulkID, geneID, cellType))
	
	theta <- matrix(NA, 
					nrow = N, ncol= K, 
					dimnames = list(bulkID, cellType))
	
	theta.cv <- matrix(NA, 
					   nrow = N, ncol= K, 
					   dimnames = list(bulkID, cellType))
	
	for (n in 1:N) {
		Z[n,,] <- gibbs.list[[n]]$Z_n
		theta[n,] <- gibbs.list[[n]]$theta_n
	}
	
	if(!is.null(gibbs.list[[1]]$theta.cv_n)){
		for (n in 1:N) theta.cv[n,] <- gibbs.list[[n]]$theta.cv_n
	}
	else theta.cv <- matrix()
	
	constant <- sum(unlist(lapply(gibbs.list, '[[', "gibbs.constant")))
	 
	new("jointPost", Z = Z, theta = theta, theta.cv = theta.cv, constant = constant)
}


#' constructor of the class jointPost
#'
#' @param mixtureID, a character vector to represent identifiyers of bulk samples
#' @param cellType, a character vector to represent identifiyers of cell types
#' @param geneID, a character vector to represent identifiyers of genes
#' @param gibbs.list, list containing the posterior means of gibbs sampling, with each element containing Z_n and theta_n 
#'
#' @return a jointPost S4 object
newThetaPost <- function(bulkID,
						 cellType,
						 gibbs.list){
		
	N <- length(bulkID)
	K <- length(cellType)
	
	stopifnot(length(gibbs.list) == N)
						 		
	theta <- matrix(NA, 
					nrow = N, ncol= K, 
					dimnames = list(bulkID, cellType))
	
	theta.cv <- matrix(NA, 
					   nrow = N, ncol= K, 
					   dimnames = list(bulkID, cellType))
	
	for (n in 1:N) {
		theta[n,] <- gibbs.list[[n]]$theta_n
		theta.cv[n,] <- gibbs.list[[n]]$theta.cv_n
	}
	
	new("thetaPost", theta = theta, theta.cv = theta.cv)
}


#' function to marginalize K, e.g. cell states within each cell type
#' 
#' @param jointPost.obj a jointPost object
#' @param map a list of the format list(cell.type1=c(cell.stateA, cell.stateB), ...)
#' @return a new jointPost object 
mergeK <- function(jointPost.obj,
				   map){
								
	bulkID <- dimnames(jointPost.obj@Z)[[1]]
	geneID <- dimnames(jointPost.obj@Z)[[2]]
	cellType <- dimnames(jointPost.obj@Z)[[3]]
	cellType.merged <- names(map)
	
	N <- length(bulkID)
	G <- length(geneID)
	K <- length(cellType)
	K_merged <- length(cellType.merged)

	stopifnot(length(unlist(map)) == K)
					 	
	Z <- array(NA, 
			   dim = c(N, G, K_merged),
			   dimnames=list(bulkID, geneID, cellType.merged))
	
	theta <- matrix(NA, 
					nrow = N, ncol= K_merged, 
					dimnames = list(bulkID, cellType.merged))
	
	#merge across cellType.merged
	for(k in 1:K_merged){
		cellType.merged.k <- names(map)[k]
		cellTypes.k <- map[[k]]
		if(length(cellTypes.k)==1) {
			#skipping summation. assign value directly
			Z[,,cellType.merged.k] <- jointPost.obj@Z[,,cellTypes.k, drop=F]
			theta[,cellType.merged.k] <- jointPost.obj@theta[,cellTypes.k, drop=F]
		}
		else {
			#marginalize cellTypes.k
			Z[,, cellType.merged.k] <- rowSums(jointPost.obj@Z[,,cellTypes.k, drop=F], dims=2)
			theta[,cellType.merged.k] <- rowSums(jointPost.obj@theta[,cellTypes.k, drop=F])
		}
	}
	 
	new("jointPost", Z = Z, theta = theta)
	
}


# #' function to marginalize N, e.g. cell states within each cell type
# #' 
# #' @param jointPost.obj a jointPost object
# #' return a matrix of dimension G*K
# mergeN<- function(jointPost.obj){
								
	# geneID <- dimnames(jointPost.obj@Z)[[2]]
	# cellType <- dimnames(jointPost.obj@Z)[[3]]
	
	# G <- length(geneID)
	# K <- length(cellType)
	
	# K_merged <- length(cellType.merged)

	# stopifnot(length(unlist(map)) == K)
					 	
	# Z <- array(NA, 
			   # dim = c(N, G, K_merged),
			   # dimnames=list(bulkID, geneID, cellType.merged))
	

	
	# #merge across cellType.merged
	# for(k in 1:K_merged){
		# cellType.merged.k <- names(map)[k]
		# cellTypes.k <- map[[k]]
		# if(length(cellTypes.k)==1) {
			# #skipping summation. assign value directly
			# Z[,,cellType.merged.k] <- jointPost.obj@Z[,,cellTypes.k]
			# theta[,cellType.merged.k] <- jointPost.obj@theta[,cellTypes.k]
		# }
		# else {
			# #marginalize cellTypes.k
			# Z[,, cellType.merged.k] <- rowSums(jointPost.obj@Z[,,cellTypes.k], dims=2)
			# theta[,cellType.merged.k] <- rowSums(jointPost.obj@theta[,cellTypes.k])
		# }
	# }
	 
	# new("jointPost", Z = Z, theta = theta)
	
# }











 

