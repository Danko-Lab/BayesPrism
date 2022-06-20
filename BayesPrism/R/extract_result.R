


#' function to extract poterior mean of cell type / cell state fraction matrix
#' @param bp a BayesPrism object
#' @param which.theta a character variable to denote to extract first or final Gibbs 
#' @param state.or.type a character variable to denote if cell type or cell state 
#' @return a N*K matrix 
get.fraction <- function(bp,
					   	 which.theta,
					   	 state.or.type){
					   	
	stopifnot(is(bp,"BayesPrism"))
	stopifnot(length(which.theta)==1 & which.theta %in% c("first", "final"))
	stopifnot(length(state.or.type)==1 & state.or.type %in% c("state", "type"))
		
	if(which.theta=="first" & state.or.type=="state")
		return(bp@posterior.initial.cellState@theta) 
	if(which.theta=="first" & state.or.type=="type")
		return(bp@posterior.initial.cellType@theta)  
	if(which.theta=="final"){
		if(state.or.type=="state")
			warning("Warning: only cell type is available for updated Gibbs. Returning cell type info.")
		return(bp@posterior.theta_f@theta) 
	}
	
}

#' function to extract poterior mean of sample specific gene expression profile for the given cell type
#' @param bp a BayesPrism object
#' @param state.or.type a character variable to denote if cell type or cell state 
#' @param cell.name a character variable to denote name of cell type/state 
#' @return a N*G matrix for the cell type/state of interest
get.exp <- function(bp,
					state.or.type,
					cell.name){
					   	
	stopifnot(is(bp,"BayesPrism"))
	stopifnot(length(state.or.type)==1 & state.or.type %in% c("state", "type"))
	
	if(state.or.type=="state")
		return(bp@posterior.initial.cellState@Z[,,cell.name])
	
	if(state.or.type=="type")
		return(bp@posterior.initial.cellType@Z[,,cell.name])
	
}

