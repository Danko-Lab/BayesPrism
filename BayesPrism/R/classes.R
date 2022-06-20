#define classes for Gibbs sampling


#' S4 class to store (non-malignant) reference matrix phi or psi (if no malignant cells)
#' 
#' @slot phi a matrix of dimension K*G, 
#'		rownames are cell state/type names; colnames are gene IDs/names
#' @slot pseudo.min the desired minimum value used to normalize phi.
setClass("refPhi",
         slots = c(
           phi = "matrix",
           pseudo.min ="numeric"
         ),
         prototype = list(
           phi = matrix(),
           pseudo.min = NA_real_
         ),
         validity = function(object){
         	errors <- character()
         	
         	if(length(object@pseudo.min) != 1 | !is.numeric(object@pseudo.min)){
  				msg <- paste("invalid pseudo.min")
  				errors <- c(errors, msg)
  			}
  			
  			phi.min <- min(object@phi)
         	if(phi.min < 0){
				msg <- paste("reference contain negative values")
				errors <- c(errors, msg)
         	}
         	if(!is.na(object@pseudo.min) & phi.min != object@pseudo.min)
         		warning("Warning: pseudo.min does not match min(phi)")
         	if (length(errors) == 0) TRUE else errors
         }
)


#' S4 class to store reference matrix psi for malignant cells
#' 
#' @slot psi_mal a matrix of dimension N*G, to denote updated sample-specific profiles (for malignant cells) 
#'		rownames are bulk sample IDs; colnames are gene IDs/names
#' @slot psi_env a matrix of dimension (K-1)*G, to denote updated shared profiles (for non-malignant cells) 
#'		rownames are non-malignant cell types; colnames are gene IDs/names
#' @slot key a character variable to denote the names for malignant cells
#' @slot pseudo.min the desired minimum value used to normalize phi.
setClass("refTumor",
         slots = c(
           psi_mal = "matrix",
           psi_env = "matrix",
           key = "character",
           pseudo.min ="numeric"
         ),
         prototype = list(
           psi_mal = matrix(),
           psi_env = matrix(),
           key = NA_character_,
           pseudo.min = NA_real_
         ),
         validity = function(object){
         	errors <- character()
  			
  			if(length(object@key) != 1 | !is.character(object@key)){
  				msg <- paste("invalid key")
  				errors <- c(errors, msg)
  			}
  			  					
  			psi_mal_genes <- colnames(object@psi_mal)
  			psi_env_genes <- colnames(object@psi_env)
  			if(!identical(psi_mal_genes, psi_env_genes)){
  				msg <- paste("Gene names of psi_mal and psi_env do not match")
  				errors <- c(errors, msg)
  			}
  				
  			if(length(object@pseudo.min) != 1 | !is.numeric(object@pseudo.min)){
  				msg <- paste("invalid pseudo.min")
  				errors <- c(errors, msg)
  			}
  			
  			phi.min <- min(object@psi_mal, object@psi_env)
         	if(phi.min < 0){
				msg <- paste("reference contain negative values")
				errors <- c(errors, msg)
         	}
         	if(!is.na(object@pseudo.min) & phi.min != object@pseudo.min)
         		warning("Warning: pseudo.min does not match min(phi)")
         	if (length(errors) == 0) TRUE else errors
         }

)


#' virtual S4 class
setClassUnion("reference", c("refPhi", "refTumor", "NULL"))


#' S4 class to represent the input to run run.prism. Contains the reference and mixture. 
#' 
#' @slot phi_cellState a matrix of a matrix of dimension K(of cell state)*G, 
#'		rownames are cell state names; colnames are gene IDs/names
#' @slot phi_cellType a matrix of a matrix of dimension K(of cell state)*G, 
#'		rownames are cell state names; colnames are gene IDs/names
#' @slot map a list to store the correspondence between cell states and cell types.
#"		of the format list(cell.type1=c(cell.stateA, cell.stateB), ...)
#' @slot key a character variable to denote the names for malignant cells
#' @slot pseudo.min the desired minimum value used to normalize phi.
#' @slot mixture the bulk RNA-seq matrix (#of samples * # of genes). 
#'		rownames are bulk sample IDs; colnames are gene IDs/names
#' @export
setClass("prism",
         slots = c(
           phi_cellState = "refPhi",
           phi_cellType = "refPhi",
           map = "list",
           key = "character",
           mixture = "matrix"
         ),
         prototype = list(
           phi_cellState = new("refPhi"),
           phi_cellType = new("refPhi"),
           map = list(),
           key = NA_character_,
           mixture = matrix()
         ),
         validity = function(object){
         	errors <- character()
  			
  			mykey <- object@key
  			mymap <- object@map
  			
  			if(!is.na(mykey) & 
  			   (length(mykey) != 1 | !is.character(mykey) | ! mykey %in% names(mymap))){
  				msg <- paste("invalid key")
  				errors <- c(errors, msg)
  			}
  			
  			if(!is.na(mykey) & min(object@phi_cellState@pseudo.min, object@phi_cellType@pseudo.min)==0){
  				msg <- paste("psuedo.min needs to be strictly positive when running under tumor mode.")
  				errors <- c(errors, msg)
  			}	
  			  					
  			cs_genes <- colnames(object@phi_cellState@phi)
  			ct_genes <- colnames(object@phi_cellType@phi)
  			bk_genes <- colnames(object@mixture)
  			if(!identical(cs_genes, ct_genes) | !identical(cs_genes, bk_genes)){
  				msg <- paste("Gene names do not match")
  				errors <- c(errors, msg)
  			}
  				
			if(!identical(rownames(object@phi_cellType@phi), names(mymap))){
  				msg <- paste("cell types between map and phi_cellType do not match")
  				errors <- c(errors, msg)
  			}
			
			if(!setequal(rownames(object@phi_cellState@phi), unlist(mymap))){
  				msg <- paste("cell states between map and phi_cellState do not match")
  				errors <- c(errors, msg)
  			}
						
         	if (length(errors) == 0) TRUE else errors
         }
)



#' An S4 class to represent the posterior mean of theta, i.e. the output, of Gibbs sampling p( theta | X, phi; alpha),
#'
#' @slot Z an array of dimension N*G*K to denote sample-specific cell type/state-specific gene expression profile
#'		dimnames are bulk IDs, cell type/state names, and gene IDs/names

#' @slot theta a matrix of dimension N*K to denote the cell type fraction in each bulk
#'		rownames are bulk IDs; colnames are gene cell type/state
#'
#' @export 
setClass("thetaPost",
         slots = c(
           theta = "matrix",
           theta.cv = "matrix"
         ),
         prototype = list(
           theta = matrix(),
           theta.cv = matrix()
         ),
         validity = function(object){
         	errors <- character()
  			
  			if(!identical(dimnames(object@theta), dimnames(object@theta.cv))){
  				msg <- paste("dimnames do not match")
  				errors <- c(errors, msg)
  			}
  			  			  											
         	if (length(errors) == 0) TRUE else errors
         }
)




#' An S4 class to represent the joint posterior mean, i.e. the output, of Gibbs sampling p( Z, theta | X, phi; alpha),
#'
#' @slot Z an array of dimension N*G*K to denote sample-specific cell type/state-specific gene expression profile
#'		dimnames are bulk IDs, cell type/state names, and gene IDs/names

#' @slot theta a matrix of dimension N*K to denote the cell type fraction in each bulk
#'		rownames are bulk IDs; colnames are gene cell type/state
#'
#' @export 
setClass("jointPost",
         slots = c(
           Z = "array",
           theta = "matrix",
           theta.cv = "matrix",
           constant = "numeric"
         ),
         prototype = list(
           Z = array(),
           theta = matrix(),
           theta.cv = matrix(),
           constant = numeric()
         ),
         validity = function(object){
         	errors <- character()
  			
  			if(!identical(dimnames(object@Z)[[1]], rownames(object@theta))){
  				msg <- paste("sample ID names do not match")
  				errors <- c(errors, msg)
  			}
  			
  			if(!identical(dimnames(object@Z)[[3]], colnames(object@theta))){
  				msg <- paste("cell type/state names do not match")
  				errors <- c(errors, msg)
  			}
  			  											
         	if (length(errors) == 0) TRUE else errors
         }
)







#' An S4 class to represent the input, X, phi; alpha and controls of Gibbs sampling 
#'
#' @slot phi an array of dimension K*G to denote reference matrix
#' @slot X a matrix of dimension N*G to denote bulk matrix
#' @slot theta a matrix of dimension K*N to denote the cell type fraction in each bulk
#' @slot alpha a numeric value to denote the symmetric Dirichlet prior, fixed to 1E-8 
#' @slot gibbs.control a list containing parameters of the Gibbs sampler
setClass("gibbsSampler",
         slots = c(
           reference = "reference",
           X = "matrix",
           gibbs.control = "list"
         ),
         prototype = list(
           reference = NULL,
           X = matrix(),
           gibbs.control = list(chain.length = NULL,
           					    burn.in = NULL,
           					    thinning = NULL,
           					    n.cores = NULL,
           					    seed = NULL,
           					    alpha = NULL)
         ),
         validity = function(object){
         	errors <- character()
  			
  			if(is(object@reference, "refPhi"))
  				ref.gene <- colnames(object@reference@phi)
  			if(is(object@reference, "refTumor"))
  				ref.gene <- colnames(object@reference@psi_mal)
  			
  			if(!identical(ref.gene, colnames(object@X))){
  				msg <- paste("Gene names do not match")
  				errors <- c(errors, msg)
  			}
  			  											
         	if (length(errors) == 0) TRUE else errors
         }
)


setClassUnion("theta_f", c("thetaPost", "NULL"))


#' An S4 class to represent the output of run.prism
#'
#' @slot input.initial.cellState an S4 object of class gibbsSampler 
#'		to store the input of initial Gibbs sampling over cell states 
#' @slot posterior.initial.cellState an S4 object of class jointPost 
#'		to store the results of initial Gibbs sampling over cell states
#' @slot posterior.initial.cellType an S4 object of class jointPost 
#'		to store the results over cell types merged from posterior.initial.cellState
#' @slot reference.update an S4 object of class reference 
#'		to store the updated reference
#' @slot posterior.theta_f matrix to represent the cell type fraction from the final Gibbs sampling
#' @slot control_param a list to store all 
#'
#' @export
setClass("BayesPrism",
         slots = c(
           prism = "prism",
           posterior.initial.cellState = "jointPost",
           posterior.initial.cellType = "jointPost",
           reference.update = "reference",
           posterior.theta_f = "theta_f",
           control_param = "list"
         ),
         prototype = list(
           prism = new("prism"),
           posterior.initial.cellState = new("jointPost"),
           posterior.initial.cellType = new("jointPost"),
           reference.update = NULL,
           posterior.theta_f = NULL,
           control_param = list(gibbs.control = list(chain.length = NULL,
           					    					  burn.in = NULL,
           					    					  thinning = NULL,
           					    					  n.cores = NULL,
           					    					  seed = NULL,
           					    					  alpha=NULL),
					  				opt.control = list(trace = NULL, 
					  								    maxit = NULL,
					  								    optimizer = NA_character_, 
					  								    n.cores = NULL),
					  				map = list(),
					  				key = NA_character_,
					  				update.gibbs = logical())
         )
)




#' An S4 class to represent the output of run.prism.st 
#'
#' @slot input.initial.cellState an S4 object of class gibbsSampler 
#'		to store the input of initial Gibbs sampling over cell states 
#' @slot posterior.initial.cellState an S4 object of class jointPost 
#'		to store the results of initial Gibbs sampling over cell states
#' @slot posterior.initial.cellType an S4 object of class jointPost 
#'		to store the results over cell types merged from posterior.initial.cellState
#' @slot reference.update an S4 object of class reference 
#'		to store the updated reference
#' @slot posterior.theta_f matrix to represent the cell type fraction from the final Gibbs sampling
#' @slot control_param a list to store all 
#'
#' @export
setClass("BayesPrismST",
         slots = c(
           prism = "prism",
           posterior.cellState = "jointPost",
           posterior.cellType = "jointPost",
           reference.update = "reference",
           control_param = "list"
         ),
         prototype = list(
           prism = new("prism"),
           posterior.cellState = new("jointPost"),
           posterior.cellType = new("jointPost"),
           control_param = list(gibbs.control = list(chain.length = NULL,
           					    					  burn.in = NULL,
           					    					  thinning = NULL,
           					    					  n.cores = NULL,
           					    					  seed = NULL,
           					    					  alpha=NULL),
					  				opt.control = list(trace = NULL, 
					  								    maxit = NULL,
					  								    optimizer = NA_character_, 
					  								    n.cores = NULL),
					  				map = list())
         )
)

