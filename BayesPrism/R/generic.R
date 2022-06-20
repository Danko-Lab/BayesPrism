#define generic functions for S4 classes

#' generic show function for an S4 object prism
#' @exportMethod
setMethod("show", "prism",
	function(object) {
		cat("Cell states in each cell type: \n")
		print(object@map)
		cat("\n")
		
		cat("Identifier of the malignant cell type: ", 
			object@key, "\n")
		
		cat("Number of cell states: ", 
			nrow(object@phi_cellState@phi), "\n")
		
		cat("Number of cell types: ",
			nrow(object@phi_cellType@phi), "\n")
		
		cat("Number of mixtures: ",
			nrow(object@mixture), "\n")
		
		cat("Number of genes: ",
			ncol(object@mixture), "\n")
		 
   }
)

#' generic show function for an S4 object BayesPrism
#' @exportMethod
setMethod("show", "BayesPrism",
	function(object) {
		cat("Input prism info: \n")
		show(object@prism)
		cat("\n")
		
		cat("Initial cell type fractions: \n")
		theta.summary <- apply(object@posterior.initial.cellType@theta, 2, summary)	
		print(round(theta.summary,3))
		
		if(!is.null(object@posterior.theta_f)){
			cat("Updated cell type fractions: \n")
			theta.summary <- apply(object@posterior.theta_f@theta, 2, summary)	
			print(round(theta.summary,3))
		}
		 
   }
)


#' generic show function for an S4 object BayesPrismST
#' @exportMethod
setMethod("show", "BayesPrismST",
	function(object) {
		cat("Input prism info: \n")
		show(object@prism)
		cat("\n")
		
		cat("Cell state fractions: \n")
		theta.summary <- apply(object@posterior.cellState@theta, 2, summary)	
		print(round(theta.summary,3))
		
		cat("Cell type fractions: \n")
		theta.summary <- apply(object@posterior.cellType@theta, 2, summary)	
		print(round(theta.summary,3))	 
   }
)







