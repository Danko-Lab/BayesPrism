#' function to normalize expression matrix, s.t. it sum up to one for each row, with the zero entries = pseudo.min
#' if no zero entries, return direct normalization (no pseudo.min returned)
#' @param ref a unnormalized matrix of dimension K*G (with rownames and colnames supplied)
#' @param pseudo.min the desired min values to replace zero after normalization
#' return a normalized matrix of the same dimension
norm.to.one <- function(ref, 
					    pseudo.min){
	
	G <- ncol(ref)
		
	phi <- ref/rowSums(ref) * (1-pseudo.min*G) + pseudo.min
	
	#if the minimum value is greater than zero. simply normalize by total depth
	min.value <- apply(ref,1,min)
	which.row <- min.value>0
	if(any(which.row)){
		#cat("One or more cell types have all genes with non-zero expression. pseudo.min is not applied to these cell types. \n")
		phi[which.row,] <- ref[which.row,,drop=F]/rowSums(ref[which.row,,drop=F])
	}
		
	
	return(phi)
	
}


#' function to summing up reads for each level in labels
#' @param ref a matrix of N*G
#' @param labels a character vector of length N
#' return a matrix of K*G, with K=number of unique levels in labels
collapse <- function(ref, labels){
	
	stopifnot(nrow(ref) == length(labels))
	
	#remove NA in labels
	non.na.idx <- !is.na(labels)
	if(sum(!non.na.idx)>0) print("Warning: NA found in the cell type/state labels. These cells will be excluded!")
	labels <- labels[non.na.idx]
	ref <- ref[non.na.idx,]
	
	labels.uniq <- unique(labels)
	
	ref.collapsed <- do.call(rbind,
							 lapply(labels.uniq,
							 		function(label.i) 
							 			colSums(ref[labels==label.i,,drop=F])
							 		)
							 )
	
	rownames(ref.collapsed) <- labels.uniq
	
	return(ref.collapsed)
}


#' function to validate if input contains negative, non-numeric, NA (stop), or is normalized (warning) or log-transformed (warning)
#' @param input count.matrix, GEP or mixture
validate.input <- function(input){
	#check if referece is non-log transformed
	if(max(input)<=1){
		warning("Warning: input seems to be normalized.") 
	}
	else {
		if(max(input) < 20) 
			warning("Warning: input seems to be log-transformed. Please double check your input. Log transformation should be avoided")
	}
	
	if(min(input)<0)
		stop(" Error: input contains negative values. 
			   Please make sure your input is untransformed raw count. \n")
	
	if (any(!is.finite(input)) | any(is.na(input)) )
		stop(" Error: input contains NaN or NA values.\n")
	
	if(is.null(colnames(input)))
		stop("Error: please specify the colnames of mixture / reference using gene identifiers!")
	
	if(!(is.matrix(input) | inherits(input, "dgCMatrix")))
		stop("Error: the type of mixture and reference need to be matrix!")
	
	NULL
}


#' function to filter bulk outliers (integrated into new.prism)
#' @param mixture the bulk RNA-seq matrix (#of samples * # of genes). 
#' @param outlier.cut & outlier.fraction: Filter genes in mixture whose expression fraction is greater than outlier.cut 
#'		in more than outlier.fraction of bulk data. Removal of outlier genes will ensure that the inference will not be dominated by outliers.
filter.bulk.outlier <- function(mixture,
							  	outlier.cut,
							  	outlier.fraction){
	
	mixture.norm <- mixture / rowSums(mixture)
		
	outlier.idx <- colSums(mixture.norm > outlier.cut) / nrow(mixture.norm) > outlier.fraction
	mixture <- mixture[, !outlier.idx, drop=F]
	cat("Number of outlier genes filtered from mixture =", sum(outlier.idx),"\n")
	
	return(mixture)
}





#' function to contruct a prism object from user provided data
#' @param reference input of scRNA-seq count matrix (#of cells * # of genes) or a GEP matrix (#of cell states/types * # of genes)
#'		rownames(reference) are either cell IDs (if input.type="count.matrix") or names of cell states/types (if input.type="GEP")
#'		colnames(reference) are genes identifiers (e.g. gene symbols or EMSEMBL IDs)
#‘ @param input.type A character string. ="count.matrix" if using the raw count matrix of scRNA-seq as input; 
#'		="GEP" if using a gene expression profile, preferably summed raw count(un-normalized), as input.
#' @param cell.type.labels a character vector to denote cell types of each cell (if input.type="count.matrix") 
#'		or each row of GEP (if input.type="GEP")
#' @param cell.state.labels a character vector to denote cell state of each cell (if input.type="count.matrix") 
#'		or each row of GEP (if input.type="GEP")
#' @param key a string that correponds to malignant cells in cell.type.labels
#'		set to NULL if no malignant cells is to be modeled (all cell types will hence be treated equally)
#' @param pseudo.min the desired min values to replace zero after normalization. Default=1E-8.
#' @param mixture the bulk RNA-seq matrix (#of samples * # of genes). 
#'		colnames(mixture) need to match the annotation of colnames(reference)
#'		new.prism will intersect colnames(reference) and colnames(mixture) and perform deconvolution using the intersected set
#' @param outlier.cut & outlier.fraction: Filter genes in mixture whose expression fraction is greater than outlier.cut (Default=0.01) 
#'		in more than outlier.fraction (Default=0.1) of bulk data. Removal of outlier genes will ensure that the inference will not be dominated by outliers.
#' @param pseudo.min the desired min values to replace zero after normalization
#' @return a prism object
#' @export
new.prism <- function(reference,
					  input.type,
					  cell.type.labels,
					  cell.state.labels,
					  key,
					  mixture,
					  outlier.cut=0.01,
					  outlier.fraction=0.1,
					  pseudo.min=1E-8){
	
	#assign default values
	if(is.null(cell.state.labels)) cell.state.labels <- cell.type.labels
	#force to character
	cell.type.labels <- as.character(cell.type.labels)
	cell.state.labels <- as.character(cell.state.labels)
	
	cat("number of cells in each cell state \n")
	print(sort(table(cell.state.labels)))
	if(min(table(cell.state.labels))<20)
		cat("recommend to have sufficient number of cells in each cell state \n")
	
	if(is.null(key))
		key <- NA_character_
	if(is.na(key))
		cat("No tumor reference is speficied. Reference cell types are treated equally. \n")
	
	#check arguments
	if(length(cell.type.labels) != length(cell.state.labels))
		stop("Error: length of cell.type.labels and cell.state.labels do not match!")
	if(length(cell.type.labels) != nrow(reference))
		stop("Error: length of cell.type.labels and nrow(reference) do not match!")	
	#creat mapping between cell type and phenotype (cell type is a superset of phenotype)
	type.to.state.mat <- unique(cbind(cell.type.labels, cell.state.labels))
	if (max(table(type.to.state.mat[,"cell.state.labels"]))>1) 
		stop("Error: one or more cell states belong to multiple cell types!")
	if (length(unique(cell.type.labels)) > length(unique(cell.state.labels))) 
		stop("Error: more cell types than states!")
	
	if(is.data.frame(mixture)) mixture <- as.matrix(mixture)
	if(is.data.frame(reference)) reference <- as.matrix(reference)	
	
	if(is.null(dim(mixture))) mixture <- matrix(mixture, nrow=1, dimnames=list("mixture-1", names(mixture)))
	
	if(is.null(rownames(mixture))) rownames(mixture) <- paste("mixture", 1:nrow(mixture), sep="-")
	
	#check input
	validate.input(reference)
	validate.input(mixture)
	
	#filter outliers from bulk
	mixture <- filter.bulk.outlier (mixture=mixture,
							  		outlier.cut=outlier.cut,
							  		outlier.fraction=outlier.fraction)
	
	#remove zero genes in reference
	reference <- reference[,colSums(reference)>0]
	
	#check gene annotation between reference and mixture
	gene.shared <- intersect(colnames(reference), colnames(mixture))
	if(length(gene.shared)==0)
		stop("Error: gene names of reference and mixture do not match!") 
	if(length(gene.shared)<100)
		warning("Warning: very few gene from reference and mixture match! Please double check your gene names.")

	#collapse reference
	ref.cs <- collapse(ref = reference, labels = cell.state.labels)
	ref.ct <- collapse(ref = reference, labels = cell.type.labels)

	#align reference and mixture
	cat("Aligning reference and mixture... \n")
	ref.ct <- ref.ct[, gene.shared, drop=F]
	ref.cs <- ref.cs[, gene.shared, drop=F]
	mixture <- mixture[, gene.shared, drop=F]
	
	#normalize
	cat("Normalizing reference... \n")
	ref.cs <- norm.to.one (ref = ref.cs, pseudo.min = pseudo.min)
	ref.ct <- norm.to.one (ref = ref.ct, pseudo.min = pseudo.min)
	
	#build the mapping dictionary(list)
	map <- lapply(rownames(ref.ct),function(t)unique(cell.state.labels[cell.type.labels==t]))
	names(map) <- rownames(ref.ct)
	
	#contruct a prism object
	new("prism", 
		phi_cellState = new("refPhi", phi = ref.cs, pseudo.min = pseudo.min),
		phi_cellType = new("refPhi", phi = ref.ct, pseudo.min = pseudo.min),
		map = map,
		key = key,
		mixture = mixture)
	
}




