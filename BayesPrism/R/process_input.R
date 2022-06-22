


#' function that map a character vector of genes to each category on the watchlist
#' @param input.genes a character vector of genes, either gene symbol of ENSEMBL IDs.
#' @param species a character variable to denote if genes are human ("mm") or mouse ("hs")
#' @return a logical matrix of dimension #of genes * # of categories from the watchlist 
#'		to denote whether each gene belongs to a particular category
assign.category <- function(input.genes,
							species=c("mm","hs")){
	
	stopifnot(length(species)==1 & species %in% c("hs","mm"))
	
	#load gene list
	if(species=="hs") gene.list <- read.table(system.file("extdata", "genelist.hs.new.txt", package="BayesPrism"),sep="\t",header=F,stringsAsFactors=F)
	if(species=="mm") gene.list <- read.table(system.file("extdata", "genelist.mm.new.txt", package="BayesPrism"),sep="\t",header=F,stringsAsFactors=F)
	
	
	#detect if EMSEMBLE ID (starts with ENS) or gene symbol is used
	if( sum(substr(input.genes,1,3)=="ENS")> length(input.genes)*0.8  ){
		# use 80% of colnames to avoid sometimes there are manually added gene name such as GFP-XXX.
		#use EMSEMBLE ID
		#strip the "." from ENSXXX.X
		cat("EMSEMBLE IDs detected.\n")
		input.genes.short <- unlist(lapply(input.genes, function(gene.id) strsplit(gene.id,split="\\.")[[1]][1]))
		gene.df <- gene.list[,c(1,2)]
	}
	else{
		#use gene symbols
		cat("Gene symbols detected. Recommend to use EMSEMBLE IDs for more unique mapping.\n")
		input.genes.short <- input.genes
		gene.df <- gene.list[,c(1,3)]
	}
		
	gene.group.matrix <- do.call(cbind.data.frame, lapply(unique(gene.df[,1]), 
						  function(gene.group.i) input.genes.short %in% gene.df[gene.df[,1]== gene.group.i,2]))
	colnames(gene.group.matrix) <- unique(gene.df[,1])
	rownames(gene.group.matrix) <- input.genes
	return(gene.group.matrix)
}

#' function to compute the maximum specificity of each gene (specificity scores)
#'		it provides a fast way of measuring how informative each gene is for deconvolution 
#'		(without the need to go through differential expression analysis).
#'		the idea is from the cell2location tutorial
#' @param input.matrix a collpased expression matrix of dimension K*G
#' @param pseudo.min the desired minimum value after normalization.
#' @returns a vector of length of #of genes to denote the maximum specificity score of each gene
compute.specificity <- function(input.matrix,
								pseudo.min = 1E-8){
	
	ref.ct <- norm.to.one (ref = input.matrix, pseudo.min = pseudo.min)

	#compute gene expression specificity score
	exp.spec <- t(ref.ct) /  colSums(ref.ct)
	max.spec <- apply(exp.spec,1,max)
	
	return(max.spec)
}



#' function to filter genes by gene group (Rb, Mrp, chrM, low.spec.high.exp, etc.)
#' @param input: a M*G count matrix. rownames are cell IDs, while colnames are gene IDs/names.
#' 		 Or a K*G profile matrix (count scale). rownames are cell state/type names; colnames are gene IDs/names.
#' @param input.type A character string. ="count.matrix" if using the raw count matrix of scRNA-seq as input; 
#'		="GEP" if using a gene expression profile, preferably summed raw count(un-normalized), as input.
#' @param species a character variable to denote if genes are human ("mm") or mouse ("hs")
#' @param gene.group a character vector to input gene groups to be removed, must be one or more elements from
#'	c("other_Rb","chrM","chrX","chrY","Rb","Mrp","act","hb","MALAT1","low.spec.high.exp")
#' @param exp.cells genes expressed in number of cells fewer than this will be excluded. Default=1.
#' @return a matrix of the same nrow as input but only keep ncol of number of genes retained
#' @export
cleanup.genes <- function (input,
						   input.type, 
						   species, 
						   gene.group,
						   exp.cells=1){
	
	stopifnot(species %in% c("hs","mm"))
	stopifnot(all(gene.group %in% c("other_Rb","chrM","chrX","chrY","Rb",
									"Mrp","act","hb","MALAT1")))
	
	if(! input.type %in% c("scRNA","count.matrix")) stop("Error: please specify the correct input.type!")
	if(input.type=="GEP"){
		exp.cells <- min(exp.cells,1)
		print("As the input is a collpased GEP, exp.cells is set to min(exp.cells,1)")
	}
	
	category.matrix <- assign.category(input.genes = colnames(input), 
									   species= species)
	category.matrix <- category.matrix[, gene.group, drop=F]
	cat("number of genes filtered in each category: \n")
	print(colSums(category.matrix), na.rm=TRUE)
	
	exclude.idx <- rowSums(category.matrix)>0	
	cat("A total of ", sum(exclude.idx)," genes from", gene.group, " have been excluded","\n")
	input.filtered <- input[, ! exclude.idx] 
	
	if(exp.cells>0) {
		exclude.lowexp.idx <- colSums(input.filtered>0) >= exp.cells
		cat("A total of ", sum(!exclude.lowexp.idx)," gene expressed in fewer than ", exp.cells, " cells have been excluded","\n")
		input.filtered <- input.filtered[, exclude.lowexp.idx]
	}
	else{
		cat("A total of 0 lowly expressed genes have been excluded","\n")
	}
	input.filtered
}




#apply to bulk TCGA
#' function to select genes by gene type (e.g. "protein_coding", "pseudogene", "lincRNA" )
#'	recommend for bulk TCGA data as the annotation was gencode.v22
#'	only works for human genes! For other species please filter manually if needed.
#' @param input: a N*G count matrix. rownames are bulk sample IDs, while colnames are gene IDs/names.
#' @param gene.type a character vector to input gene groups to be retained, must be one or more elements from
#'	c("protein_coding", "pseudogene", "lincRNA")
#' @return a matrix of the same nrow as input but only keep ncol of number of genes retained
select.gene.type <- function(input,
							 gene.type){
	
	stopifnot(all(gene.type %in% c("protein_coding", "pseudogene", "lincRNA")))
	
	input.genes <- colnames(input)
					 	
	gene.tab.path <- system.file("extdata","gencode.v22.broad.category.txt", package="BayesPrism")	
	gene.list <- read.table(gene.tab.path, sep="\t",header=F,stringsAsFactors=F)
	
	#decide if emsembl ID or gene symbol is used
	if( sum(substr(input.genes,1,3)=="ENS")> length(input.genes)*0.8  ){
		# use 80% of colnames to avoid sometimes there are manually added gene name such as GFP-XXX.
		#use EMSEMBLE ID
		#strip the "." from ENSXXX.X
		cat("EMSEMBLE IDs detected.\n")
		input.genes <- unlist(lapply(input.genes, function(gene.id) strsplit(gene.id,split="\\.")[[1]][1]))
		gene.df <- gene.list[match(input.genes, gene.list[,8]),c(8,9)]
	}
	else{
		#use gene symbols
		cat("Gene symbols detected. Recommend to use EMSEMBLE IDs for more unique mapping.\n")
		gene.df <- gene.list[match(input.genes, gene.list[,5]),c(5,9)]
	}
	colnames(gene.df) <- c("gene.name", "category")
	
	selected.gene.idx <- gene.df $category %in% gene.type

	cat("number of genes retained in each category: \n")
	print(table(gene.df[selected.gene.idx,"category"]))
					 	
	input.filtered <- input[, selected.gene.idx]
	
	return(input.filtered)				 	
}



#' function to perform differential expression test 
#' @param sc.dat a cell-by-gene count matrix
#' @param cell.type.labels a character vector to denote cell types of each cell 
#' @param cell.state.labels a character vector to denote cell state of each cell 
#' @param psuedo.count a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
#' @param cell.count.cutoff a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50
#' @return a list of dataframes with p.value and logFC
#' @import scran
#' @import BiocParallel
#' @export 
get.exp.stat <- function(sc.dat,
						 cell.type.labels,
						 cell.state.labels,
						 psuedo.count=0.1,
						 cell.count.cutoff=50,
						 n.cores=1){
		
	ct.to.cst <- unique(cbind(cell.type=cell.type.labels, cell.state=cell.state.labels))
	cst.count.table <- table(cell.state.labels)
	low.count.cst <- names(cst.count.table)[cst.count.table < cell.count.cutoff]

	#normalize ref.dat to prepare input for findMarker	
	lib.size <- rowSums(sc.dat)
	lib.size <- lib.size / median(lib.size)
	dat.tmp <- sc.dat/lib.size
	dat.tmp <- log2(dat.tmp + psuedo.count) - log2(psuedo.count)

	#pairwise t test
	fit.up <- pairwiseTTests(x= t(dat.tmp), 
					   groups= cell.state.labels, 
					   direction="up",
					   BPPARAM = MulticoreParam(n.cores))
	
		
	#filter out comparisons for cell states from the same cell type
	pairs.celltype.first <- ct.to.cst[match(fit.up$pairs$first, ct.to.cst[,"cell.state"]),"cell.type"]
	pairs.celltype.second <- ct.to.cst[match(fit.up$pairs$second, ct.to.cst[,"cell.state"]),"cell.type"]
	
	filter.idx <- pairs.celltype.first != pairs.celltype.second & ! fit.up$pairs$second %in% low.count.cst
	
	fit.up[[1]] <- fit.up[[1]][filter.idx]
	fit.up[[2]] <- fit.up[[2]][filter.idx,]
	
	#get the maxmimum pvalue
	output.up <- combineMarkers(fit.up$statistics, fit.up$pairs, pval.type="all", min.prop=NULL, 
        log.p.in=F, log.p.out=F, full.stats=F, pval.field="p.value", 
        effect.field="logFC", sorted=F)
								     
	all.ct <- unique(ct.to.cst[,"cell.type"])

	#loop over each cell type, and then get statistics over its associated subtypes
	ct.stat.list <- lapply(all.ct,function(ct.i){
		
		#subset on the subtypes associated with celltype i (ct.i)
		cst.i <- ct.to.cst[ct.to.cst[,"cell.type"]==ct.i,"cell.state"]
		output.up.i <- output.up[cst.i]
		
		#take the minimum pvalue over all cst.i
		pval.up.i <- do.call(cbind,lapply(output.up.i, '[', "p.value"))
		pval.up.min.i <- apply(pval.up.i,1,min)
		
		#take the max lfc over the min lfc of cst.i over cst.j in other ct.j (same as the pvalue=min over "all" type)
		lfc.i <- apply(do.call(cbind,lapply(output.up.i, function(output.up.i.j) {
			apply(output.up.i.j[,grepl("logFC",colnames(output.up.i.j)),drop=F],1,min)
		} )),1,max)
		
		data.frame(pval.up.min = pval.up.min.i, 
				   min.lfc = lfc.i)	
	})		

	names(ct.stat.list) <- all.ct
	return(ct.stat.list)
}



#' function to select markers for each cell type
#' @param sc.dat a cell-by-gene count matrix
#' @param stat a list of dataframes with p.value and logFC, outputted from get.exp.stat
#' @param pval.max maximum p value. Default=0.01.
#' @param lfc.min mimimum log2 fold change. Dafult=0.1
#' @return a count matrix with ncol of marker genes.
#' @export 
select.marker <- function(sc.dat,
						  stat,
						  pval.max=0.01,
						  lfc.min=0.1){
	
	cat("number of markers selected for each cell type: \n")
	stat.filtered <- lapply(1:length(stat), function(i){
		stat.i <- stat[[i]]
		stat.i.filtered <- stat.i[stat.i $pval.up.min < pval.max & 
							  	  stat.i $min.lfc > lfc.min,,drop=F]
		cat(names(stat)[i],": ", nrow(stat.i.filtered), "\n")
		return(stat.i.filtered)
	})

	markers.all <- unique(unlist(lapply(stat.filtered,rownames)))
	
	return(sc.dat[, colnames(sc.dat) %in% markers.all])
}




