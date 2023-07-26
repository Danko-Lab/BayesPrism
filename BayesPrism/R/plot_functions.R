#' function to visualize outliers in scRNA-seq data or reference matrix 
#'	by plotting specificty score vs log (mean normalized expression) 
#' 	This function will also recommend outliers (high expression and low specificty) for removal
#' 	Outliers are determined using ribosomal protein coding genes. 
#'	Users are recommended to add back ribosomal genes if they have been previosuly removed during preprocessing.
#' @param input: a M*G count matrix. rownames are cell IDs, while colnames are gene IDs/names.
#' 		 Or a K*G profile matrix (count scale). rownames are cell state/type names; colnames are gene IDs/names.
#' @param cell.type.labels a character vector to denote cell types of each cell (if input is a count.matrix) 
#'		or each row of GEP (if input is a "GEP")
#' @param species a character variable to denote if genes are human ("hs") or mouse ("mm")

#' @param pdf.prefix a character variable for the prefix of the pdf name. Default=NULL (plot directly).
#' @return a logical matrix of dimension #of genes * # of categories from the watchlist 
#'	to denote whether each gene belongs to a particular category
#'	an additional column low.spec.high.exp, determined based on ribosomal protein genes as described previously,
#'	will be appended to the logical matrix.  	
#' @export 
plot.scRNA.outlier <- function(input,
							   cell.type.labels,
							   species,
							   pdf.prefix=NULL,
							   return.raw=FALSE){
	
	input <- collapse(ref = input, labels = cell.type.labels)
	ref.ct <- norm.to.one (ref = input, pseudo.min = 1E-8)

	#compute mean expression
	exp.mean <- colMeans(ref.ct)
	exp.mean.log <- log(exp.mean)

	#compute gene expression specificity score
	max.spec <- compute.specificity(input)

	category.matrix <- assign.category(input.genes = colnames(ref.ct), 
									   species= species)
	
	#make plots	
	if(!is.null(pdf.prefix)) pdf(paste(pdf.prefix,"_scRNA_outlier.pdf",sep=""),pointsize=8,useDingbats=FALSE )
	
	cex=0.5
	plot(y= max.spec, x= exp.mean.log, cex=0.2, pch=16, col="grey", 
		 xlab="log mean expression", ylab="maximum expression specificity")
	points(y= max.spec[category.matrix[,"Rb"]], x= exp.mean.log[category.matrix[,"Rb"]],cex= cex,col="red")
	points(y= max.spec[category.matrix[,"Mrp"]], x= exp.mean.log[category.matrix[,"Mrp"]],cex= cex,col="blue")
	points(y= max.spec[category.matrix[,"other_Rb"]], x= exp.mean.log[category.matrix[,"other_Rb"]],cex= 0.3,col="pink")
	
	points(y= max.spec[category.matrix[,"chrM"]], x= exp.mean.log[category.matrix[,"chrM"]],cex= cex,col="skyblue")

	points(y= max.spec[category.matrix[,"act"]], x= exp.mean.log[category.matrix[,"act"]],cex= cex,col="green")
	points(y= max.spec[category.matrix[,"hb"]], x= exp.mean.log[category.matrix[,"hb"]],cex= cex,col="purple")

	if(species=="hs")
		points(y= max.spec[category.matrix[,"MALAT1"]], x= exp.mean.log[category.matrix[,"MALAT1"]],cex=1,col="darkgreen")
	
	legend(min(exp.mean.log)+(max(exp.mean.log)-min(exp.mean.log))*0.75, 
		   max(max.spec), 
			legend=c("Ribosomal protein", "Mt Ribosomal protein", 
					 "Ribosomal pseudogenes", "chrM", "Actin","Hemoglobin","MALAT1"),
       		col=c("red", "blue", "pink", "skyblue", "green", "purple", "darkgreen"), 
       		pch=1, cex=c(rep(0.7,7)))

	if(!is.null(pdf.prefix)) dev.off()

	if(return.raw) return(cbind.data.frame(exp.mean.log, max.spec, category.matrix))
}




#' function to visualize outliers in bulk by plotting specificty score (determined from scRNA-seq referece) 
#'	vs log (mean normalized expression) determined from bulk
#' 	This function will also recommend outliers (high expression in bulk and low specificty in scRNA) for removal
#' 	Outliers are determined using ribosomal protein coding genes. 
#'	Users are recommended to add back ribosomal genes if they have been previosuly removed during preprocessing.
#' @param bulk.input: bulk RNA-seq input: a N*G matrix . rownames are bulk sample IDs, while colnames are gene IDs/names.
#' @param sc.input: scRNA-seq input: a M*G count matrix. rownames are cell IDs, while colnames are gene IDs/names.
#' 		 Or a K*G profile matrix (count scale). rownames are cell state/type names; colnames are gene IDs/names.
#' @param cell.type.labels a character vector to denote cell types of each cell (if input is a count.matrix) 
#'		or each row of GEP (if input is a "GEP")
#' @param species a character variable to denote if genes are human ("mm") or mouse ("hs")
#' @param pdf.prefix a character variable for the prefix of the pdf name. Default=NULL (plot directly).
#' @return a logical matrix of dimension #of genes * # of categories from the watchlist 
#'	to denote whether each gene belongs to a particular category
#'	an additional column low.spec.high.exp, determined based on ribosomal protein genes as described previously,
#'	will be appended to the logical matrix.  	
#' @export 
plot.bulk.outlier <- function(bulk.input,
							  sc.input,
							  cell.type.labels,
							  species,
							  pdf.prefix=NULL,
							  return.raw=FALSE){
	
	bulk.norm <- norm.to.one (ref = bulk.input, pseudo.min = 1E-8)

	#compute mean expression
	exp.mean <- colMeans(bulk.norm)
	exp.mean.log <- log(exp.mean)

	#compute gene expression specificity score
	sc.input <- collapse(ref = sc.input, labels = cell.type.labels)
	max.spec <- compute.specificity(sc.input)
		
	max.spec <- max.spec[match(colnames(bulk.norm),names(max.spec))]

	category.matrix <- assign.category(input.genes = colnames(bulk.norm), 
									   species= species)

	
	#make plots	
	if(!is.null(pdf.prefix)) pdf(paste(pdf.prefix,"_bulk_outlier.pdf",sep=""),pointsize=8,useDingbats=FALSE )
	
	cex=0.5
	plot(y= max.spec, x= exp.mean.log, cex=0.2, pch=16, col="grey", 
		 xlab="log mean expression", ylab="maximum expression specificity")
	points(y= max.spec[category.matrix[,"Rb"]], x= exp.mean.log[category.matrix[,"Rb"]],cex= cex,col="red")
	points(y= max.spec[category.matrix[,"Mrp"]], x= exp.mean.log[category.matrix[,"Mrp"]],cex= cex,col="blue")
	points(y= max.spec[category.matrix[,"other_Rb"]], x= exp.mean.log[category.matrix[,"other_Rb"]],cex= 0.3,col="pink")
	
	points(y= max.spec[category.matrix[,"chrM"]], x= exp.mean.log[category.matrix[,"chrM"]],cex= cex,col="skyblue")

	points(y= max.spec[category.matrix[,"act"]], x= exp.mean.log[category.matrix[,"act"]],cex= cex,col="green")
	points(y= max.spec[category.matrix[,"hb"]], x= exp.mean.log[category.matrix[,"hb"]],cex= cex,col="purple")

	if(species=="hs")
		points(y= max.spec[category.matrix[,"MALAT1"]], x= exp.mean.log[category.matrix[,"MALAT1"]],cex=1,col="darkgreen")

	legend(min(exp.mean.log)+(max(exp.mean.log)-min(exp.mean.log))*0.75, 
		   max(max.spec,na.rm=T), 
			legend=c("Ribosomal protein", "Mt Ribosomal protein", 
					 "Ribosomal pseudogenes", "chrM", "Actin","Hemoglobin","MALAT1"),
       		col=c("red", "blue", "pink", "skyblue", "green", "purple", "darkgreen"), 
       		pch=1, cex=c(rep(0.7,7)))
	#hist(exp.mean.log, breaks=50, col="skyblue")
	#hist(max.spec, breaks=50, col="skyblue")
	if(!is.null(pdf.prefix)) dev.off()
	
	if(return.raw) return(cbind.data.frame(exp.mean.log, max.spec, category.matrix))
}



#' function to make scatter plot while also show pearson, spearman correlation coef and MSE 
#' @param y.value: a numeric vector
#' @param x.value: a numeric vector
#' @param cex: dot size. Default=0.2
#' @param pch: dot shape. Default=16
#' @param col: dot color. Default=adjustcolor("gray",alpha.f=0.5)
#' @param title: plot title. Default=NULL.
#' @param x.lim: range of x axis. Default=NULL.
#' @param y.lim: range of y axis. Default=NULL.
#' @param ... other parameters passed to plot
plot.each <- function(y.value, 
					  x.value,
					  cex=0.2,
					  pch=16,
					  col=adjustcolor("gray",alpha.f=0.5),
					  title=NULL,
					  x.lim=NULL,
					  y.lim=NULL,
					  ...){
	
	if(is.null(x.lim)) x.lim <- range(x.value, na.rm=T)
	if(is.null(y.lim)) y.lim <- range(y.value, na.rm=T)
	
	plot(x= x.value, y= y.value, cex= cex,col= col,pch=16, 
		 ylim=y.lim,xlim=x.lim, main= title, ...)
	abline(a=0,b=1,lty=2,col="red")
		
	cor.sp <- cor(x.value, y.value, method="spearman", use="complete")
	cor.ps <- cor(x.value, y.value, use="complete" )
	mse <-  mean((x.value- y.value)^2, na.rm=T)
		
	text(x=0.75*(x.lim[2]-x.lim[1])+x.lim[1], y=0.40*(y.lim[2]-y.lim[1])+y.lim[1], labels=paste("R=", round(cor.ps,3) ,sep=""))
	text(x=0.75*(x.lim[2]-x.lim[1])+x.lim[1], y=0.30*(y.lim[2]-y.lim[1])+y.lim[1], labels=paste("rho=", round(cor.sp,3),sep=""))
	text(x=0.75*(x.lim[2]-x.lim[1])+x.lim[1], y=0.20*(y.lim[2]-y.lim[1])+y.lim[1], labels=paste("MSE=", signif(mse,3),sep=""))
	
	NULL		   	
}




#' functions that make the bulk vs scRNA scatter plot and color by gene category
#' only implemented for human TCGA 
#' @param bulk.input: bulk RNA-seq input: a N*G matrix . rownames are bulk sample IDs, while colnames are gene IDs/names.
#' @param sc.input: scRNA-seq input: a M*G count matrix. rownames are cell IDs, while colnames are gene IDs/names.
#' 		 Or a K*G profile matrix (count scale). rownames are cell state/type names; colnames are gene IDs/names.
#' @param pdf.prefix a character variable for the prefix of the pdf name. Default=NULL (plot directly).
#' @param return.value a character variable for the prefix of the pdf name. Default=NULL (plot directly).

#' @return a logical matrix of dimension #of genes * # of categories from the watchlist 
#'	to denote whether each gene belongs to a particular category
#'	an additional column low.spec.high.exp, determined based on ribosomal protein genes as described previously,
#'	will be appended to the logical matrix.  	
#' @export 
plot.bulk.vs.sc <- function(sc.input,
							bulk.input,
							pdf.prefix=NULL,
							return.value=FALSE){
								
	gene.shared <- intersect(colnames(sc.input), colnames(bulk.input))
	if(length(gene.shared)==0)
		stop("Error: gene names of reference and mixture do not match!") 
	if(length(gene.shared) < 100)
		warning("Warning: very few gene from reference and mixture match! Please double check your gene names.")

	#align reference and mixture
	bulk.input <- bulk.input[, gene.shared]
	sc.input <- sc.input[, gene.shared]
							
	bulk.norm <- norm.to.one (ref = matrix(colSums(bulk.input),nrow=1),
				 			  pseudo.min = 1E-8)[1,]
	sc.norm <- norm.to.one (ref = matrix(colSums(sc.input),nrow=1),
				 			  pseudo.min = 1E-8)[1,]
	
	gene.tab.path <- system.file("extdata","gencode.v22.broad.category.txt", package="BayesPrism")	
	gene.list <- read.table(gene.tab.path, sep="\t",header=F,stringsAsFactors=F)
	
	#decide if emsembl ID or gene symbol is used
	if( sum(substr(gene.shared,1,3)=="ENS")> length(gene.shared)*0.8  ){
		# use 80% of colnames to avoid sometimes there are manually added gene name such as GFP-XXX.
		#use EMSEMBLE ID
		#strip the "." from ENSXXX.X
		cat("EMSEMBLE IDs detected.\n")
		gene.shared <- unlist(lapply(gene.shared, function(gene.id) strsplit(gene.id,split="\\.")[[1]][1]))
		gene.df <- gene.list[match(gene.shared, gene.list[,8]),c(8,9)]
	}
	else{
		#use gene symbols
		cat("Gene symbols detected. Recommend to use EMSEMBLE IDs for more unique mapping.\n")
		gene.df <- gene.list[match(gene.shared, gene.list[,5]),c(5,9)]
	}
	colnames(gene.df) <- c("gene.name", "category")
	
	#plotting three groups: lincRNA, protein_coding and pseudogene
	
	plot.df <- data.frame(log2.bulk = log2(bulk.norm),
						  log2.sc = log2(sc.norm),
						  gene.df)
	
	selected.gene.types <- c("lincRNA", "protein_coding", "pseudogene")
	plot.df <- plot.df[plot.df$category %in% selected.gene.types,]
	
	#make plots	
	if(!is.null(pdf.prefix)) pdf(paste(pdf.prefix,"_cor_by_geneType.pdf",sep=""),pointsize=8,useDingbats=FALSE )
	
	pc.idx <- plot.df$category=="protein_coding"
	lnc.idx <- plot.df$category=="lincRNA"
	psg.idx <- plot.df$category=="pseudogene"
	
	par(mfrow = c(1,3), pty = "s")
	
	xlim <- range(plot.df[,"log2.sc"])
	ylim <- range(plot.df[,"log2.bulk"])
	plot.each(y.value = plot.df[pc.idx,"log2.bulk"], 
		 	  x.value = plot.df[pc.idx,"log2.sc"], 
		 	  xlab ="log2 total expression in scRNA-seq", 
		 	  ylab ="log2 total expression in bulk",
		 	  title = "protein coding", 
		 	  x.lim = xlim, y.lim = ylim)

	plot.each(y.value = plot.df[lnc.idx,"log2.bulk"], 
		 	  x.value = plot.df[lnc.idx,"log2.sc"], 
		 	  xlab ="log2 total expression in scRNA-seq", 
		 	  ylab ="log2 total expression in bulk",
		 	  title = "lncRNA", 
		 	  x.lim = xlim, y.lim = ylim)

	plot.each(y.value = plot.df[psg.idx,"log2.bulk"], 
		 	  x.value = plot.df[psg.idx,"log2.sc"], 
		 	  xlab ="log2 total expression in scRNA-seq", 
		 	  ylab ="log2 total expression in bulk",
		 	  title = "pseudogene", 
		 	  x.lim = xlim, y.lim = ylim)
		 	  
	if(!is.null(pdf.prefix)) dev.off()
	
	if(return.value) return(plot.df)
}




library(gplots)

#' functions that plots the heatmap for the correlation bewteen each cell type from the reference 
#' @param input: a M*G count matrix. rownames are cell IDs, while colnames are gene IDs/names.
#' 		 Or a K*G profile matrix (count scale). rownames are cell state/type names; colnames are gene IDs/names.
#' @param input.labels a character vector to denote cell stat/type of each cell (if input is a count.matrix) 
#'		or each row of GEP (if input is a "GEP")
#' @param min.exp genes expressed in at least this number of cells will be used. Default=3.
#' @param pseudo.min the desired min values to replace zero after normalization
#' @param my_palette the color palette for heatmap. Default to blue white red
#' @param title
#' @symkey symkey parameter of heatmap.2
#' @symbreaks symbreaks parameter of heatmap.2
#' @param pdf.prefix a character variable for the prefix of the pdf name. Default=NULL (plot directly).
#' @param ... other parameters for heatmap.2
#'
#' @import gplots
#' @return 
#' @export
plot.cor.phi <- function(input, 
					     input.labels,
					     min.exp=3,
					     pseudo.min=1E-8, 
					     my_palette= colorRampPalette(c("blue", "white", "red"))(n = 254), 
					     title="",
					     symkey=F,
					     symbreaks=F,
					     pdf.prefix=NULL,
					     ...){
	
	input <- input[,colSums(input) >= min.exp]
	
	ref.ct <- collapse(ref = input, labels = input.labels)
	ref.ct <- norm.to.one (ref = ref.ct, pseudo.min = pseudo.min)

	ref.ct <- scale(log2(ref.ct),center=T,scale=F)
	ref.ct <- t(ref.ct)
	cor.mat <- cor(ref.ct)
	
	hc.cols <- hclust(as.dist(1-cor.mat),method="ward.D2")
	hc.rows <- hc.cols
	
	if(!is.null(pdf.prefix)) pdf(paste(pdf.prefix,"_cor_phi.pdf",sep=""),pointsize=8,useDingbats=FALSE)
    heatmap.2(cor.mat, 
				col=my_palette, 
				density.info="none", 
				trace="none",
				dendrogram='both', 
				symm=TRUE, 
				symkey= symkey, 
				symbreaks= symbreaks,
				scale="none" , 
				Rowv=as.dendrogram(hc.rows), 
				Colv=as.dendrogram(hc.cols),
				...)

	if(!is.null(pdf.prefix)) dev.off()	
}

