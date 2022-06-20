\name{plot.bulk.vs.sc}
\alias{plot.bulk.vs.sc}
\title{
Make the bulk vs scRNA scatter plot by gene type (coding, lncRNA and psuedogenes)
}
\description{

This function check the concordance of gene expression for different types of genes. As bulk and single cell data are usually collected by different experimental protocols, they may have different sensitivity to different types of genes. Selecting gene types of higher concordance between two assays may help reduce batch effects and improve the accruacy of deconvolution.

Note that this function only works for human data (using GENCODE v22 to match TCGA). For other species, you are advised to make plots by yourself.

}
\usage{
plot.bulk.vs.sc (sc.input,
                 bulk.input,
                 pdf.prefix=NULL,
                 return.value=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sc.input}{A cell-by-gene count matrix. rownames are cell IDs, while colnames are gene IDs/names. Or a cell type-by-gene profile matrix (non-log scale). rownames are cell state/type names; colnames are gene IDs/names.}
  \item{bulk.input}{bulk RNA-seq input: a sample-by-gene count matrix. rownames are bulk sample IDs, while colnames are gene IDs/names.}
  \item{pdf.prefix}{A character for the prefix of the name of pdf outputs. Default=NULL (plot directly).}
  \item{return.value}{Logical to denote if return the data matrix for plotting. Default=FALSE.} 
}
\details{
Three genes types are plotted. R: Pearson's correlation; Rho=spearman's correlation; MSE: mean squared error. 
}
\value{
Returns a dataframe. 
Each row is a gene (shared between bulk and scRNA). 
Colnames are log2.bulk: log2 normalized bulk expression; log2.sc: log2 normalized scRNA-seq expression.
}
