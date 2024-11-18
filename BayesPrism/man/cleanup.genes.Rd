\name{cleanup.genes}
\alias{cleanup.genes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to filter genes by gene category and to exlcude lowly expressed genes.
}
\description{
A function to filter genes by gene category (Rb, Mrp, chrM, etc.) and to exlcude lowly expressed genes.
}
\usage{
cleanup.genes (input,
						   input.type, 
						   species, 
						   gene.group,
						   exp.cells=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{Either mixture or scRNA-seq reference count matrix. If input is a reference, colnames are gene IDs/names, and rownames are cell IDs or cell state/type names if input is a GEP. If the input is bulk , colnames are gene IDs/names, and rownames are sample IDs.}
  \item{input.type}{A character string. ="count.matrix" if using the raw count matrix of scRNA-seq as input; 	="GEP" if using a gene expression profile, preferably summed raw count(un-normalized), as input.}
  \item{species}{A character variable to denote if genes are human ("hs") or mouse ("mm").}
  \item{gene.group}{a character vector to input gene groups to be removed, must be one or more elements from c("other_Rb","chrM","chrX","chrY","Rb","Mrp","act","hb","MALAT1")}
  \item{exp.cells}{Genes expressed in number of cells fewer than this will be excluded. Default=1. If the input is GEP, gene will be selected by automatically setting exp.cells is set to min(exp.cells,1). As a result genes expressed in at least 0 or 1 cell type will be retained.}  
}

\value{
Return a matrix of the same nrow as input but only keep ncol of number of genes retained.
}

