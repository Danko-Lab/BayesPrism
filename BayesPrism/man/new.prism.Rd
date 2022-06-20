\name{new.prism}
\alias{new.prism}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct a prism object
}
\description{
Function to contruct a prism object from user provided data.
}
\usage{
new.prism (reference,
           input.type,
           cell.type.labels,
           cell.state.labels,
           key,
           mixture,
           outlier.cut=0.01,
					 outlier.fraction=0.1,
           pseudo.min=1E-8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reference}{input of scRNA-seq count matrix (cell-by-gene) or a GEP matrix (cell state/type-by-gene).	rownames(reference) are either cell IDs (if input.type="count.matrix") or names of cell states/types (if input.type="GEP"). colnames(reference) are genes identifiers (e.g. gene symbols or EMSEMBL IDs). See details. }
  \item{input.type}{A character string. ="count.matrix" if using the raw count matrix of scRNA-seq as input; ="GEP" if using a cell type/state-by-gene gene expression profile (GEP) matrix. See details.}
  \item{cell.type.labels}{a character vector to denote cell types of each cell (if input.type="count.matrix") or each row of GEP (if input.type="GEP"). See details.}
  \item{cell.state.labels}{a character vector to denote cell state of each cell (if input.type="count.matrix") or each row of GEP (if input.type="GEP"). See details.}
  \item{key}{The character string that correponds to malignant cells in cell.type.labels. Set to NULL if there are no malignant cells or the malignant cells between reference and mixture are from matched sample, in which case all cell types will be treated equally.} 
  \item{mixture}{The bulk RNA-seq matrix (sample-by-gene). colnames(mixture) need to match the annotation of colnames(reference). new.prism will intersect colnames(reference) and colnames(mixture) and perform deconvolution using the intersected genes.}
    \item{outlier.cut & outlier.fraction}{Filter genes in X whose expression fraction is greater than outlier.cut (Default=0.01) in more than outlier.fraction (Default=0.1) of bulk data. Typically for dataset with reasonable quality control, very few genes will be filtered. Removal of outlier genes will ensure that the inference will not be dominated by outliers, which sometimes may be resulted from poor QC in mapping.}  
  \item{pseudo.min}{The desired min value to replace zero after normalization. Default=1E-8. From v2.0 pseudo.min=0 is also supported (may need to increase alpha to prevent underflow). }
}
\details{

Can I use normalized or transformed reference / mixture?

We recommend the use of unnormalized and untransformed count data. When raw count is not available, linear normalization, such as TPM, RPM, RPKM, FPKM, is also acceptable, as BayesPrism is robust to linear multiplicative difference between the reference and mixture. Ideally, if using normalized data, it is the best to supply reference and bulk transformed by the same method. log transformation should be avoided.


Should I use count.matrix or GEP? 

When using scRNA-seq count matrix as the input (recommended), user needs to specify input.type = "count.matrix". The other option for input.type is "GEP" (gene expression profile) which is a cell state-by-gene matrix. The option of GEP is used when using reference derived from other assays, such as sorted bulk data. When using GEP, summed raw count (un-normalized) for each cell state is preferred. Although normalized (non-log scale) is also acceptable. new.prism will internally take the mean of each gene across cell states within each cell type to generate a cell type-by-gene GEP, i.e., \eqn{\phi} of cell types. 


Why use and what to supply for What to supply for cell.state.labels and cell.type.labels?

The definition of cell type and cell state can be somewhat arbitrary (similar to the issue of assigning cell types for scRNA-seq) and depends on the question of interest. Their definitions depend on the granularity we aim at and the confidence of the cell.type.labels in scRNA-seq data. Usually, a good rule of thumb is as follows. 1) Define cell types as the cluster of cells having a sufficient number of significantly differentially expressed genes than other cell types, e.g., greater than 50 or even 100. For clusters that are too similar in transcription, we recommend treating them as cell states, which will be summed up before the final Gibbs sampling. Therefore, cell states are often suitable for cells that form a continuum on the phenotypic manifold rather than distinct clusters. 2) Define multiple cell states for cell types of significant heterogeneity, such as malignant cells, and of interest to deconvolve their transcription.

}
\value{
Return a prism object of the following slots:

phi_cellState: an S4 object of the class refPhi that stores the reference profile of cell states.

phi_cellType: an S4 object of the class refPhi that stores the reference profile of cell types.

map: a named list that stores correspondence between cell type (name of each element) and cell states (content of each element). 

key: The character string that correponds to malignant cells in cell.type.labels. =_NA_value if is set to NULL.

mixture: The raw count matrix of bulk.
}

