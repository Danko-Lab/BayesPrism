\name{plot.scRNA.outlier}
\alias{plot.scRNA.outlier}
\title{
Visualize and determine outliers in scRNA-seq reference
}
\description{
Gene expressed at high magnitude, such as ribosomal protein genes and mitochondrial genes, may dominate the distribution and bias the inference. These genes are often not informative in distinguishing cell types and can be a source of large spurious variance. As a result, they can be detrimental to deconvolution. We recommend the removal of these genes.

This function plots the cell type specificity score vs. log (mean normalized expression) for each gene in scRNA-seq reference.

Users may apply additional filtering using the statistics returned this function (by specifying return.raw=TRUE) if necessary. 

}
\usage{
plot.scRNA.outlier (input,
                    cell.type.labels,
                    species,
                    pdf.prefix=NULL,
                    return.raw=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{A cell-by-gene count matrix. rownames are cell IDs, while colnames are gene IDs/names. Or a cell type-by-gene profile matrix (non-log scale). rownames are cell state/type names; colnames are gene IDs/names.}
  \item{cell.type.labels}{A character vector to denote cell types of each cell (if input is a count.matrix) or each row of GEP (if input is a GEP)}
  \item{species}{A character variable to denote if genes are human ("mm") or mouse ("hs"). Currently only these two species are supported.}
  \item{pdf.prefix}{A character for the prefix of the name of pdf outputs. Default=NULL (plot directly).}
  \item{return.raw}{A logical to denote if return the data frame of computed statistics. Default=FALSE. .}
}
\details{
Several categories of genes are visualized on the plot. Users can inspect their distribution and determine the cutoffs to use.
}
\value{
	
A data frame containing	log of mean expression across cell types (exp.mean.log) and the maximum specificity score of each gene, followed by a logical matrix of dimension gene-by-category to represent if each gene belongs to the putative outliers of each category. These categories are:

"other_Rb": curated ribosome-related genes. They contain mostly ribosomal psuedogenes. We suggest to exclude these genes.

"chrM": genes from mitochondrial chromosome. They usually represent techincal artifacts, such as cell quality. We suggest to exclude these genes.

"chrX and chrY": genes from sex chromosomes. If reference and mixture are from different sexes, we suggest to exclude these gene.

"Mrp": mitochondrial ribosomal protein genes. We suggest to exclude these gene.

"act": actin genes. We suggest to exclude the outliers in these genes (as defined by the low.spec.high.exp column).

"hb": actin genes. We suggest to exclude the outliers in these genes (as defined by the low.spec.high.exp column).

"MALAT1": The MALAT1 gene. This genes are sometimes transcribed to an extremely high value in tumor samples. We suggest to exclude it if it is marked as an outlier (as defined by the low.spec.high.exp column).

}
