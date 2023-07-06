\name{get.exp.stat}
\alias{get.exp.stat}
\title{
The function to perform differential expression test for the selection of marker genes.
}
\description{
The function to perform differential expression test using pairwise t-tests between every pair of cell states from different cell types. For each gene, it first summarizes the maximum p value and mimimum log2foldchange for each cell state (when compared with cell states from other cell types), and then take the mimimum p value and maximum log2foldchange 
for cell states within each cell type.  Used if needs to run BayesPrism over marker genes. 
}
\usage{
get.exp.stat <- function(sc.dat,
						 cell.type.labels,
						 cell.state.labels,
						 pseudo.count=0.1,
						 cell.count.cutoff=50,
						 n.cores=1)
}
\arguments{
  \item{sc.dat}{A cell-by-gene count matrix. rownames are cell state/type names; colnames are gene IDs/names.}
  \item{cell.type.labels}{A character vector to denote cell types of each cell}
  \item{cell.state.labels}{A character vector to denote cell state of each cell}
  \item{pseudo.count}{A numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.}
  \item{cell.count.cutoff}{A numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50}
  \item{n.cores}{Number of threads to be used for the pairwiseTTests function. Default=1}
}
\details{

When to use signature genes/marker genes?

Optionally, in cases where cell types are defined in a way that some of them show very similar transcription or severe batch effects exist, e.g., reference and mixture are from very different assays (ribo-depleted RNA-seq vs poly-A tail RNA-seq or PRO-seq (nascent RNA-seq) vs RNA-seq (steady state RNA)), selecting signature genes can be beneficial. This is because the selection of signature genes can enrich for genes informative for deconvolution while reducing the impact of noise caused by technical batch effects.

Users may also use signature genes derived from other differential expression methods, such as MAST.


}
\value{
A list of dataframes. Each element of the list represents a cell type. Each row of a dataframe represents a gene. Each dataframe contains two columns p.value and logFC. The output will be used as input of select.marker for selecting marker genes.
}

