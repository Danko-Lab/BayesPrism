\name{select.marker}
\alias{select.marker}
\title{
The function to perform differential expression test for the selection of marker genes.
}
\description{
The function to perform differential expression test using pairwise t-tests between every pair of cell states from different cell types. For each gene, it first summarizes the maximum p value and mimimum log2foldchange for each cell state (when compared with cell states from other cell types), and then take the mimimum p value and maximum log2foldchange 
for cell states within each cell type.  Used if needs to run BayesPrism over marker genes. 
}
\usage{
select.marker (sc.dat,
						  stat,
						  pval.max=0.01,
						  lfc.min=0.1)
}
\arguments{
  \item{sc.dat}{A cell-by-gene count matrix. rownames are cell state/type names; colnames are gene IDs/names.}
  \item{stat}{A list of dataframes with p.value and logFC, outputted from get.exp.stat.}
  \item{pval.max}{maximum p value. Default=0.01.}
  \item{lfc.min}{lfc.min mimimum log2 fold change. Dafult=0.1}
  \item{cell.count.cutoff}{A numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50}
}
\details{

Users are adviced to inspect the number of differentially expressed genes (DEG) in each cell type (e.g. at the default cutoff). Having sufficient number of DEG, e.g ~50, for each cell type suggests good confidence of posterior cell type fraction.

}
\value{
Return a count matrix with ncol of marker genes.
}

