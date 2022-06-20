\name{collapse}
\alias{collapse}
\title{
Summing up reads across cells from each cell type/state
}
\description{
A utility function used to sum up reads across cells from each cell type/state. The output is a gene profile matrix (GEP) of the dimension cell type/state-by-gene. Usually used in conjunction with get.high.spec.genes to get signature genes. 
}
\usage{
collapse (ref, 
          labels)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ref}{A cell-by-gene count matrix. rownames are cell IDs, while colnames are gene IDs/names.}
  \item{labels}{A character vector of length equal to nrow(ref) to denote the cell type/states of each cell.}
}
\value{
A gene profile matrix (GEP) of the dimension cell type/state-by-gene.
}
\seealso{
get.high.spec.genes
}
