\name{select.gene.type}
\alias{select.gene.type}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to select genes by gene type (e.g. "protein_coding", "pseudogene", "lincRNA" ).
}
\description{
A function to select genes by gene type, as informed by the plot from plot.bulk.vs.sc. Only works for human data, as the annotation was for TCGA (GENCODE v22). For other species please filter manually if needed.

Applying to either bulk or scRNA-seq is sufficient. As the contruction of prism object will automatically take the intersection of genes from mixture and reference.
}
\usage{
select.gene.type (input,
                  gene.type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{Either mixture or scRNA-seq reference count matrix. If input is a reference, colnames are gene IDs/names, and rownames are cell IDs or cell state/type names if input is a GEP. If the input is bulk , colnames are gene IDs/names, and rownames are sample IDs.}
  \item{gene.type}{A character vector to input gene groups to be retained, must be one or more elements from	c("protein_coding", "pseudogene", "lincRNA")}
}
\value{
Return a matrix of the same nrow as input but only keep ncol of number of genes selected.
}

