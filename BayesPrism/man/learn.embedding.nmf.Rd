\name{learn.embedding.nmf}
\alias{learn.embedding.nmf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Embedding learning
}
\description{
Run the embedding learning module from the output of run.prism (when malignant cells are specified by the "key" argument). 
}
\usage{
learn.embedding.nmf (bp,
                     K,
                     cycle = 50,
                     gibbs.control=list(),
                     opt.control=list(),
                     nmf.control=list(),
                     EM.res = NULL,
                     compute.elbo = F)

learn.embedding (bp,
                 eta_prior,
                 cycle = 50,
                 gibbs.control=list(),
                 opt.control=list(),
                 EM.res = NULL,
                 compute.elbo = F)								     

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bp}{A BayesPrism object.}
  \item{eta_prior}{A program-by-gene matrix (preferably on a non-log scale). Normalization will be carried out under the hood. Number of genes in eta_prior needs to match that in bp. See details.}
  \item{K}{Number of malignant gene programs (embeddings). See details.}
  \item{cycle}{Number of EM cycles. Default=50. EM typically converges in 50 cycles.}
  \item{gibbs.control}{A named list containing parameters of the Gibbs sampler. Default uses bp@control_param$gibbs.control. 
                       chain.length: length of MCMC chain. Default=1000;
                       
                       burn.in: length of burn in period. Default=500;
                       
                       thinning: retain every # of MCMC samples after the burn in period to reduce auto-correlation. Default=2;
                       
                       n.cores: number of cores to use. Default uses n.cores in the main argument;
                       
                       seed: seed number to use for repoducibility. Default = 123. Set to NULL if use pseudo-random;
                       
                       alpha: a numeric vector to represent the parameter of dirichlet distribution. }
  \item{opt.control}{A list containing parameters for the optimization step. Default uses bp@control_param$opt.control.
                     maxit: maximum number of cycles to interate. Default=100000;
                     
                     sigma: hyper-parameter of the prior if optimizer="paper". Default=2;
                     
                     other arguments passed to Rcgmin.}
  \item{nmf.control}{A list containing arguments for the nmf function. Default uses nrun=200, seed=123, and opt="vp.n.cores" (with n.cores being the min(gibbs.control$n.cores, opt.control$n.cores)).}
  \item{EM.res}{A previous EM result containing eta and omega. eta from EM.res will be used to initate additional EM runs if supplied. Default=NULL.}
  \item{compute.elbo}{A logical to denote if compute ELBO. Default=FALSE. which will slightly speed up the computation. }
}
\details{
When no prior of eta is available. Users can run learn.embedding.nmf which will initiate eta by running nmf on normalized deconvovled expression of malignant cells, i.e., t(bp.res@reference.update@psi_mal). To determine the optimum K, users may first scan a range of K using the nmf function from the package NMF over t(bp.res@reference.update@psi_mal), and then plot the metric (refer to vignette tutorial.embedding.learning.html for an example).  

In cases where users have their own priors about malignant gene programs, e.g., those derived from developmental processes or other biological processes of interest, they may run learn.embedding using a user-defined eta_prior. Ideally, eta_prior should be on a non-log scale. Normalization will be carried out under the hood. 
}
\value{
A list of eta, omega and ELBO to represent the EM results
}
\references{
Gaujoux R, Seoighe C (2010). “A flexible R package for nonnegative matrix factorization.” BMC Bioinformatics, 11(1), 367. ISSN 1471-2105, doi: 10.1186/1471-2105-11-367, https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-367.
}

