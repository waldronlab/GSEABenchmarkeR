% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ea.R
\name{runEA}
\alias{runEA}
\title{Application of enrichment methods to multiple datasets}
\usage{
runEA(
  exp.list,
  methods,
  gs,
  perm = 1000,
  parallel = NULL,
  save2file = FALSE,
  out.dir = NULL,
  ...
)
}
\arguments{
\item{exp.list}{Experiment list.  A \code{list} of datasets, each being of
class \code{\linkS4class{SummarizedExperiment}}. In case of just one dataset
a single \code{\linkS4class{SummarizedExperiment}} is also allowed. 
See the documentation of \code{\link{sbea}} for required minimal annotations.}

\item{methods}{Methods for enrichment analysis.  This can be either \itemize{ 
\item a character vector with method names chosen from \code{\link{sbeaMethods}}
and \code{\link{nbeaMethods}}, 
\item a user-defined function implementing a method for enrichment analysis, or
\item a named list, containing pre-defined and/or user-defined enrichment methods.
See examples.}}

\item{gs}{Gene sets, i.e. a list of character vectors of gene IDs.}

\item{perm}{Number of permutations of the sample group assignments. 
Defaults to 1000. Can also be an integer vector matching
the length of \code{methods} to assign different numbers of permutations for
different methods.}

\item{parallel}{Parallel computation mode.  An instance of class
\code{\linkS4class{BiocParallelParam}}.  See the vignette of the
\code{BiocParallel} package for switching between serial, multi-core, and
grid execution.  Defaults to \code{NULL}, which then uses the first element
of \code{BiocParallel::registered()} for execution.  If not changed by the
user, this accordingly defaults to multi-core execution on the local host.}

\item{save2file}{Logical. Should results be saved to file for subsequent
benchmarking?  Defaults to \code{FALSE}.}

\item{out.dir}{Character.  Determines the output directory where results are
saved to.  Defaults to \code{NULL}, which then writes to
\code{tools::R_user_dir("GSEABenchmarkeR")} in case \code{save2file}
is set to \code{TRUE}.}

\item{...}{Additional arguments passed to the selected enrichment methods.}
}
\value{
A list with an entry for each method applied.  Each method entry is
a list with an entry for each dataset analyzed.  Each dataset entry is a
list of length 2, with the first element being the runtime and the second
element being the gene set ranking, as obtained from applying the respective
method to the respective dataset.
}
\description{
This function applies selected methods for enrichment analysis to selected
datasets of a compendium.
}
\examples{

    # loading three datasets from the GEO2KEGG compendium
    geo2kegg <- loadEData("geo2kegg", nr.datasets=3)

    # only considering the first 1000 probes for demonstration
    geo2kegg <- lapply(geo2kegg, function(d) d[1:1000,]) 

    # preprocessing and DE analysis for two of the datasets
    geo2kegg <- maPreproc(geo2kegg[2:3])
    geo2kegg <- runDE(geo2kegg)

    # getting a subset of human KEGG gene sets
    gs.file <- system.file("extdata/hsa_kegg_gs.gmt", package="EnrichmentBrowser")
    kegg.gs <- EnrichmentBrowser::getGenesets(gs.file)

    # applying two methods to two datasets 
    res <- runEA(geo2kegg, methods=c("ora", "camera"), gs=kegg.gs, perm=0)

    # applying a user-defined enrichment method
    dummySBEA <- function(se, gs) 
    {
         sig.ps <- sample(seq(0, 0.05, length=1000), 5)
         nsig.ps <- sample(seq(0.1, 1, length=1000), length(gs)-5)
         ps <- sample(c(sig.ps, nsig.ps), length(gs))
         names(ps) <- names(gs)
         return(ps)
    }  
    res <- runEA(geo2kegg, methods=dummySBEA, gs=kegg.gs)

    # applying a mix of pre-defined and user-defined methods
    methods <- list(camera = "camera", dummySBEA = dummySBEA)
    res <- runEA(geo2kegg, methods, gs=kegg.gs, perm=0)   


}
\seealso{
\code{\link{sbea}} and \code{\link{nbea}}
for carrying out set- and network-based enrichment analysis.

\code{\linkS4class{BiocParallelParam}} and \code{\link{register}} for
configuration of parallel computation.
}
\author{
Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
