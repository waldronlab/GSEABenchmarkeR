\name{eaApply}
\alias{eaApply}
\title{ 
    Application of enrichment methods to multiple datasets
}
\description{
    This function applies selected methods for enrichment analysis
    to selected datasets of a compendium.
}
\usage{
    eaApply( exp.list, methods, gs, 
                parallel=NULL, save2file=FALSE, out.dir=NULL, ... )
}
\arguments{
    \item{exp.list}{
        Experiment list.
        A \code{list} of datasets, each being of class 
        \code{\linkS4class{ExpressionSet}}.
    }
    \item{methods}{
        Methods for enrichment analysis.
        A character vector with method names chosen from 
        \code{EnrichmentBrowser::sbea.methods} and 
        \code{EnrichmentBrowser::nbea.methods}, or user-defined 
        functions implementing methods for enrichment analysis.
    }
    \item{gs}{
        Gene sets, i.e. a list of character vectors of gene IDs.
    }
    \item{parallel}{
        Parallel computation mode.
        An instance of class \code{\linkS4class{BiocParallelParam}}.
        See the vignette of the \code{BiocParallel} package for switching 
        between serial, multi-core, and grid execution.
        Defaults to \code{NULL}, which then uses the first element of 
        \code{BiocParallel::registered()} for execution.
        If not changed by the user, this accordingly defaults to multi-core
        execution on the local host. 
    }
    \item{save2file}{
        Logical. Should results be saved to file for subsequent benchmarking?
        Defaults to \code{FALSE}.
    } 
    \item{out.dir}{
        Character.
        Determines the output directory where results are saved to. 
        Defaults to \code{NULL}, which then writes to 
        \code{rappdirs::user_data_dir("GSEABenchmarkeR")} in case
        \code{save2file} is set to \code{TRUE}.
    }
    \item{...}{
        Additional arguments passed to the selected enrichment methods.
    }
}
\value{
    A list with an entry for each method applied.
    Each method entry is a list with an entry for each dataset analyzed.
    Each dataset entry is a list of length 2, with the first element being
    the runtime and the second element being the gene set ranking, as obtained 
    from applying the respective method to the respective dataset.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
\seealso{
    \code{EnrichmentBrowser::sbea} and \code{EnrichmentBrowser::nbea}
    for carrying out set- and network-based enrichment analysis.

    \code{\linkS4class{BiocParallelParam}} and \code{BiocParallel::register}
    for configuration of parallel computation.
}
\examples{
    # loading two datasets from the GEO2KEGG compendium
    geo2kegg <- loadEData("geo2kegg", nr.datasets=2)

    # preprocessing and DE analysis
    geo2kegg <- maPreprocApply(geo2kegg)
    geo2kegg <- deApply(geo2kegg)

    # getting a subset of human KEGG gene sets
    gs.file <- system.file("extdata/hsa_kegg_gs.gmt", package="EnrichmentBrowser")
    kegg.gs <- EnrichmentBrowser::parse.genesets.from.GMT(gs.file)

    # applying two methods to two datasets 
    res <- eaApply(geo2kegg, methods=c("ora", "camera"), gs=kegg.gs)
    
}   

