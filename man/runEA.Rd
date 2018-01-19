\name{runEA}
\alias{runEA}
\title{ 
    Execution of enrichment analysis
}
\description{
    This is a wrapper function to carry out specific enrichment 
    methods and save results for subsequent benchmarking. 
}
\usage{
    runEA( se, method, gs, save2file=FALSE, out.dir=NULL, ... )
}
\arguments{
    \item{se}{
        A \code{\linkS4class{SummarizedExperiment}} storing the expression
        dataset. See the documentation of \code{EnrichmentBrowser::sbea} for
        required minimal annotations. 
    }
    \item{method}{
        The method for enrichment analysis.
        This should be either a character vector of length 1, when using
        a method from \code{EnrichmentBrowser::sbea.methods} or 
        \code{EnrichmentBrowser::nbea.methods}, or a user-defined 
        function implementing a method for enrichment analysis.
    }
    \item{gs}{
        Gene sets, i.e. a list of character vectors of gene IDs.
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
        Additional arguments passed to the selected enrichment method.
    }
}
\value{
    A list with 2 elements: (1) the runtime of executing the chosen method on
    the given dataset, and (2) a \code{\linkS4class{DataFrame}} storing the 
    resulting gene set ranking.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
\seealso{
    \code{EnrichmentBrowser::sbea} and \code{EnrichmentBrowser::nbea}
}
\examples{
    # (1) expression data: 
    # simulated expression values of 100 genes
    # in two sample groups of 6 samples each
    eset <- EnrichmentBrowser::make.example.data(what="eset")
    eset <- EnrichmentBrowser::de.ana(eset) 

    # (2) gene sets:
    # draw 10 gene sets with 15-25 genes
    gs <- EnrichmentBrowser::make.example.data(what="gs", gnames=names(eset))

    # (4) performing the enrichment analysis
    ea.res <- runEA(eset, method="ora", gs=gs)
}   

