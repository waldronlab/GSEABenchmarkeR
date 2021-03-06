% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{readResults}
\alias{readResults}
\title{Reading results of enrichment analysis}
\usage{
readResults(
  data.dir,
  data.ids,
  methods,
  type = c("runtime", "ranking", "typeI")
)
}
\arguments{
\item{data.dir}{Character.  The data directory where results have been saved
to.}

\item{data.ids}{A character vector of dataset IDs.}

\item{methods}{Methods for enrichment analysis.  A character vector with
method names typically chosen from \code{\link{sbeaMethods}}
and \code{\link{nbeaMethods}}, or user-defined functions
implementing methods for enrichment analysis.}

\item{type}{Character. Type of the result. Should be one out of 'runtime',
'ranking', or 'typeI'.}
}
\value{
A result list with an entry for each method applied.  Each entry
stores corresponding runtimes (\code{type="runtime"} ), gene set
rankings (\code{type="ranking"}), or type I error rates (\code{type="typeI"})
as obtained from applying the respective method to the given datasets.
}
\description{
These functions read results obtained from the application of
enrichment methods to multiple datasets for subsequent assessment.
}
\examples{

    # simulated setup: 
    # 1 methods & 1 datasets
    methods <- paste0("m", 1:2)
    data.ids <- paste0("d", 1:2)

    # result directory
    res.dir <- tempdir()
    sdirs <- file.path(res.dir, methods)
    for(d in sdirs) dir.create(d)
    
    # store runtime & rankings 
    for(m in 1:2)
    {
        rt <- runif(5, min=m, max=m+1)
        for(d in 1:2)
        {
            # runtime
            out.file <- paste(data.ids[d], "txt", sep=".")
            out.file <- file.path(sdirs[m], out.file)
            cat(rt[d], file=out.file) 

            # ranking
            out.file <- sub("txt$", "rds", out.file)
            r <- EnrichmentBrowser::makeExampleData("ea.res") 
            r <- EnrichmentBrowser::gsRanking(r, signif.only=FALSE)
            saveRDS(r, file=out.file)   
        }
    }

    # reading runtime & rankings
    rts <- readResults(res.dir, data.ids, methods, type="runtime")
    rkgs <- readResults(res.dir, data.ids, methods, type="ranking")

}
\seealso{
\code{runEA} to apply enrichment methods to multiple datasets.
}
\author{
Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
