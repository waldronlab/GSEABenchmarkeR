\name{readResults}
\alias{readRankings}
\alias{readRuntime}
\title{ 
    Reading results of enrichment analysis
}
\description{
    These functions read results obtained from the application of enrichment 
    methods to multiple datasets for subsequent assessment.
}
\usage{
    readRuntime( data.dir, data.ids, methods )
    
    readRankings( data.dir, data.ids, methods )
}
\arguments{
    \item{data.dir}{
        Character.
        The data directory where results have been saved to. 
    }
    \item{data.ids}{
        A character vector of dataset IDs.
    }
    \item{methods}{
        Methods for enrichment analysis.
        A character vector with method names typically chosen from 
        \code{EnrichmentBrowser::sbea.methods} and 
        \code{EnrichmentBrowser::nbea.methods}, or user-defined 
        functions implementing methods for enrichment analysis.
    }
}
\value{
    A result list with an entry for each method applied. 
    Each entry stores corresponding
    runtimes (for \code{readRuntimes}) and gene set rankings 
    (for \code{readRankings}) as obtained from applying the 
    respective method to the given datasets.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
\seealso{
    \code{eaApply} to apply enrichment methods to multiple datasets.
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
            r <- EnrichmentBrowser::make.example.data("ea.res") 
            r <- EnrichmentBrowser::gs.ranking(r, signif.only=FALSE)
            saveRDS(r, file=out.file)   
        }
    }

    # reading runtime & rankings
    rts <- readRuntime(res.dir, data.ids, methods)
    rkgs <- readRankings(res.dir, data.ids, methods)
}   

