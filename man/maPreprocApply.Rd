\name{maPreprocApply}
\alias{maPreprocApply}
\title{ 
    Preprocessing of microarray expression data
}
\description{
    This function prepares datasets of the GEO2KEGG microarray compendium 
    for further analysis. This includes summarization of probe level 
    expression to gene level expression as well as annotation of required 
    colData slots for sample grouping.
}
\usage{
    maPreprocApply( data.ids, parallel=NULL )
}
\arguments{
    \item{data.ids}{
        A character vector of dataset IDs.
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
}
\value{
    None, the function manipulates the datasets loaded in the workspace.  
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
\seealso{
    \code{loadEData} to load a specified expression data compendium.
}
\examples{
    # reading user-defined expression data from file
    data.ids <- loadEData("geo2kegg")

    # preprocessing the first two datasets
    maPreprocApply(data.ids[1:2])
    se <- getDataset(data.ids[1])
}   

