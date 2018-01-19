\name{getDataset}
\alias{getDataset}
\title{ 
    Obtaining an expression dataset by name
}
\description{
    This is a convenience function to obtain specified
    expression datasets from either one of the pre-defined 
    compendia (GEO2KEGG or TCGA) or a user-defined compendium.
    It requires that the respective compendium has been loaded 
    into the user's workspace.
}
\usage{
    getDataset( data.id )
}
\arguments{
    \item{data.id}{
        Dataset ID. 
        A character vector of length 1.
    }
}
\value{
    A \code{\linkS4class{SummarizedExperiment}}.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
\seealso{
    \code{loadEData} to load a specified expression data compendium.
}
\examples{
    # reading user-defined expression data from file
    data.dir <- system.file("extdata/myEData", package="GSEABenchmarkeR")
    data.ids <- loadEData(data.dir)

    # obtaining the first dataset
    se <- getDataset(data.ids[1])
}   

