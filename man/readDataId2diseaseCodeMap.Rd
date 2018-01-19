\name{readDataId2diseaseCodeMap}
\alias{readDataId2diseaseCodeMap}
\title{ 
    Read a mapping between dataset ID and disease code
}
\description{
    When assessing enrichment analysis results for phenotype relevance,
    it is assumed that each analyzed dataset investigates a certain phenotype 
    such as a disease. 
    This function reads a mapping between dataset IDs and assigned disease codes.
}
\usage{
    readDataId2diseaseCodeMap( map.file )
}
\arguments{
    \item{map.file}{
        Character. The path to the mapping file. 
    }
}
\value{
    A named character vector where each element of the vector is a disease code 
    and the names are the dataset IDs.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
\seealso{
    \code{evalApply} for evaluating phenotype relevance of gene set rankings.
}
\examples{
    data.dir <- system.file("extdata", package="GSEABenchmarkeR")
    d2d.file <- file.path(data.dir, "malacards", "GseId2Disease.txt")
    d2d.map <- readDataId2diseaseCodeMap(d2d.file)
}   

