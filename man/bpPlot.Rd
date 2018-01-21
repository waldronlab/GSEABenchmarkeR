\name{bpPlot}
\alias{bpPlot}
\title{ 
    Customized boxplot visualization of benchmark results
}
\description{
    This is a convenience function to create customized boxplots
    for specific benchmark criteria such as runtime, statistical
    significance and phenotype relevance.
}
\usage{
    bpPlot( data, what=c("runtime", "sig.sets", "rel.sets") )
}
\arguments{
    \item{data}{
        Numeric matrix or list of numeric vectors. 
        In case of a matrix, column names are assumed to be method names
        and rownames are assumed to be dataset IDs.
        In case of a list, names are assumed to be methods names and each
        element corresponds to a numeric vector with names assumed to be
        dataset IDs.
    }
    \item{what}{
        Character.
        Determines how the plot is customized.
        One of
            \itemize{
                \item runtime: displays runtime of methods across datasets,
                \item sig.sets: displays percentage of significant gene sets,
                \item rel.sets: displays phenotype relevance scores.
            }
    }
}
\value{
    None. Plots to a graphics device.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
\seealso{
    \code{evalNrSigSets} to evaluate fractions of significant gene sets, 
    
    \code{evalRanks} to evaluate phenotype relevance of gene set rankings.
}
\examples{
    # simulated setup:
    # 3 methods & 5 datasets
    methods <- paste0("m", 1:3)
    data.ids <- paste0("d", 1:5)

    # runtime data
    rt <- vapply(1:3, function(m) runif(5, min=m, max=m+1), numeric(5))
    rownames(rt) <- data.ids 
    colnames(rt) <- methods

    # plot
    bpPlot(rt, what="runtime")
}   

