\name{evalNrSigSets}
\alias{evalNrSigSets}
\alias{evalNrSets}
\title{ 
    Evaluating gene set rankings for the number of (significant) sets
}
\description{
    These functions evaluate gene set rankings obtained from 
    applying enrichment methods to multiple datasets.
    This allows to assess resulting rankings for granularity 
    (how many gene sets have a unique p-value?) and statistical significance 
    (how many gene sets have a p-value below a significance threshold?). 
}
\usage{
    evalNrSets( ea.ranks, uniq.pval=TRUE, perc=TRUE )
    
    evalNrSigSets( ea.ranks, alpha=0.05, padj="none", perc=TRUE )
}
\arguments{
    \item{ea.ranks}{
        Enrichment analysis rankings.
        A list with an entry for each enrichment method applied.
        Each entry is a list that stores for each dataset analyzed
        the resulting gene set ranking as obtained from applying the 
        respective method to the respective dataset.
    }
    \item{uniq.pval}{
        Logical.
        Should the number of gene sets with a unique p-value or the total 
        number of gene sets per ranking be returned? 
        Defaults to \code{TRUE}.
    }
    \item{perc}{
        Logical.
        Should the percentage or absolute number of gene sets be returned?
        Percentage is typically more useful for comparison between rankings
        with a potentially different total number of gene sets. 
        Defaults to \code{TRUE}.
    }
    \item{alpha}{
        Statistical significance level. Defaults to 0.05.
    }
    \item{padj}{
        Character.
        Method for adjusting p-values to multiple testing.  For
        available methods see the man page of the stats
        function \code{p.adjust}.  Defaults to \code{BH}.
    }
}
\value{
    A list of numeric vectors storing for each method the number of 
    (significant) gene sets for each dataset analyzed.
    If each element of the resulting list is of equal length 
    (corresponds to successful application of each enrichment method to each
    dataset), the list is automatically simplified to a numeric matrix 
    (rows = datasets, columns = methods).
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
\seealso{
    \code{eaApply} to apply enrichment methods to multiple datasets;
    \code{readRankings} to read saved rankings as an input for the eval-functions.
}
\examples{
    # simulated setup:
    # 3 methods & 5 datasets
    methods <- paste0("m", 1:3)
    data.ids <- paste0("d", 1:5)

    # simulate gene set rankings
    ea.ranks <- sapply(methods, function(m) 
            sapply(data.ids, 
                function(d)
                {
                    r <- EnrichmentBrowser::make.example.data("ea.res") 
                    r <- EnrichmentBrowser::gs.ranking(r, signif.only=FALSE)
                    return(r)
                }, simplify=FALSE),
                simplify=FALSE)

    # evaluate
    evalNrSets(ea.ranks)
    evalNrSigSets(ea.ranks)
}   

