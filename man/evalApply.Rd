\name{evalApply}
\alias{evalApply}
\title{ 
    Evaluating phenotype relevance of gene set rankings 
}
\description{
    This function evaluates gene set rankings obtained from 
    the application of enrichment methods to multiple datasets - where
    each dataset investigates a certain phenotype such as a disease.
    Given pre-defined phenotype relevance scores for the gene sets, indicating
    how important a gene set is for the investigated phenotype (as e.g. judged by
    evidence from the literature), this allows to assess whether enrichment 
    methods produce gene set rankings in which phenotype-relevant gene sets 
    accumulate at the top.
}
\usage{
    evalApply( ea.ranks, rel.ranks, 
                data2pheno, perc=TRUE, top=0, rand=FALSE )
}
\arguments{
    \item{ea.ranks}{
        Enrichment analysis rankings.
        A list with an entry for each enrichment method applied.
        Each entry is a list that stores for each dataset analyzed
        the resulting gene set ranking obtained from applying the respective 
        method to the respective dataset.
        Resulting gene set rankings are assumed to be of class 
        \code{\linkS4class{DataFrame}} in which gene sets (required column 
        named \code{GENE.SET}) are ranked according to a ranking measure such as
        a gene set p-value (required column named \code{P.VALUE}). 
        See \code{EnrichmentBrowser::gs.ranking} for an example.
    }
    \item{rel.ranks}{
        Relevance score rankings.
        A list with an entry for each phenotype investigated.
        Each entry should be a \code{\linkS4class{DataFrame}} in which gene sets
        (rownames are assumed to be gene set IDs) are ranked according to a 
        phenotype relevance score (required column \code{REL.SCORE}).
    }
    \item{data2pheno}{
        A named character vector where the names correspond to dataset IDs and 
        the elements of the vector to the corresponding phenotypes investigated.
    }
    \item{perc}{
        Logical.
        Should observed scores be returned as-is or as a *perc*entage of the 
        respective optimal score. Percentages of the optimal score are typically
        easier to interpret and are comparable between datasets / phenotypes.
        Defaults to \code{TRUE}.   
    }
    \item{top}{
        Integer.
        If \code{top} is non-zero, the evaluation will be restricted to the first 
        \code{top} gene sets of each enrichment analysis ranking. 
        Defaults to \code{0}, which will then evaluate the full ranking.
    }
    \item{rand}{
        Logical.        
        Should gene set rankings be randomized to assess how likely it is to 
        observe a score equal or greater than the respective obtained score?
        Defaults to \code{FALSE}.   
    }
}
\value{
    A numeric matrix (rows = datasets, columns = methods) storing in each cell
    the relevance score sum obtained from applying the respective method to the 
    respective dataset. 
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
\seealso{
    \code{eaApply} to apply enrichment methods to multiple datasets;
    \code{readRankings} to read saved rankings as an input for the eval-functions;
    \code{evalRanks} to evaluate the phenotype relevance of one gene set ranking
        at a time;
    \code{compOpt} to compute the theoretical optimal score of a given relevance 
        score ranking; and
    \code{compRand} to compute scores for randomly drawn gene set rankings.
}
\examples{
    # simulated setup:
    # 2 methods & 2 datasets
    methods <- paste0("m", 1:2)
    data.ids <- paste0("d", 1:2)

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

    # simulate a mapping from datasets to disease codes
    d2d <- c("ALZ", "BRCA")
    names(d2d) <- data.ids

    # simulate relevance score rankings
    rel.ranks <- lapply(ea.ranks[[1]],
        function(rr)
        {
            rr[,2] <- runif(nrow(rr), min=1, max=100)
            colnames(rr)[2] <- "REL.SCORE"
            rownames(rr) <- rr[,"GENE.SET"]
            ind <- order(rr[,"REL.SCORE"], decreasing=TRUE)
            rr <- rr[ind,]
            return(rr)
        })
    names(rel.ranks) <- unname(d2d)

    # evaluate
    evalApply(ea.ranks, rel.ranks, d2d)
}   

