\name{evalRanks}
\alias{evalRanks}
\alias{compOpt}
\alias{compRand}
\title{ 
    Evaluating phenotype relevance of gene set rankings 
}
\description{
    These functions evaluate a gene set ranking obtained from 
    the application of an enrichment method to an expression dataset - that 
    investigates a certain phenotype such as a disease.
    Given pre-defined phenotype relevance scores for the gene sets, indicating
    how important a gene set is for the investigated phenotype (as e.g. judged by
    evidence from the literature), this allows to assess whether the ranking 
    accumulates phenotype-relevant gene sets at the top.
}
\usage{
    evalRanks( ea.ranks, rel.ranks, top=0 )

    compOpt( rel.ranks, gs, top=0 )

    compRand( rel.ranks, gs, n=1000 )    
}
\arguments{
    \item{ea.ranks}{
        Enrichment analysis ranking.
        A \code{\linkS4class{DataFrame}} in which gene sets (required column 
        named \code{GENE.SET}) are ranked according to a ranking measure such as
        a gene set p-value (required column named \code{P.VALUE}). 
        See \code{EnrichmentBrowser::gs.ranking} for an example.
    }
    \item{rel.ranks}{
        Relevance score ranking.
        A \code{\linkS4class{DataFrame}} in which gene sets (rownames are assumed
        to be gene set IDs) are ranked according to a phenotype relevance score
        (required column \code{REL.SCORE}).
    }
    \item{top}{
        Integer.
        If \code{top} is non-zero, the evaluation will be restricted to the first 
        \code{top} gene sets of the enrichment analysis ranking. 
        Defaults to \code{0}, which will then evaluate the full ranking.
    }
    \item{gs}{
        Character vector of gene set IDs on which enrichment analysis has
        been carried out. 
    }
    \item{n}{
        Integer. Number of permutations.
    }
}
\value{
    For \code{evalRanks} and \code{compOpt}, a numeric value corresponding to the
    observed and optimal score, respectively.
    For \code{compRand}, a numeric vector storing the \code{n} random scores. 
}
\details{
    The function \code{evalRanks} evaluates the similarity of a gene set ranking 
    obtained from enrichment analysis and a gene set ranking based on phenotype
    relevance scores. Therefore, the function first transforms the ranks 'r' from 
    the enrichment analysis to weights 'w' in [0,1] via w = 1 - r / N; where 'N'
    denotes the total number of gene sets on which the enrichment analysis has 
    been carried out. These weights are then multiplied with the corresponding 
    relevance scores and summed up.

    The function \code{compOpt} applies \code{evalRanks} to the 
    theoretically optimal case in which the enrichment analysis ranking is 
    identical to the relevance score ranking. The ratio between observed and
    optimal score is useful for comparing observed scores between datasets / 
    phenotypes.

    The function \code{compRand} repeatedly applies \code{evalRanks} to randomly
    drawn gene set rankings to assess how likely it is to observe a score equal 
    or greater than the one obtained. 
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
\seealso{
    \code{runEA} to execute a gene set enrichment analysis.
}
\examples{
    # simulate gene set ranking
    ea.ranks <- EnrichmentBrowser::make.example.data("ea.res")
    ea.ranks <- EnrichmentBrowser::gs.ranking(ea.ranks, signif.only=FALSE)

    # simulated relevance score ranking
    rel.ranks <- ea.ranks
    rel.ranks[,2] <- runif(nrow(ea.ranks), min=1, max=100)
    colnames(rel.ranks)[2] <- "REL.SCORE"
    rownames(rel.ranks) <- rel.ranks[,"GENE.SET"]
    ind <- order(rel.ranks[,"REL.SCORE"], decreasing=TRUE)
    rel.ranks <- rel.ranks[ind,]

    # evaluate
    evalRanks(ea.ranks, rel.ranks)    
    compOpt(rel.ranks, ea.ranks[,"GENE.SET"])
    compRand(rel.ranks, ea.ranks[,"GENE.SET"], n=3)
}   

