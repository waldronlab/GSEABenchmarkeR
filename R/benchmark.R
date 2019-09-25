############################################################
# 
# author: Ludwig Geistlinger
# date: 2016-07-29 16:49:12
# 
# descr: core benchmark functionality 
# 
############################################################

#
# iterate over expression datasets
#
.iter <- function(exp.list, f, ..., parallel=NULL)
{
    is.windows <- is.null(parallel) && .Platform$OS.type == "windows"
    if(is.windows) parallel <- BiocParallel::SerialParam()
    if(is.null(parallel)) parallel <- BiocParallel::bpparam()

    fail <- TRUE
    trials <- 0
    while(fail && trials < 5)
    {
        res <- try(
            BiocParallel::bplapply(exp.list, f, ..., BPPARAM=parallel),
            silent=TRUE)
        fail <- is(res, "try-error")
        trials <- trials + 1
    }
    if(fail) stop(res)
    return(res)
}

#' Evaluation of the type I error rate of enrichment methods
#' 
#' This function evaluates the type I error rate of selected methods for 
#' enrichment analysis when applied to one or more expression datasets.
#' 
#' 
#' @param exp.list Experiment list.  A \code{list} of datasets, each being of
#' class \code{\linkS4class{SummarizedExperiment}}.
#' @param methods Methods for enrichment analysis.  A character vector with
#' method names chosen from \code{\link{sbeaMethods}} and
#' \code{\link{nbeaMethods}}, or user-defined functions
#' implementing methods for enrichment analysis.
#' @param gs Gene sets, i.e. a list of character vectors of gene IDs.
#' @param alpha Numeric. Statistical significance level. Defaults to 0.05.
#' @param ea.perm Integer. Number of permutations of the sample group assignments 
#' during enrichment analysis. Defaults to 1000. Can also be an integer vector 
#' matching the length of 'methods' to assign different numbers of permutations
#' for different methods.
#' @param tI.perm Integer. Number of permutations of the sample group assignments 
#' during type I error rate evaluation. Defaults to 1000. Can also be an integer
#' vector matching the length of \code{methods} to assign different numbers of
#' permutations for different methods.
#' @param perm.block.size Integer. When running in parallel, splits \code{tI.perm}
#' into blocks of the indicated size. Defaults to -1, which indicates to not  
#' partition \code{tI.perm}.
#' @param summarize Logical. If \code{TRUE} (default) applies \code{\link{summary}}
#' to the vector of type I error rates across \code{tI.perm} permutations of the 
#' sample labels. Use \code{FALSE} to return the full vector of type I error rates.
#' @param save2file Logical. Should results be saved to file for subsequent
#' benchmarking?  Defaults to \code{FALSE}.
#' @param out.dir Character.  Determines the output directory where results are
#' saved to.  Defaults to \code{NULL}, which then writes to
#' \code{rappdirs::user_data_dir("GSEABenchmarkeR")} in case \code{save2file}
#' is set to \code{TRUE}.
#' @param verbose Logical. Should progress be reported? Defaults to \code{TRUE}.
#' @param ...  Additional arguments passed to the selected enrichment methods.
#' @return A list with an entry for each method applied.  Each method entry is
#' a list with an entry for each dataset analyzed.  Each dataset entry is either
#' a summary (\code{summarize=TRUE}) or the type I error rates itself 
#' (\code{summarize=FALSE}) across \code{tI.perm} permutations of the sample labels. 
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\link{sbea}} and \code{\link{nbea}}
#' for carrying out set- and network-based enrichment analysis.
#' 
#' \code{\linkS4class{BiocParallelParam}} and \code{\link{register}} for
#' configuration of parallel computation.
#' @examples
#' 
#'     # loading three datasets from the GEO2KEGG compendium
#'     geo2kegg <- loadEData("geo2kegg", nr.datasets=3)
#'
#'     # only considering the first 1000 probes for demonstration
#'     geo2kegg <- lapply(geo2kegg, function(d) d[1:1000,]) 
#'
#'     # preprocessing and DE analysis for two of the datasets
#'     geo2kegg <- maPreproc(geo2kegg[2:3])
#'     geo2kegg <- runDE(geo2kegg)
#' 
#'     # getting a subset of human KEGG gene sets
#'     gs.file <- system.file("extdata", package="EnrichmentBrowser")
#'     gs.file <- file.path(gs.file, "hsa_kegg_gs.gmt") 
#'     kegg.gs <- EnrichmentBrowser::getGenesets(gs.file)
#' 
#'     # evaluating type I error rate of two methods on two datasets
#'     # NOTE: using a small number of permutations for demonstration;
#'     #       for a meaningful evaluation tI.perm should be >= 1000   
#'     res <- evalTypeIError(geo2kegg, methods=c("ora", 
#'              "camera"), gs=kegg.gs, ea.perm=0, tI.perm=3)
#'     
#' 
#' @export evalTypeIError
evalTypeIError <- function(methods, exp.list, gs, alpha=0.05, 
    ea.perm=1000, tI.perm=1000, perm.block.size=-1, summarize=TRUE, 
    save2file=FALSE, out.dir=NULL, verbose=TRUE, ...)
{
    # singleton call?
    if(is(exp.list, "SummarizedExperiment"))
    {
        if(length(methods) == 1)
        {
            res <- .evalTypeI(methods, exp.list, gs, alpha, 
                                ea.perm[1], tI.perm[1], 
                                perm.block.size=perm.block.size, 
                                summarize=summarize,    
                                save2file=save2file, out.dir=out.dir, ...)
            return(res)
        }
        exp.list <- list(exp.list)
    }
    
    # setup
    .eaPkgs(methods)
    GRP.COL <- EnrichmentBrowser::configEBrowser("GRP.COL")
    BLK.COL <- EnrichmentBrowser::configEBrowser("BLK.COL")
  
    # remove blocking
    for(i in seq_along(exp.list))
    {
        nind <- colnames(colData(exp.list[[i]])) != BLK.COL
        colData(exp.list[[i]]) <- colData(exp.list[[i]])[,nind]
    } 
    
    # different number of permutations for different methods?
    nr.meth <- length(methods)
    if(length(ea.perm) != nr.meth) ea.perm <- rep(ea.perm[1], nr.meth)   
    if(length(tI.perm) != nr.meth) tI.perm <- rep(tI.perm[1], nr.meth)
    if(length(perm.block.size) != nr.meth) 
        perm.block.size <- rep(perm.block.size[1], nr.meth)
    
    names(ea.perm) <- names(tI.perm) <- names(perm.block.size) <- methods   
    show.progress <- verbose && length(exp.list) > 2
    # iterate one method over multiple datasets 
    .iterD <- function(se, m, pb)
    {
        id <- metadata(se)$dataId
        if(show.progress) setTxtProgressBar(pb, match(id, names(exp.list)))
            
        .evalTypeI(m, se, gs, alpha, ea.perm[m], tI.perm[m], 
                    perm.block.size=perm.block.size, summarize=summarize,
                    save2file=save2file, out.dir=out.dir, ...)
    }

    # iterate over methods
    .iterM <- function(m)
    {
        if(verbose) message(m)
        pb <- NULL
        if(show.progress) pb <- txtProgressBar(0, length(exp.list), style=3)

        res <- lapply(exp.list, .iterD, m=m, pb=pb)
        names(res) <- names(exp.list)
        if(show.progress) close(pb)
        return(res)
    }

    res <- lapply(methods, .iterM) 
    names(res) <- methods
    return(res)
}

# for one method and one dataset
.evalTypeI <- function(method, se, gs, alpha=0.05, 
    ea.perm=1000, tI.perm=1000, perm.mat=NULL, perm.block.size=-1, 
    summarize=TRUE, save2file=FALSE, out.dir=NULL, ...)
{
    GRP.COL <- EnrichmentBrowser::configEBrowser("GRP.COL")
    PVAL.COL <- EnrichmentBrowser::configEBrowser("PVAL.COL")

    if(length(tI.perm) > 1) tI.perm <- tI.perm[method]
    if(length(ea.perm) > 1) ea.perm <- ea.perm[method]
    if(length(perm.block.size) > 1) perm.block.size <- perm.block.size[method]

    # is permutation matrix given as arg? 
    if(is.null(perm.mat)) perm.mat <- .getPermMat(se[[GRP.COL]], tI.perm)
    else if(tI.perm < ncol(perm.mat)) perm.mat <- perm.mat[,seq_len(tI.perm)]

    grid <- seq_len(ncol(perm.mat))
    .calcFPR <- function(i)
    {
        se[[GRP.COL]] <- perm.mat[,i]    
        # TODO: argument requireDE (for methods such as ora, ebm, ...)
        if(method %in% c("ora", "ebm")) 
            se <- EnrichmentBrowser::deAna(se, padj.method="none")
        res <- runEA(se, method, gs, ea.perm, ...)
        res <- res$ranking
        res <- mean(res[[PVAL.COL]] < alpha)
        return(res)
    }

    # parallel: one (or more) permutation per core
    serial <- perm.block.size < 0 || ncol(perm.mat) < 10
    if(serial) res <- vapply(grid, .calcFPR, numeric(1))
    else if(perm.block.size > 1)
    {
        # split into blocks of defined size 
        blocks <- seq(1, ncol(perm.mat), by=perm.block.size)
        bdiff <- perm.block.size - 1
        last <- blocks[length(blocks)]
        blocks <- lapply(blocks, function(b) 
            c(b, ifelse(b == last, ncol(perm.mat), b + bdiff)))
        .f <- function(b) vapply(b[1]:b[2], .calcFPR, numeric(1))
        res <- BiocParallel::bplapply(blocks, .f)
        res <- unlist(res)
    }
    else
    { 
        # each permutation in parallel
        res <- BiocParallel::bplapply(grid, .calcFPR) 
        res <- unlist(res)
    }

    if(summarize) res <- summary(res)
    if(save2file) .save2file(res, out.dir, method, metadata(se)$dataId)
    return(res)
}

.getPermMat <- function(grp, perm=1000)
{
    perm.mat <- replicate(perm, sample(grp))

    # any permutation identical to observed setup?
    ind.id <- apply(perm.mat, 2, function(x) identical(grp, x))
    if(any(ind.id)) perm.mat <- perm.mat[,!ind.id]

    # any duplicated permutations?
    ind.dup <- duplicated(data.frame(t(perm.mat)))
    if(any(ind.dup)) perm.mat <- perm.mat[,!ind.dup]

    d <- perm - ncol(perm.mat) 
    if(d) message(paste0("Removed ", d, 
                " identical / duplicated permutation",
                ifelse(d > 1, "s", "")))
    
    return(perm.mat)
}

#' Evaluation of enrichment methods on random gene sets
#' 
#' This function evaluates the proportion of rejected null hypotheses 
#' (= the fraction of significant gene sets) of an enrichment method  
#' when applied to random gene sets of defined size
#' 
#' 
#' @param method Enrichment analysis method.  A character scalar chosen 
#' from \code{\link{sbeaMethods}} and \code{\link{nbeaMethods}}, or a user-defined
#' function implementing a method for enrichment analysis.
#' @param se An expression dataset of class \code{\linkS4class{SummarizedExperiment}}.
#' @param nr.gs Integer. Number of random gene sets. Defaults to 100. 
#' @param set.size Integer. Gene set size, i.e. number of genes in each random gene set. 
#' @param alpha Numeric. Statistical significance level. Defaults to 0.05.
#' @param padj Character. Method for adjusting p-values to multiple testing.
#' For available methods see the man page of the stats function
#' \code{\link{p.adjust}}. Defaults to \code{"none"}.
#' @param perc Logical.  Should the percentage (between 0 and 100, default)
#' or the proportion (between 0 and 1) of significant gene sets be returned?
#' @param reps Integer. Number of replications. Defaults to 100.
#' @param rep.block.size Integer. When running in parallel, splits \code{reps}
#' into blocks of the indicated size. Defaults to -1, which indicates to not  
#' partition \code{reps}.
#' @param summarize Logical. If \code{TRUE} (default) returns the mean (\code{\link{mean}})
#' and the standard deviation (\code{\link{sd}}) of the proportion of significant
#' gene sets across \code{reps} replications. Use \code{FALSE} to return the full 
#' vector storing the proportion of significant gene sets for each replication.
#' @param save2file Logical. Should results be saved to file for subsequent
#' benchmarking?  Defaults to \code{FALSE}.
#' @param out.dir Character.  Determines the output directory where results are
#' saved to.  Defaults to \code{NULL}, which then writes to
#' \code{rappdirs::user_data_dir("GSEABenchmarkeR")} in case \code{save2file}
#' is set to \code{TRUE}.
#' @param ...  Additional arguments passed to the selected enrichment method.
#' @return A named numeric vector of length 2 storing mean and standard deviation
#' of the proportion of significant gene sets across \code{reps} replications 
#' (\code{summarize=TRUE}); or a numeric vector of length \code{reps} storing the
#' the proportion of significant gene sets for each replication itself 
#' (\code{summarize=FALSE}). 
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\link{sbea}} and \code{\link{nbea}}
#' for carrying out set- and network-based enrichment analysis.
#' 
#' \code{\linkS4class{BiocParallelParam}} and \code{\link{register}} for
#' configuration of parallel computation.
#' @examples
#' 
#'     # loading three datasets from the GEO2KEGG compendium
#'     geo2kegg <- loadEData("geo2kegg", nr.datasets=3)
#'
#'     # only considering the first 1000 probes for demonstration
#'     geo2kegg <- lapply(geo2kegg, function(d) d[1:1000,]) 
#'
#'     # preprocessing and DE analysis for two of the datasets
#'     geo2kegg <- maPreproc(geo2kegg[2:3])
#'     geo2kegg <- runDE(geo2kegg)
#' 
#'     evalRandomGS("camera", geo2kegg[[1]], reps=3)
#'     
#' 
#' @export evalRandomGS
evalRandomGS <- function(method, se, nr.gs=100, set.size=5, 
    alpha=0.05, padj = "none", perc=TRUE, reps=100, rep.block.size=-1, 
    summarize=TRUE, save2file=FALSE, out.dir=NULL, ...)
{
    PVAL.COL <- EnrichmentBrowser::configEBrowser("PVAL.COL")
    
    .eval <- function(i)
    {
        gs <- replicate(nr.gs, sample(names(se), set.size), simplify=FALSE)
        names(gs) <- paste0("gs", seq_len(nr.gs))
        res <- runEA(se, method, gs, ...)
        res <- res$ranking
        res <- p.adjust(res[,PVAL.COL], method=padj)
        res <- mean(res < alpha)
        if(perc) res <- round(res * 100, digits=2)
    }

    serial <- rep.block.size < 0 || reps < 10
    if(serial) res <- replicate(reps, .eval())
    else if(rep.block.size > 1)
    {
        # split into blocks of defined size 
        blocks <- seq(1, reps, by=rep.block.size)
        bdiff <- rep.block.size - 1
        last <- blocks[length(blocks)]
        blocks <- lapply(blocks, function(b) 
            c(b, ifelse(b == last, reps, b + bdiff)))
        .f <- function(b) vapply(b[1]:b[2], .eval, numeric(1))
        res <- BiocParallel::bplapply(blocks, .f)
        res <- unlist(res)
    }
    else
    { 
        # each permutation in parallel
        res <- BiocParallel::bplapply(seq_len(reps), .eval) 
        res <- unlist(res)
    }

    if(summarize) res <- c(mean=mean(res), sd=sd(res))
    if(save2file) .save2file(res, out.dir, method, paste0("gs", set.size))
    return(res)
}


#' Evaluating gene set rankings for the number of (significant) sets
#' 
#' These functions evaluate gene set rankings obtained from applying enrichment
#' methods to multiple datasets.  This allows to assess resulting rankings for
#' granularity (how many gene sets have a unique p-value?) and statistical
#' significance (how many gene sets have a p-value below a significance
#' threshold?).
#' 
#' 
#' @param ea.ranks Enrichment analysis rankings.  A list with an entry for each
#' enrichment method applied.  Each entry is a list that stores for each
#' dataset analyzed the resulting gene set ranking as obtained from applying
#' the respective method to the respective dataset.
#' @param uniq.pval Logical.  Should the number of gene sets with a unique
#' p-value or the total number of gene sets per ranking be returned?  Defaults
#' to \code{TRUE}.
#' @param perc Logical.  Should the percentage or absolute number of gene sets
#' be returned?  Percentage is typically more useful for comparison between
#' rankings with a potentially different total number of gene sets.  Defaults
#' to \code{TRUE}.
#' @param alpha Statistical significance level. Defaults to 0.05.
#' @param padj Character.  Method for adjusting p-values to multiple testing.
#' For available methods see the man page of the stats function
#' \code{\link{p.adjust}}.  Defaults to \code{BH}.
#' @return A list of numeric vectors storing for each method the number of
#' (significant) gene sets for each dataset analyzed.  If each element of the
#' resulting list is of equal length (corresponds to successful application of
#' each enrichment method to each dataset), the list is automatically
#' simplified to a numeric matrix (rows = datasets, columns = methods).
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\link{runEA}} to apply enrichment methods to multiple
#' datasets; \code{\link{readResults}} to read saved rankings as an input for
#' the eval-functions.
#' @examples
#' 
#'     # simulated setup:
#'     # 2 methods & 2 datasets
#'     methods <- paste0("m", 1:2)
#'     data.ids <- paste0("d", 1:2)
#' 
#'     # simulate gene set rankings
#'     ea.ranks <- sapply(methods, function(m) 
#'             sapply(data.ids, 
#'                 function(d)
#'                 {
#'                     r <- EnrichmentBrowser::makeExampleData("ea.res") 
#'                     r <- EnrichmentBrowser::gsRanking(r, signif.only=FALSE)
#'                     return(r)
#'                 }, simplify=FALSE),
#'                 simplify=FALSE)
#' 
#'     # evaluate
#'     evalNrSets(ea.ranks)
#'     evalNrSigSets(ea.ranks)
#' 
#' @export evalNrSigSets
evalNrSigSets <- function(ea.ranks, alpha=0.05, padj="none", perc=TRUE)
{
    pcol <- EnrichmentBrowser::configEBrowser("PVAL.COL")  
    res <- sapply(ea.ranks, 
        function(rs) vapply(rs, 
            function(r)
            { 
                nr <- sum(p.adjust(r[,pcol], method=padj) < alpha) 
                if(perc) nr <- nr / nrow(r) * 100
                return(nr)
            }, numeric(1)))
    return(res)
}

#' @rdname evalNrSigSets
#' @export
evalNrSets <- function(ea.ranks, uniq.pval=TRUE, perc=TRUE)
{ 
    pcol <- EnrichmentBrowser::configEBrowser("PVAL.COL")  
    res <- sapply(ea.ranks, 
		function(rs) vapply(rs, 
			function(r) 
            {
                nr <- ifelse(uniq.pval, length(unique(r[,pcol])), nrow(r)) 
                if(perc) nr <- nr / nrow(r) * 100 
                return(nr)
            }, numeric(1)))
	return(res)
}


#' Evaluating phenotype relevance of gene set rankings
#' 
#' This function evaluates gene set rankings obtained from the application of
#' enrichment methods to multiple datasets - where each dataset investigates a
#' certain phenotype such as a disease.  Given pre-defined phenotype relevance
#' scores for the gene sets, indicating how important a gene set is for the
#' investigated phenotype (as e.g. judged by evidence from the literature),
#' this allows to assess whether enrichment methods produce gene set rankings
#' in which phenotype-relevant gene sets accumulate at the top.
#' 
#' The function \code{evalRelevance} evaluates the similarity of a gene set ranking
#' obtained from enrichment analysis and a gene set ranking based on phenotype
#' relevance scores. Therefore, the function first transforms the ranks 'r'
#' from the enrichment analysis to weights 'w' in [0,1] via w = 1 - r / N;
#' where 'N' denotes the total number of gene sets on which the enrichment
#' analysis has been carried out. These weights are then multiplied with the
#' corresponding relevance scores and summed up.
#' 
#' The function \code{compOpt} applies \code{evalRelevance} to the theoretically
#' optimal case in which the enrichment analysis ranking is identical to the
#' relevance score ranking. The ratio between observed and optimal score is
#' useful for comparing observed scores between datasets / phenotypes.
#' 
#' The function \code{compRand} repeatedly applies \code{evalRelevance} to randomly
#' drawn gene set rankings to assess how likely it is to observe a score equal
#' or greater than the one obtained.
#' 
#' @param ea.ranks Enrichment analysis rankings.  A list with an entry for each
#' enrichment method applied.  Each entry is a list that stores for each
#' dataset analyzed the resulting gene set ranking obtained from applying the
#' respective method to the respective dataset.  Resulting gene set rankings
#' are assumed to be of class \code{\linkS4class{DataFrame}} in which gene sets
#' (required column named \code{GENE.SET}) are ranked according to a ranking
#' measure such as a gene set p-value (required column named \code{P.VALUE}).
#' See \code{\link{gsRanking}} for an example.
#' @param rel.ranks Relevance score rankings.  A list with an entry for each
#' phenotype investigated.  Each entry should be a
#' \code{\linkS4class{DataFrame}} in which gene sets (rownames are assumed to
#' be gene set IDs) are ranked according to a phenotype relevance score
#' (required column \code{REL.SCORE}).
#' @param data2pheno A named character vector where the names correspond to
#' dataset IDs and the elements of the vector to the corresponding phenotypes
#' investigated.
#' @param perc Logical.  Should observed scores be returned as-is or as a
#' *perc*entage of the respective optimal score. Percentages of the optimal
#' score are typically easier to interpret and are comparable between datasets
#' / phenotypes.  Defaults to \code{TRUE}.
#' @param top Integer.  If \code{top} is non-zero, the evaluation will be
#' restricted to the first \code{top} gene sets of each enrichment analysis
#' ranking.  Defaults to \code{0}, which will then evaluate the full ranking.
#' @param rand Logical.  Should gene set rankings be randomized to assess how
#' likely it is to observe a score equal or greater than the respective
#' obtained score?  Defaults to \code{FALSE}.
#' @param perm Integer. Number of permutations if \code{rand} set to \code{TRUE}.
#' @param gs.ids Character vector of gene set IDs on which enrichment analysis 
#' has been carried out.
#' @return A numeric matrix (rows = datasets, columns = methods) storing in
#' each cell the relevance score sum obtained from applying the respective
#' method to the respective dataset.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{runEA} to apply enrichment methods to multiple datasets;
#' \code{readResults} to read saved rankings as an input for the eval-functions;
#' @examples
#' 
#'     #
#'     # (1) simulated setup: 1 enrichment method applied to 1 dataset
#'     #  
#'
#'     # simulate gene set ranking
#'     ea.ranks <- EnrichmentBrowser::makeExampleData("ea.res")
#'     ea.ranks <- EnrichmentBrowser::gsRanking(ea.ranks, signif.only=FALSE)
#' 
#'     # simulated relevance score ranking
#'     rel.ranks <- ea.ranks
#'     rel.ranks[,2] <- runif(nrow(ea.ranks), min=1, max=100)
#'     colnames(rel.ranks)[2] <- "REL.SCORE"
#'     rownames(rel.ranks) <- rel.ranks[,"GENE.SET"]
#'     ind <- order(rel.ranks[,"REL.SCORE"], decreasing=TRUE)
#'     rel.ranks <- rel.ranks[ind,]
#' 
#'     # evaluate
#'     evalRelevance(ea.ranks, rel.ranks)    
#'     compOpt(rel.ranks, ea.ranks[,"GENE.SET"])
#'     compRand(rel.ranks, ea.ranks[,"GENE.SET"], perm=3)
#' 
#'     # 
#'     # (2) simulated setup: 2 methods & 2 datasets
#'     #
#'        
#'     methods <- paste0("m", 1:2)
#'     data.ids <- paste0("d", 1:2)
#' 
#'     # simulate gene set rankings
#'     ea.ranks <- sapply(methods, function(m) 
#'             sapply(data.ids, 
#'                 function(d)
#'                 {
#'                     r <- EnrichmentBrowser::makeExampleData("ea.res") 
#'                     r <- EnrichmentBrowser::gsRanking(r, signif.only=FALSE)
#'                     return(r)
#'                 }, simplify=FALSE),
#'                 simplify=FALSE)
#' 
#'     # simulate a mapping from datasets to disease codes
#'     d2d <- c("ALZ", "BRCA")
#'     names(d2d) <- data.ids
#' 
#'     # simulate relevance score rankings
#'     rel.ranks <- lapply(ea.ranks[[1]],
#'         function(rr)
#'         {
#'             rr[,2] <- runif(nrow(rr), min=1, max=100)
#'             colnames(rr)[2] <- "REL.SCORE"
#'             rownames(rr) <- rr[,"GENE.SET"]
#'             ind <- order(rr[,"REL.SCORE"], decreasing=TRUE)
#'             rr <- rr[ind,]
#'             return(rr)
#'         })
#'     names(rel.ranks) <- unname(d2d)
#' 
#'     # evaluate
#'     evalRelevance(ea.ranks, rel.ranks, d2d)
#' 
#' @export evalRelevance
evalRelevance <- function( ea.ranks, rel.ranks, 
                            data2pheno, perc=TRUE, top=0, rand=FALSE )
{
    # singleton call?
    is.singleton <- is(ea.ranks, "DataFrame") && is(rel.ranks, "DataFrame")
    if(is.singleton) return( .relScore(ea.ranks, rel.ranks, top) )

    # iterating over enrichment methods included
    res <- lapply(ea.ranks, 
        function(mranks)
        {   
            # iterating over datasets included
            meval.res <- vapply( names(mranks), 
                function(d) .deployScore(d, mranks, rel.ranks, 
                    data2pheno, top, type="rel"), numeric(1) )
            return(meval.res)
        }) 

    # postprocess missing
    for(i in seq_along(res)) res[[i]] <- res[[i]][names(data2pheno)]
    x <- do.call(cbind, res)

    # get analyzed gene set IDs
    ind <- which.max(lengths(ea.ranks))
    .fspl <- function(n) unlist(strsplit(n, "_"))[1]
    gs.ids <- lapply(ea.ranks[[ind]], function(r) 
        vapply(r[,"GENE.SET"], .fspl, character(1), USE.NAMES=FALSE))
    
    names(gs.ids) <- names(ea.ranks[[ind]])
    opt <- compOpt(rel.ranks, gs.ids, data2pheno, top)
    opt <- opt[names(data2pheno)]

    # percentage of optimal score ...
    if(perc) x <- x / opt * 100

    # ... or absolute score?
    else
    {
        x <- cbind(x, opt)
        colnames(x)[ncol(x)] <- "opt"
    }

    # include random scores?
    if(rand)
    {
        rscores <- compRand(rel.ranks, gs.ids, data2pheno)
        if(perc) rscores <- rscores / opt * 100
        x <- cbind(x, rscores)
        colnames(x)[ncol(x)] <- "rand"
    }
    return(x)
}

#evalROC <- function(ea.ranks, rel.ranks, gs, top=10)
#{
#    # create random gene sets matching the rel sets
#    rel.sets <- rownames(rel.ranks)
#    gs.ids <- vapply(names(gs), 
#                        function(s) unlist(strsplit(s, "_"))[1],
#                        character(1), USE.NAMES=FALSE)
#    rgs <- relist(sample(unlist(gs)), skeleton=gs)    
#    rgs <- rgs[gs.ids %in% rel.sets]
#    names(rgs) <- paste0("random", seq_along(rgs))
#
#    # perform EA on the intermixture of real and random GS 
#    res <- sbea("ora", se=geo2kegg$GSE1297, gs=c(gs,rgs), perm=0)
#    res <- res$res.tbl
#    gs.ids <- vapply(res$GENE.SET, 
#                        function(s) unlist(strsplit(s, "_"))[1],
#                        character(1), USE.NAMES=FALSE)
#    ind <- gs.ids %in% c(rel.sets, names(rgs))
#    res <- res[ind,]
#    gs.ids <- gs.ids[ind]
#
#    # handmade ROC
#    fp <- function(i) sum(!(gs.ids[seq_len(i)] %in% rel.sets[seq_len(top)])) / (nrow(res)-top)
#    tp <- function(i) sum(gs.ids[seq_len(i)] %in% rel.sets[seq_len(top)]) / top
#    fps <- vapply(seq_len(nrow(res)), fp, numeric(1))
#    tps <- vapply(seq_len(nrow(res)), tp, numeric(1))
#    o <- order(fps)
#    plot(fps[o], tps[o], type="l")
#
#    # using ROCR
#    rscores <- vapply(ea.methods, 
#        function(m) evalRelevance(ma.kegg.ranks[[m]]$GSE1297, mala.kegg$ALZ), 
#        numeric(1))
#
#    auc <- function(ea.ranks, rel.ranks)
#    {
#        rel.sets <- rownames(rel.ranks)
#        r <- 1 - EnrichmentBrowser:::.getRanks(ea.ranks) / 100
#        gs.ids <- vapply(ea.ranks$GENE.SET, 
#                        function(s) unlist(strsplit(s, "_"))[1],
#                        character(1), USE.NAMES=FALSE)
#        labels <- gs.ids %in% rel.sets[1:10]
#        pr <- ROCR::prediction(r, ifelse(labels,1,0))
#        auc <- ROCR::performance(pr, "auc")
#        auc <- unlist(auc@y.values)    
#        return(auc)
#    }
#
#    aucs <- vapply(ea.methods, 
#        function(m) auc(ma.kegg.ranks[[m]]$GSE1297, mala.kegg$ALZ), 
#        numeric(1))
#
#}

.relScore <- function(ea.ranks, rel.ranks, top=0)
{
    ids <- vapply(ea.ranks[, "GENE.SET"], 
                    function(n) unlist(strsplit(n, "_"))[1],
                    character(1), 
                    USE.NAMES=FALSE)
    
    r <- 1 - EnrichmentBrowser:::.getRanks(ea.ranks) / 100
    names(r) <- ids
    if(top) r <- r[seq_len(top)]
    r <- r[rownames(rel.ranks)]
    score <- sum(r * rel.ranks[, "REL.SCORE"], na.rm=TRUE)
    return(score)
}

.deployScore <- function(d, mdat, rel.ranks, 
    data2pheno, add, type=c("rel", "opt", "rand"))
{
    type <- match.arg(type)     
     
    ear <- mdat[[d]]
    d <- data2pheno[d]
    mar <- rel.ranks[[d]]

    if(type == "rel") r <- .relScore(ear, mar, add)
    else if(type == "opt") r <- .optScore(mar, ear, add)
    else r <- median( .randScore(mar, ear, add) )
    return(r)
}


# compute optimal score for relevance rankings
#' @rdname evalRelevance
#' @export
compOpt <- function(rel.ranks, gs.ids, data2pheno=NULL, top=0)
{
    # singleton call?
    if(is(rel.ranks, "DataFrame"))
        return( .optScore(rel.ranks, gs.ids, top) )

    opt <- vapply( names(gs.ids), function(d)
            .deployScore(d, gs.ids, rel.ranks, 
                            data2pheno, top, type="opt"), 
            numeric(1) )
    
    return(opt)
}

# compute optimal score for one particular relevance ranking, eg ALZ
.optScore <- function(rel.ranks, gs.ids, top=0)
{
    # compute opt
    dummy.gsp <- seq_along(gs.ids) / length(gs.ids)
    scol <- EnrichmentBrowser::configEBrowser("GS.COL")
    pcol <- EnrichmentBrowser::configEBrowser("PVAL.COL")

    mar <- rel.ranks
    marn <- rownames(mar)[rownames(mar) %in% gs.ids]
    optGS <- c(marn, gs.ids[!(gs.ids %in% marn)])
    optR <- DataFrame(optGS, dummy.gsp)
    colnames(optR) <- c(scol, pcol)
    r <- .relScore(optR, mar, top)
    return(r)    
}

# compute random score for relevance rankings
#' @rdname evalRelevance
#' @export
compRand <- function(rel.ranks, gs.ids, data2pheno=NULL, perm=1000)
{
    # singleton call?
    if(is(rel.ranks, "DataFrame"))
        return( .randScore(rel.ranks, gs.ids, perm) )
 
    rand <- vapply( names(gs.ids), 
                function(d)
                {
                    message(d) 
                    .deployScore(d, gs.ids, rel.ranks, 
                                    data2pheno, perm, type="rand")
                }, numeric(1) )

    return(rand)
}
 
# compute random score for one particular relevance ranking, eg ALZ
.randScore <- function(rel.ranks, gs.ids, perm=1000)
{
    gs.ids <- sort(gs.ids)
    dummy.gsp <- seq_along(gs.ids) / length(gs.ids)
    scol <- EnrichmentBrowser::configEBrowser("GS.COL")
    pcol <- EnrichmentBrowser::configEBrowser("PVAL.COL")

    f <- function()
    {
        randGS <- sample(gs.ids)
        randR <- DataFrame(randGS, dummy.gsp)
        colnames(randR) <- c(scol, pcol)
        r <- .relScore(randR, rel.ranks)
        return(r)
    }

    rand <- replicate(perm, f())
    return(rand) 
}


# analyse pheno label shuffling for datasets
.shuffleD2D <- function(ea.ranks, rel.ranks)
{
    gs <- vapply(ea.ranks[,"GENE.SET"], 
                    function(n) unlist(strsplit(n, "_"))[1],
                    character(1), USE.NAMES=FALSE)

    res <- vapply(rel.ranks, 
        function(m)
        {
            x <- evalRelevance(ea.ranks, m)
            opt <- compOpt(m, gs)
            return(x/opt)
        }, numeric(1))

    return(res)
}

# check how well relevance rankings recover each other
.rel2rel <- function(rel.ranks, gs)
{
    dummy.gsp <- seq_along(gs) / length(gs)
    scol <- EnrichmentBrowser::configEBrowser("GS.COL")
    pcol <- EnrichmentBrowser::configEBrowser("PVAL.COL")

    m2m <- vapply(names(rel.ranks),
        function(m)
        {
            curr <- rel.ranks[[m]]
            opt <- compOpt(curr, gs)
            res <- vapply(names(rel.ranks),
                function(n)
                {
                    mar <- rel.ranks[[n]]
                    marn <- rownames(mar)[rownames(mar) %in% gs]
                    optGS <- c(marn, gs[!(gs %in% marn)])
                    optR <- DataFrame(optGS, dummy.gsp)
                    colnames(optR) <- c(scol, pcol)
                    r <- evalRelevance(optR, curr)
                    return(r / opt)    
                }, numeric(1))
            return(res)
        }, numeric(length(rel.ranks))) 
}

#
# MAIN
#
benchmark <- function(
    method, 
    edata=c("geo2kegg", "tcga"),
    dataType=c("ma", "rseq"), 
    deMethod=c("limma", "edgeR", "DESeq"),
    gs=c("kegg", "go_bp", "go_mf"),
    grn=c("kegg", "encode"), 
    eaType=c("sbea", "nbea", "both"), 
    org="hsa", 
    parallel=NULL,
    ...)
{
    dataType <- match.arg(dataType)
    eaType <- match.arg(eaType)

    # expression data
    message("Loading expression data ...")
    edata <- loadEData(edata[1], ...)

    message("DE analysis ...")
    edata <- runDE(edata, deMethod, ...)

    # gene sets
    if(is.character(gs))
    {
        message("Loading gene sets ...")
        gs <- match.arg(gs)
        if(gs == "kegg") gs <- EnrichmentBrowser::getGenesets(org, db=gs)
        else EnrichmentBrowser::getGenesets(org, go.onto=toupper(substring(gs,4,5)))
    }

    # enrichment analysis
    message("Executing EA ...")

    # method to execute
    methods <- switch(eaType,
                        sbea = EnrichmentBrowser::sbeaMethods(),
                        nbea = EnrichmentBrowser::nbeaMethods(),
                        c(EnrichmentBrowser::sbeaMethods(), 
                            EnrichmentBrowser::nbeaMethods())) 

    if(!missing(method)) methods <- union(method, methods)   

    if(any(methods %in% EnrichmentBrowser::nbeaMethods()))
    {
        # gene regulatory network
        message("Loading gene regulatory network")
        if(grn == "kegg") grn <- EnrichmentBrowser::compileGRN(org, db="kegg")
        # TODO: else
    }

    res <- runEA(edata, methods, gs, ...)      
    return(res)
}


