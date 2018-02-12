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
.iter <- function(exp.list, f, parallel=NULL)
{
    if(.Platform$OS.type == "windows") 
        parallel <- BiocParallel::SerialParam()
    if(!is.null(parallel)) BiocParallel::register(parallel)
    
    res <- BiocParallel::bplapply(exp.list, f)
    return(res)
}

#
# SIGNIFICANT SETS   
#

# total number of gene sets evaluated by enrichment methods across datasets  
evalNrSets <- function(ea.ranks, uniq.pval=TRUE, perc=TRUE)
{ 
    pcol <- EnrichmentBrowser::config.ebrowser("GSP.COL")  
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

# number of significant gene sets returned by enrichment methods across datasets
evalNrSigSets <- function(ea.ranks, alpha=0.05, padj="none", perc=TRUE)
{
    pcol <- EnrichmentBrowser::config.ebrowser("GSP.COL")  
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


#
# DISEASE RELEVANCE: MALACARDS
#

# if <what>=="score":
#   returns weighted relevance relevance score sum
#   for 1 enrichment method applied to 1 dataset
#
# alternatively <what> can be in [0,1]
#   then it will be how many wrong hits will be ranked before 
#   <what>*100% of the rel-ranking  
evalRanks <- function(ea.ranks, rel.ranks, top=0)
{
    ids <- vapply(ea.ranks[, "GENE.SET"], 
                    function(n) unlist(strsplit(n, "_"))[1],
                    character(1), 
                    USE.NAMES=FALSE)
    
    # eval rel score
    #if(what== "score")
    #{
        r <- 1 - EnrichmentBrowser:::.getRanks(ea.ranks) / 100
        names(r) <- ids
        if(top) r <- r[seq_len(top)]
        score <- sum(r[rownames(rel.ranks)] * rel.ranks[, "REL.SCORE"], na.rm=TRUE)
        return(score)
    #}
    # count wrongs
    #else
    #{
    #    is.rel <- ids %in% rownames(rel.ranks)
    #    cs.rel <- cumsum(is.rel)
    #    nr.rights <- round(what * nrow(rel.ranks))
    #    nr.wrongs <- which(cs.rel == nr.rights)[1] - 1 
    #    return(nr.wrongs)
    #}
}


# returns weighted relevance relevance score sum
# for a list of enrichment methods applied to a list of datasets
# data2pheno: a mapping from datasetID (eg GSExxx) to phenoID (eg ALZ)
evalApply <- function( ea.ranks, rel.ranks, 
                        data2pheno, perc=TRUE, top=0, rand=FALSE)
{
    # iterating over enrichment methods included
    res <- lapply(ea.ranks, 
        function(mranks)
        {   
            # iterating over datasets included
            meval.res <- vapply(names(mranks), 
                function(d)
                {
                    ear <- mranks[[d]]
                    d <- data2pheno[d]
                    mar <- rel.ranks[[d]]
                    eval.res <- evalRanks(ear, mar, top=top)
                    return(eval.res)
                }, numeric(1))
            return(meval.res)
        }) 

    # postprocess missing
    for(i in seq_along(res)) res[[i]] <- res[[i]][names(data2pheno)]
    x <- do.call(cbind, res)

    # 
    ind <- which.max(lengths(ea.ranks))
    gs <- lapply(ea.ranks[[ind]],
            function(r) vapply(r[,"GENE.SET"], 
                function(n) unlist(strsplit(n, "_"))[1], 
                character(1), USE.NAMES=FALSE))
    names(gs) <- names(ea.ranks[[ind]])
    opt <- optApply(gs, data2pheno, rel.ranks, top)
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
        rscores <- randApply(gs, data2pheno, rel.ranks)
        if(perc) rscores <- rscores / opt * 100
        x <- cbind(x, rscores)
        colnames(x)[ncol(x)] <- "rand"
    }
    return(x)
}

# compute optimal score for one particular relevance ranking, eg ALZ
compOpt <- function(rel.ranks, gs, top=0)
{
    # compute opt
    dummy.gsp <- seq_along(gs) / length(gs)
    scol <- EnrichmentBrowser::config.ebrowser("GS.COL")
    pcol <- EnrichmentBrowser::config.ebrowser("GSP.COL")

    mar <- rel.ranks
    marn <- rownames(mar)[rownames(mar) %in% gs]
    optGS <- c(marn, gs[!(gs %in% marn)])
    optR <- DataFrame(optGS, dummy.gsp)
    colnames(optR) <- c(scol, pcol)
    r <- evalRanks(optR, mar, top=top)
    return(r)    
}


# compute optimal score for all datasets and corresponding relevance rankings
optApply <- function(gs, data2pheno, rel.ranks, top)
{
    opt <- vapply(names(gs), 
        function(n) 
        {
            s <- gs[[n]]
            n <- data2pheno[n]
            mar <- rel.ranks[[n]]
            r <- compOpt(mar, s, top)    
            return(r)    
        }, numeric(1))
    return(opt)
}

# compute random score for one particular relevance ranking, eg ALZ
compRand <- function(rel.ranks, gs, n=1000)
{
    gs <- sort(gs)
    dummy.gsp <- seq_along(gs) / length(gs)
    scol <- EnrichmentBrowser::config.ebrowser("GS.COL")
    pcol <- EnrichmentBrowser::config.ebrowser("GSP.COL")

    f <- function()
    {
        randGS <- sample(gs)
        randR <- DataFrame(randGS, dummy.gsp)
        colnames(randR) <- c(scol, pcol)
        r <- evalRanks(randR, rel.ranks)
        return(r)
    }
    rand <- replicate(n,f())
    return(rand) 
}

# compute random score for all datasets and corresponding relevance rankings
randApply <- function(gs, data2pheno, rel.ranks, n=1000)
{
    rand <- vapply(names(gs), 
        function(n) 
        {
            message(n)
            s <- gs[[n]]
            n <- data2pheno[n]
            mar <- rel.ranks[[n]]
            r <- compRand(mar, s, n)    
            return(median(r))    
        }, numeric(1))
    return(rand)
}
 
# analyse pheno label shuffling for datasets
shuffleD2D <- function(ea.ranks, rel.ranks)
{
    gs <- vapply(ea.ranks[,"GENE.SET"], 
                    function(n) unlist(strsplit(n, "_"))[1],
                    character(1), USE.NAMES=FALSE)

    res <- vapply(rel.ranks, 
        function(m)
        {
            x <- evalRanks(ea.ranks, m)
            opt <- compOpt(m, gs)
            return(x/opt)
        }, numeric(1))

    return(res)
}

# check how well relevance rankings recover each other
rel2rel <- function(rel.ranks, gs)
{
    dummy.gsp <- seq_along(gs) / length(gs)
    scol <- EnrichmentBrowser::config.ebrowser("GS.COL")
    pcol <- EnrichmentBrowser::config.ebrowser("GSP.COL")

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
                    r <- evalRanks(optR, curr)
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
    data.ids <- loadEData(edata)

    message("DE analysis ...")
    deApply(data.ids, deMethod)

    # gene sets
    if(is.character(gs))
    {
        message("Loading gene sets ...")
        gs <- match.arg(gs)
        gs <- ifelse(gs == "kegg",
            EnrichmentBrowser::get.kegg.genesets(org),
            EnrichmentBrowser::get.go.genesets(org, toupper(substring(gs,4,5))))
    }

    # enrichment analysis
    message("Executing EA ...")

    # method to execute
    methods <- switch(eaType,
                        sbea = EnrichmentBrowser::sbea.methods(),
                        nbea = EnrichmentBrowser::nbea.methods(),
                        c(EnrichmentBrowser::sbea.methods(), 
                            EnrichmentBrowser::nbea.methods())) 

    if(!missing(method)) methods <- union(method, methods)   

    if(any(methods %in% EnrichmentBrowser::nbea.methods()))
    {
        # gene regulatory network
        message("Loading gene regulatory network")
        if(grn == "kegg") grn <- EnrichmentBrowser::compile.grn.from.kegg(org)
        # TODO: else
    }

    res <- eaApply(data.ids, methods, gs, ...)      
    return(res)
}


