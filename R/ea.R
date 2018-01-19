#
# DE analysis
#
deApply <- function(data.ids, 
    de.method=c("limma", "edgeR", "DESeq"), 
    padj.method="flexible", parallel=NULL, ...)
{
    flex <- padj.method == "flexible"
    if(flex) padj.method <- "none"
    
    ADJP.COL <- EnrichmentBrowser::config.ebrowser("ADJP.COL")
    .f <- function(i)
    { 
        se <- EnrichmentBrowser::de.ana(i, de.method=de.method, 
                                            padj.method=padj.method, ...)

        if(flex)
        {
            padj <- p.adjust(rowData(se)[,ADJP.COL], method="BH")
            fracSigP <- mean(padj < 0.05)
            if(fracSigP > 0.25) 
            padj <- p.adjust(rowData(se)[,ADJP.COL], method="bonf") 
            if(fracSigP > 0.01) rowData(se)[,ADJP.COL] <- padj
        }
        return(se)
    }
    
    esets <- .iter(.f, data.ids, parallel)
    for(d in data.ids) .EDATA[[d]] <- esets[[d]]
}

.fractDE <- function(eset, alpha=0.05, beta=1, freq=c("rel", "abs"))
{
    freq <- match.arg(freq)
    FC.COL <- EnrichmentBrowser::config.ebrowser("FC.COL")
    ADJP.COL <- EnrichmentBrowser::config.ebrowser("ADJP.COL")

    fract.p <- sum(rowData(eset)[,ADJP.COL] < alpha)
    fract.fc <- sum(abs(rowData(eset)[,FC.COL]) > beta)
    
    if(freq == "rel")
    { 
        fract.p <- fract.p / nrow(eset)
        fract.fc <- fract.fc / nrow(eset)
    }
    return(list(fp=fract.p, ffc=fract.fc))
}

#
# Enrichment analyis
#

eaApply <- function(data.ids, methods, gs, 
    parallel=NULL, save2file=FALSE, out.dir=NULL, ...)
{ 
    res <- lapply(methods, 
        function(m)
        {
            message(m)
            f <- function(i) runEA(i, 
                method=m, gs=gs, save2file=save2file, out.dir=out.dir, ...)
            r <- .iter(f, data.ids, parallel)
            return(r)
        }
    )    
    names(res) <- methods
    return(res)
}

runEA <- function(se, method, gs, save2file=FALSE, out.dir=NULL, ...)
{
   
    id <- metadata(se)$dataId
    f <- file()
    sink(file=f)

    # execute
    suppressMessages(
    if(method %in% EnrichmentBrowser::sbea.methods())
    {
        # TODO: check effect when 'perm' is given in '...' 
        perm <- ifelse(method=="ora", 0, 1000)
        ti <- system.time(
            res <- try(
                EnrichmentBrowser::sbea(method, se, gs, perm=perm, ...),
                silent=TRUE
            )
        )
    }
    else
    { 
        ti <- system.time(
            res <- try(EnrichmentBrowser::nbea(method, se, gs, ...), silent=TRUE)
        )
    })
    sink() ## undo silencing
    close(f)
    if(is(res, "try-error"))
    {
        message(paste(method, "could not be evaluated on", id))
        message(res)
        message("Returning NULL")
        return(NULL)
    }

    # output
    res <- EnrichmentBrowser::gs.ranking(res, signif.only=FALSE)

    if(save2file)
    {
        if(is.null(out.dir)) 
            out.dir <- rappdirs::user_data_dir("GSEABenchmarkeR")
        if(!file.exists(out.dir)) dir.create(out.dir)
        out.dir <- file.path(out.dir, method)
        if(!file.exists(out.dir)) dir.create(out.dir)

        out.file <- file.path(out.dir, paste0(id, ".rds"))
        saveRDS(res, file=out.file)
        cat(ti[3], file=sub("rds$", "txt", out.file))
    }
    return(list(runtime=ti[3], ranking=res))
}

