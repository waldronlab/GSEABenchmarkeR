.EDATA <- new.env(parent=emptyenv())

getDataset <- function(data.id) .EDATA[[data.id]]

#
# loading expression data 
#
loadEData <- function(edata, reload=FALSE, nr.datasets=NULL, ...)
{
    if(edata=="geo2kegg") data.ids <- .loadGEO2KEGG(reload, nr.datasets, ...)
    else if(edata=="tcga") data.ids <- .loadTCGA(reload, nr.datasets, ...)
    else if(file.exists(edata)) 
        data.ids <- .loadEDataFromFile(edata, reload, nr.datasets)
    else stop("edata must be \'geo2kegg\', \'tcga\' or an absolute file path")
    return(data.ids)
}


#
# loading expression data from file
#
.loadEDataFromFile <- function(edata, reload=FALSE, nr.datasets=NULL)
{
    data.ids <- sub(".rds$", "", list.files(edata, pattern="*.rds$"))
    if(!is.null(nr.datasets)) data.ids <- data.ids[seq_len(nr.datasets)]
    is.loaded <- all(vapply(data.ids, exists, logical(1), envir=.EDATA))
    if(!is.loaded || reload)
    {
        # message("Loading user-defined data compendium ...")
        # message(paste("... from", edata))
        if(interactive())
            pb <- txtProgressBar(1, length(data.ids), style=3, width=length(data.ids))
        for(did in data.ids)
        {
            dfile <- file.path(edata, paste(did, "rds", sep="."))
            se <- readRDS(dfile)
            metadata(se)$dataId <- did 
            .EDATA[[did]] <- se
            if(interactive()) setTxtProgressBar(pb, match(did, data.ids))
        }
        if(interactive()) close(pb)
    }
    return(data.ids)
}

#
# write differential expression per dataset
#
writeDE <- function(data.ids, out.dir=NULL)
{
    if(is.null(out.dir))
    {
        out.dir <- rappdirs::user_data_dir("GSEABenchmarkeR")
        if(!file.exists(out.dir)) 
        out.dir <- file.path(out.dir, "de")
    }
    if(!file.exists(out.dir)) dir.create(out.dir)
    .f <- function(xx)
    {
        gt.file <- paste0(metadata(xx)$dataId, ".txt")
        gt.file <- file.path(out.dir, gt.file)
        gt <- EnrichmentBrowser:::.getGeneAnno(names(xx), "hsa")
        gt <- cbind(gt, rowData(xx, use.names=TRUE))
        gt <- EnrichmentBrowser:::.sortGeneTable(gt)  
        ind <- gt[,config.ebrowser("EZ.COL")]      
        write.table(gt, file=gt.file, row.names=FALSE, quote=FALSE, sep="\t")
    }
    .iter(.f, data.ids)
}

#
# RUNTIME
#
readRuntime <- function(data.dir, data.ids, methods)
{
    rts <- lapply(methods, 
        function(m)
        {
            mdir <- file.path(data.dir, m)
            mfiles <- paste0(data.ids, ".txt")
            mfiles <- file.path(mdir, mfiles) 
            mrts <- vapply(mfiles, .readRuntimeFile, numeric(1), USE.NAMES=FALSE) 
            mrts <- mrts[!is.na(mrts)] 
            return(mrts)
        })
    names(rts) <- methods
    return(rts)
}

.readRuntimeFile <- function(f)
{
    r <- NA
    if(file.exists(f))
    { 
        r <- scan(f, quiet=TRUE)
        if(length(r) != 1 || !is.numeric(r)) r <- NA
    }
    
    if(is.na(r))
    {
        m <- basename(dirname(f))
        d <- sub(".txt$", "", basename(f))
        wmsg <- paste("Runtime of", m, "for", d, "not found")
        warning(wmsg)
    }
    return(r)
}

#
# RANKINGS
#

# read resulting rankings of a list of enrichment methods
readRankings <- function(data.dir, data.ids, methods)
{
    rts <- lapply(methods, 
        function(m)
        {
            mdir <- file.path(data.dir, m)
            mfiles <- paste0(data.ids, ".rds")
            mfiles <- file.path(mdir, mfiles) 
            r <- lapply(mfiles, 
                function(f)
                {
                    if(file.exists(f)) x <- readRDS(f)
                    else 
                    {
                        d <- sub(".rds$", "", basename(f))
                        wmsg <- paste("Ranking of", m, "for", d, "not found")
                        warning(wmsg)
                        x <- NULL
                    }
                    return(x)
                })
            names(r) <- sub(".rds$", "", basename(mfiles))
            r <- r[!vapply(r, is.null, logical(1))]
            return(r)
        })
    names(rts) <- methods
    return(rts)
}

# reads a mapping from datasetID (eg GSExxx) to diseaseID (eg ALZ)
# as eg required for the Tarca microarray compendium
readDataId2diseaseCodeMap <- function(map.file)
{
    cont <- readLines(map.file)
    abbr <- nam <- rep("", length=length(cont))
    for(i in seq_along(cont))
    {
        spl <- unlist(strsplit(cont[i], " "))
        nam[i] <- spl[1]
        abbr[i] <- spl[length(spl)]
        abbr[i] <- sub("\\[", "", abbr[i])
        abbr[i] <- sub("\\]", "", abbr[i])
    }
    names(abbr) <- nam
    return(abbr)
}


