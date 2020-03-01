#' Loading pre-defined and user-defined expression data
#' 
#' This function implements a general interface for loading the pre-defined
#' GEO2KEGG microarray compendium and the TCGA RNA-seq compendium.  It also
#' allows loading of user-defined data from file.
#' 
#' The pre-defined GEO2KEGG microarray compendium consists of 42 datasets
#' investigating a total of 19 different human diseases as collected by Tarca
#' et al. (2012 and 2013).
#' 
#' The pre-defined TCGA RNA-seq compendium consists of datasets from The Cancer
#' Genome Atlas (TCGA, 2013) investigating a total of 34 different cancer
#' types.
#' 
#' User-defined data can also be loaded, given that datasets, preferably of
#' class \code{\linkS4class{SummarizedExperiment}}, have been saved as
#' \code{RDS} files.
#' 
#' @param edata Expression data compendium.  A character vector of length 1
#' that must be either \itemize{ \item 'geo2kegg': to load the GEO2KEGG
#' microarray compendium, \item 'tcga': to load the TCGA RNA-seq compendium, or
#' \item an absolute file path pointing to a directory, in which a user-defined
#' compendium has been saved in RDS files.  } See details.
#' @param nr.datasets Integer.  Number of datasets that should be loaded from
#' the compendium.  This is mainly for demonstration purposes.
#' @param cache Logical.  Should an already cached version used if available?
#' Defaults to \code{TRUE}.
#' @param ...  Additional arguments passed to the internal loading routines of
#' the GEO2KEGG and TCGA compendia.  This currently includes for loading of the
#' GEO2KEGG compendium \itemize{ \item \code{preproc}: logical. Should probe level
#' data automatically be summarized to gene level data? Defaults to \code{FALSE}.
#' \item \code{de.only}: logical. Include only datasets in which differentially
#' expressed genes have been found? Defaults to \code{FALSE}. \item \code{excl.metac}:
#' logical. Exclude datasets for which MetaCore rather than KEGG pathways have
#' been assigned as target pathways?  Defaults to \code{FALSE}.  } And for loading of
#' the TCGA compendium \itemize{ \item \code{mode}: character, determines how TCGA 
#' RNA-seq datasets are obtained. To obtain raw read counts from GSE62944 use 
#' either \code{'ehub'} (default, via ExperimentHub) or \code{'geo'} (direct
#' download from GEO, slow). Alternatively, use \code{'cTD'} to obtain normalized 
#' log2 TPM values from curatedTCGAData.  
#' \item \code{data.dir}: character. Absolute file path
#' indicating where processed RDS files for each dataset are written to.
#' Defaults to \code{NULL}, which will then write to
#' \code{rappdirs::user_data_dir("GSEABenchmarkeR")}.
#' \item \code{min.ctrls}: integer. Minimum number of controls, i.e. adjacent
#' normal samples, for a cancer type to be included. Defaults to 9. 
#' \item \code{paired}: Logical. Should the pairing of samples (tumor and adjacent
#' normal) be taken into account? Defaults to \code{TRUE}, which reduces the
#' data for each cancer type to patients for which both sample types 
#' (tumor and adjacent normal) are available. Use \code{FALSE} to obtain all 
#' samples in an unpaired manner. 
#' \item \code{min.cpm}: integer. Minimum
#' counts-per-million reads mapped.  See the edgeR vignette for details. The
#' default filter is to exclude genes with cpm < 2 in more than half of the
#' samples.  \item \code{with.clin.vars}: logical. Should clinical variables (>500) be
#' kept to allow for more advanced sample groupings in addition to the default
#' binary grouping (tumor vs. normal)? \item \code{map2entrez}: Should human gene symbols
#' be automatically mapped to Entrez Gene IDs? Defaults to \code{TRUE}.}
#' @return A \code{list} of datasets, typically of class
#' \code{\linkS4class{SummarizedExperiment}}. 
#'
#' Note that \code{loadEData("geo2kegg", preproc = FALSE)} (the default) 
#' returns the original microarray probe level data as a list of 
#' \code{\linkS4class{ExpressionSet}} objects. Use \code{preproc = TRUE} or
#' the \code{\link{maPreproc}} function to summarize the probe level
#' data to gene level data and to obtain a \code{list} of 
#' \code{\linkS4class{SummarizedExperiment}} objects.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\linkS4class{SummarizedExperiment}}, 
#' \code{\linkS4class{ExpressionSet}}, \code{\link{maPreproc}}
#' @references Tarca et al. (2012) Down-weighting overlapping genes improves
#' gene set analysis.  BMC Bioinformatics, 13:136.
#' 
#' Tarca et al. (2013) A comparison of gene set analysis methods in terms of
#' sensitivity, prioritization and specificity.  PLoS One, 8(11):e79217.
#' 
#' The Cancer Genome Atlas Research Network (2013) The Cancer Genome Atlas
#' Pan-Cancer analysis project.  Nat Genet, 45(10):1113-20.
#' 
#' Rahman et al. (2015) Alternative preprocessing of RNA-Sequencing data in The
#' Cancer Genome Atlas leads to improved analysis results.  Bioinformatics,
#' 31(22):3666-72.
#' @examples
#' 
#'     # (1) Loading the GEO2KEGG microarray compendium
#'     geo2kegg <- loadEData("geo2kegg", nr.datasets=2)
#' 
#'     # (2) Loading the TCGA RNA-seq compendium
#'     tcga <- loadEData("tcga", nr.datasets=2)
#' 
#'     # (3) reading user-defined expression data from file
#'     data.dir <- system.file("extdata/myEData", package="GSEABenchmarkeR")
#'     edat <- loadEData(data.dir)
#' 
#' @export loadEData
loadEData <- function(edata, nr.datasets=NULL, cache=TRUE, ...)
{
    stopifnot(is.character(edata))
    if(edata == "geo2kegg") exp.list <- .loadGEO2KEGG(nr.datasets, cache, ...)
    else if(edata == "tcga") exp.list <- .loadTCGA(nr.datasets, cache, ...)
    else if(file.exists(edata)) 
        exp.list <- .loadEDataFromFile(edata, nr.datasets, cache)
    else stop("edata must be \'geo2kegg\', \'tcga\' or an absolute file path")
    return(exp.list)
}

#
# loading expression data from file
#
.loadEDataFromFile <- function(edata, nr.datasets=NULL, cache=TRUE)
{
    data.ids <- sub(".rds$", "", list.files(edata, pattern="*.rds$"))
    if(!is.null(nr.datasets)) 
    {
        nr.datasets <- min(nr.datasets, length(data.ids))
        data.ids <- data.ids[seq_len(nr.datasets)]
    }

    # should a cached version be used?
    if(cache)
    {  
        eid <- basename(edata)
        exp.list <- .getResourceFromCache(rname=eid, update.value=NA)
        if(!is.null(exp.list)) return(exp.list[intersect(data.ids, names(exp.list))])
    }
         
    if(interactive())
        pb <- txtProgressBar(0, length(data.ids), style=3)
    exp.list <- lapply(data.ids,
        function(id)
        {
            dfile <- file.path(edata, paste(id, "rds", sep="."))
            se <- readRDS(dfile)
            metadata(se)$dataId <- id 
            if(interactive()) setTxtProgressBar(pb, match(id, data.ids))
            return(se)
        })
    names(exp.list) <- data.ids
    if(interactive()) close(pb)
    return(exp.list)
}

#' Caching of a resource
#' 
#' Convenience function to flexibly save and restore an already processed 
#' expression data compendium via caching.
#' 
#' @param res Resource. An arbitrary R object. 
#' @param rname Character. Resource name. 
#' @param ucdir Character. User cache directory. Defaults to 'GSEABenchmarkeR',
#' which will accordingly use \code{rappdirs::user_cache_dir("GSEABenchmarkeR")}.
#' @return None. Stores the object in the cache by side effect.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\link{loadEData}}, \code{\link{user_cache_dir}}, 
#' \code{\linkS4class{BiocFileCache}}
#' @examples
#' 
#'      # load user-defined expression compendium
#'      data.dir <- system.file("extdata/myEData", package = "GSEABenchmarkeR")
#'      edat <- loadEData(data.dir)
#'
#'      # do some processing of the compendium
#'      edat <- lapply(edat, function(d) d[1:50,])
#'
#'      # cache it ...
#'      cacheResource(edat, "myEData")
#'
#'      # ... and restore it at a later time
#'      edat <- loadEData(data.dir, cache = TRUE)    
#' 
#' @export
cacheResource <- function(res, rname, ucdir="GSEABenchmarkeR")
{
    # are we running interactively?
    cache.dir <- ifelse(interactive(), 
                    rappdirs::user_cache_dir(ucdir),
                    tempdir())

    if(!file.exists(cache.dir)) dir.create(cache.dir)
    bfc <- BiocFileCache::BiocFileCache(cache.dir)
    
    # replace existing version if necessary 
    qgsc <-  BiocFileCache::bfcquery(bfc, rname)
    if(BiocFileCache::bfccount(qgsc)) BiocFileCache::bfcremove(bfc, qgsc$rid) 
    
    cache.file <- BiocFileCache::bfcnew(bfc, rname)
    saveRDS(res, file=cache.file)
}

.getResourceFromCache <- function(rname, 
    update.value=7, update.unit="days", ucdir="GSEABenchmarkeR")
{
    # are we running interactively?
    cache.dir <- ifelse(interactive(), 
                    rappdirs::user_cache_dir(ucdir),
                    tempdir())

    if(!file.exists(cache.dir)) return(NULL) 

    bfc <- BiocFileCache::BiocFileCache(cache.dir)
    qgsc <-  BiocFileCache::bfcquery(bfc, rname)

    # is there already a cached version?
    res <- NULL
    if(BiocFileCache::bfccount(qgsc))
    {
        # is the cached version outdated?
        dt <- difftime(Sys.time(), qgsc$create_time, units=update.unit)   
        if(is.na(update.value) || dt < update.value)
        {
            if(interactive())
                message(paste("Using cached version from", qgsc$create_time))
            res <- readRDS(BiocFileCache::bfcrpath(bfc, rname))
        }
    }
    return(res)   
}


#
# write differential expression per dataset
#
#' @rdname runDE
#' @export
writeDE <- function(exp.list, out.dir=NULL)
{
    if(is.null(out.dir))
    {
        out.dir <- rappdirs::user_data_dir("GSEABenchmarkeR")
        out.dir <- file.path(out.dir, "de")
    }
    if(!file.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
    .f <- function(xx)
    {
        gt.file <- paste0(metadata(xx)$dataId, ".txt")
        gt.file <- file.path(out.dir, gt.file)
        gt <- EnrichmentBrowser:::.getGeneAnno(names(xx), "hsa")
        gt <- cbind(gt, rowData(xx, use.names=TRUE))
        gt <- EnrichmentBrowser:::.sortGeneTable(gt)  
        EZ.COL <- EnrichmentBrowser::configEBrowser("EZ.COL")
        ind <- gt[,EZ.COL]      
        write.table(gt, file=gt.file, row.names=FALSE, quote=FALSE, sep="\t")
    }
    invisible(lapply(exp.list, .f))
}


# reads a mapping from datasetID (eg GSExxx) to diseaseID (eg ALZ)
# as eg required for the Tarca microarray compendium


#' Read a mapping between dataset ID and disease code
#' 
#' When assessing enrichment analysis results for phenotype relevance, it is
#' assumed that each analyzed dataset investigates a certain phenotype such as
#' a disease.  This function reads a mapping between dataset IDs and assigned
#' disease codes.
#' 
#' 
#' @param map.file Character. The path to the mapping file.
#' @return A named character vector where each element of the vector is a
#' disease code and the names are the dataset IDs.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{evalRelevance} for evaluating phenotype relevance of gene set
#' rankings.
#' @examples
#' 
#'     data.dir <- system.file("extdata", package="GSEABenchmarkeR")
#'     d2d.file <- file.path(data.dir, "malacards", "GseId2Disease.txt")
#'     d2d.map <- readDataId2diseaseCodeMap(d2d.file)
#' 
#' @export readDataId2diseaseCodeMap
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


#
# READ RESULTS (RUNTIME, RANKINGS, TYPE I ERROR RATE)
#

#' Reading results of enrichment analysis
#' 
#' These functions read results obtained from the application of
#' enrichment methods to multiple datasets for subsequent assessment.
#' 
#' 
#' @param data.dir Character.  The data directory where results have been saved
#' to.
#' @param data.ids A character vector of dataset IDs.
#' @param methods Methods for enrichment analysis.  A character vector with
#' method names typically chosen from \code{\link{sbeaMethods}}
#' and \code{\link{nbeaMethods}}, or user-defined functions
#' implementing methods for enrichment analysis.
#' @param type Character. Type of the result. Should be one out of 'runtime',
#' 'ranking', or 'typeI'. 
#' @return A result list with an entry for each method applied.  Each entry
#' stores corresponding runtimes (\code{type="runtime"} ), gene set
#' rankings (\code{type="ranking"}), or type I error rates (\code{type="typeI"})
#' as obtained from applying the respective method to the given datasets.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{runEA} to apply enrichment methods to multiple datasets.
#' @examples
#' 
#'     # simulated setup: 
#'     # 1 methods & 1 datasets
#'     methods <- paste0("m", 1:2)
#'     data.ids <- paste0("d", 1:2)
#' 
#'     # result directory
#'     res.dir <- tempdir()
#'     sdirs <- file.path(res.dir, methods)
#'     for(d in sdirs) dir.create(d)
#'     
#'     # store runtime & rankings 
#'     for(m in 1:2)
#'     {
#'         rt <- runif(5, min=m, max=m+1)
#'         for(d in 1:2)
#'         {
#'             # runtime
#'             out.file <- paste(data.ids[d], "txt", sep=".")
#'             out.file <- file.path(sdirs[m], out.file)
#'             cat(rt[d], file=out.file) 
#' 
#'             # ranking
#'             out.file <- sub("txt$", "rds", out.file)
#'             r <- EnrichmentBrowser::makeExampleData("ea.res") 
#'             r <- EnrichmentBrowser::gsRanking(r, signif.only=FALSE)
#'             saveRDS(r, file=out.file)   
#'         }
#'     }
#' 
#'     # reading runtime & rankings
#'     rts <- readResults(res.dir, data.ids, methods, type="runtime")
#'     rkgs <- readResults(res.dir, data.ids, methods, type="ranking")
#'
#' @export readResults
readResults <- function(data.dir, data.ids, methods, 
    type=c("runtime", "ranking", "typeI"))
{
    type <- match.arg(type)
    res <- lapply(methods, .readResultsOfMethod, 
                    data.dir=data.dir, data.ids=data.ids, type=type)
    names(res) <- methods
    return(res)
}

.readResultsOfMethod <- function(m, data.dir, data.ids, type)
{
    ext <- ifelse(type %in% c("ranking", "typeI"), ".rds", ".txt")
    
    mdir <- file.path(data.dir, m)
    mfiles <- paste0(data.ids, ext)
    mfiles <- file.path(mdir, mfiles) 
    if(type == "runtime")
    { 
        res <- vapply(mfiles, .readResultFile, 
                        numeric(1), type="runtime")
        ind <- !is.na(res)
    }
    else 
    {
        res <- lapply(mfiles, .readResultFile, type=type)
        ind <- !vapply(res, is.null, logical(1))
    }
    ext <- paste0(ext, "$")
    names(res) <- sub(ext, "", basename(mfiles))
    res <- res[ind] 
    if(type == "typeI") res <- simplify2array(res, higher=FALSE)
    return(res)
}

.readResultFile <- function(res.file, type)
{
    ext <- ifelse(type %in% c("ranking", "typeI"), ".rds", ".txt")
    m <- basename(dirname(res.file))
    d <- sub(ext, "", basename(res.file))
    id <- paste(m, d, type, sep="_")   
    
    if(file.exists(res.file))
    { 
        if(ext == ".rds") res <- readRDS(res.file)
        else res <- scan(res.file, quiet=TRUE)
    }
    else
    {
        rtype <- switch(type,
                        runtime = "Runtime",
                        ranking = "Ranking",
                        typeI = "Type I error rate")
    
        ext <- paste0(ext, "$")
        wmsg <- paste(rtype, "of", m, "for", d, "not found")
        warning(wmsg)
        res <- if(type == "rt") NA else NULL
    }
    return(res)
}

