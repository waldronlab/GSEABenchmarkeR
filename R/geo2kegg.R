############################################################
# 
# author: Ludwig Geistlinger
# date: 2016-07-29 16:49:12
# 
# descr: basic functionality for evaluating KEGGdzGEO 
# (i.e. GEO datasets, 
#   with defined presumably perturbed KEGG pathways)
# 
############################################################

#
# loading geo2kegg benchmark compendium
#
.loadGEO2KEGG <- function(nr.datasets=NULL, 
    cache=TRUE, preproc=FALSE, de.only=FALSE, excl.metac=FALSE)
{
    data.ids <- .getGeo2KeggIds(de.only, excl.metac)
    if(!is.null(nr.datasets)) 
    {
        nr.datasets <- min(nr.datasets, length(data.ids))  
        data.ids <- data.ids[seq_len(nr.datasets)]
    }
    message("Loading GEO2KEGG data compendium ...")
    
    # should a cached version be used?
    if(cache)
    {   
        # using a cached version ignores all other arguments
        # such as nr.datasets, de.only, ... 
        el <- .getResourceFromCache(rname="geo2kegg", update.value=NA)
        if(!is.null(el)) return(el)
    }       

    # if not, create from package
    geoIds1 <- data(package="KEGGdzPathwaysGEO")$results[,"Item"]
    if(interactive()) 
        pb <- txtProgressBar(0, length(data.ids), style=3)

    .EDATA <- environment()
    el <- lapply(data.ids, 
        function(id) 
        {
            # load
            add <- ifelse(id %in% geoIds1, "d", "andMetacoreD")
            pkg <- paste0("KEGG", add, "zPathwaysGEO")
            data(list=id, package=pkg, envir=.EDATA)
            if(interactive()) setTxtProgressBar(pb, match(id, data.ids))
            return(.EDATA[[id]])
        })
    if(interactive()) close(pb)
    names(el) <- data.ids

    if(preproc)
    {
        el <- maPreproc(el)
        cacheResource(el, "geo2kegg")
    }
    return(el)
}



#' Preprocessing of microarray expression data
#' 
#' This function prepares datasets of the GEO2KEGG microarray compendium for
#' further analysis. This includes summarization of probe level expression to
#' gene level expression as well as annotation of required colData slots for
#' sample grouping.
#' 
#' 
#' @param exp.list Experiment list.  A \code{list} of datasets, each being of
#' class \code{\linkS4class{ExpressionSet}}.
#' @param parallel Parallel computation mode.  An instance of class
#' \code{\linkS4class{BiocParallelParam}}.  See the vignette of the
#' \code{BiocParallel} package for switching between serial, multi-core, and
#' grid execution.  Defaults to \code{NULL}, which then uses the first element
#' of \code{BiocParallel::registered()} for execution.  If not changed by the
#' user, this accordingly defaults to multi-core execution on the local host.
#' @return A \code{list} of datasets, each being of class
#' \code{\linkS4class{SummarizedExperiment}}.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\link{loadEData}} to load a specified expression data 
#' compendium.
#' @examples
#' 
#'     # reading user-defined expression data from file
#'     geo2kegg <- loadEData("geo2kegg", nr.datasets=3)
#' 
#'     # only considering the first 100 probes for demonstration
#'     geo2kegg <- lapply(geo2kegg, function(d) d[1:100,])   
#'
#'     # preprocessing two datasets
#'     geo2kegg <- maPreproc(geo2kegg[2:3])
#' 
#' @export maPreproc
maPreproc <- function(exp.list, parallel=NULL)
{
    message("Summarizing probe level expression ...")
    anno.pkgs <- vapply(exp.list, annotation, character(1))

    anno.pkgs <- unique(anno.pkgs)
    anno.pkgs <- paste(anno.pkgs, "db", sep=".")

    # load the annotation packages
    suppressPackageStartupMessages(
        for(p in anno.pkgs) EnrichmentBrowser::isAvailable(p)
    )

    # retrieve the mappings
    suppressMessages(
        p2g.maps <- sapply(anno.pkgs, 
            function(anno.pkg)
            {
                anno.pkg <- get(anno.pkg) 
                m <- AnnotationDbi::mapIds(anno.pkg, 
                        keys=AnnotationDbi::keys(anno.pkg), 
                        keytype="PROBEID", column="ENTREZID")
                return(m)
            }, simplify=FALSE) 
    )

    # helper function doing the actual job per expression set
    GRPCOL <- EnrichmentBrowser::configEBrowser("GRP.COL")
    BLKCOL <- EnrichmentBrowser::configEBrowser("BLK.COL")
    .pre <- function(eset)
    {
        # probe 2 gene & annotation
        id <- experimentData(eset)@name
        anno <- paste0(annotation(eset), ".db")
        fData(eset)$ENTREZID <- p2g.maps[[anno]][rownames(eset)] 
        se <- EnrichmentBrowser::probe2gene(eset)
        colnames(colData(se))[2] <- GRPCOL 
        if(ncol(colData(se)) > 2) colnames(colData(se))[3] <- BLKCOL 
        se[[GRPCOL]] <- ifelse(se[[GRPCOL]] == "d", 1, 0)
        metadata(se)$dataId <- id 
        metadata(se)$dataType <- "ma"
        metadata(se)$annotation <- "hsa"
        return(se) 
    }

    exp.list <- .iter(exp.list, .pre, parallel=parallel)
    return(exp.list)
}

#
# get GEO IDs of geo2kegg benchmark compendium
# 
.getGeo2KeggIds <- function(de.only=FALSE, excl.metac=FALSE)
{
    geoIds1 <- data(package="KEGGdzPathwaysGEO")$results[,"Item"]
    geoIds2 <- data(package="KEGGandMetacoreDzPathwaysGEO")$results[,"Item"]
    geoIds <- c(geoIds1, geoIds2)

    # metacore geo ids    
    if(excl.metac)
    {
        tgs <- function(d) notes(experimentData(d))$targetGeneSets
        g2k <- .iter(tgs, geoIds)

        metacTgs <- g2k[!grepl("^[0-9]+$", g2k)]
        metacGeoIds <- names(metacTgs) 
        geoIds <- geoIds[-match(metacGeoIds, geoIds)]
    }

    # de geo ids
    if(de.only)
    {
        #de <- iter.geos(fract.de, geo.ids=geo.ids, alpha=0.05, beta=1)
        #de <- t(de)
        #non.de.geos <- rownames(de)[de[,1]==0 | de[,2]==0] 
        nonDeGeos <- paste0("GSE", c("20153", "20291", 
            "6956AA", "6956C", "781", "8762", "16759",
            "19420", "20164", "22780", "30153", "42057"))
        geoIds <- geoIds[!(geoIds %in% nonDeGeos)]
    }

    return(geoIds)
}

