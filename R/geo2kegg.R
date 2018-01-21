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
.loadGEO2KEGG <- function(reload=FALSE, nr.datasets=NULL,
                            preproc=FALSE, de.only=FALSE, excl.metac=FALSE)
{
    data.ids <- .getGeo2KeggIds(de.only, excl.metac)
    if(!is.null(nr.datasets)) data.ids <- data.ids[seq_len(nr.datasets)]
    isLoaded <- all(vapply(data.ids, exists, logical(1), envir=.EDATA))
    if(!isLoaded || reload)
    {
        message("Loading GEO2KEGG data compendium ...")
        
        # already on disk?
        data.dir <- rappdirs::user_data_dir("GSEABenchmarkeR")
        gse.dir <- file.path(data.dir, "GEO2KEGG")
        if(file.exists(gse.dir))
        {        
            message("Found GEO2KEGG compendium on disk. Loading from file ...")
            eids <- .loadEDataFromFile(gse.dir, nr.datasets=nr.datasets)
            if(all(data.ids %in% eids)) return(data.ids)
        }
        
        # if not, create from package
        geoIds1 <- data(package="KEGGdzPathwaysGEO")$results[,"Item"]
        if(interactive()) 
            pb <- txtProgressBar(1, length(data.ids), style=3, width=length(data.ids))
        for(did in data.ids) 
        {
            # load
            add <- ifelse(did %in% geoIds1, "d", "andMetacoreD")
            pkg <- paste0("KEGG", add, "zPathwaysGEO")
            data(list=did, package=pkg, envir=.EDATA)
            if(interactive()) setTxtProgressBar(pb, match(did, data.ids))
        }    
        if(interactive()) close(pb)
        if(preproc) maPreprocApply(data.ids)
    }
    return(data.ids)
}

maPreprocApply <- function(data.ids, parallel=NULL)
{
    message("Summarizing probe level expression ...")
    anno.pkgs <- vapply(data.ids, 
                    function(i) annotation(getDataset(i)), 
                    character(1))

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

    GRPCOL <- EnrichmentBrowser::config.ebrowser("GRP.COL")
    BLKCOL <- EnrichmentBrowser::config.ebrowser("BLK.COL")
    .pre <- function(eset)
    {
        # probe 2 gene & annotation
        id <- experimentData(eset)@name
        anno <- paste0(annotation(eset), ".db")
        fData(eset)$ENTREZID <- p2g.maps[[anno]][rownames(eset)] 
        eset <- EnrichmentBrowser::probe.2.gene.eset(eset)
        colnames(colData(eset))[2] <- GRPCOL 
        if(ncol(colData(eset)) > 2) colnames(colData(eset))[3] <- BLKCOL 
        eset[[GRPCOL]] <- ifelse(eset[[GRPCOL]] == "d", 1, 0)
        metadata(eset)$dataId <- id 
        metadata(eset)$dataType <- "ma"
        metadata(eset)$annotation <- "hsa"
        return(eset) 
    }
    esets <- .iter(.pre, data.ids, parallel)
    for (d in data.ids) .EDATA[[d]] <- esets[[d]] 
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

