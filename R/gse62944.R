############################################################
# 
# author: Ludwig Geistlinger
# date: 2016-06-07 13:21:22
# 
# descr: get GSE62944 (TCGA RNA-seq read counts) ready for
#           GSEA benchmark 
# 
############################################################

#
# load TCGA compendium
#
.loadTCGA <- function(  nr.datasets = NULL, 
                        cache = TRUE, 
                        mode = c("ehub", "geo", "cTD"),
                        data.dir = NULL,
                        min.ctrls = 9,
                        min.cpm = 2,
                        with.clin.vars = FALSE,
                        paired = TRUE,
                        map2entrez = TRUE)
{
    message("Loading TCGA data compendium ...")
    # should a cached version be used?
    if(cache)
    {  
        # ignore all other arguments 
        el <- .getResourceFromCache(rname="tcga", update.value = NA)
        if(!is.null(el)) return(el[seq_len(min(nr.datasets, length(el)))])
    }      
 
    # otherwise download it
    mode <- match.arg(mode)
    
    if(mode == "cTD") 
        el <- .retrieveFromCTD(nr.datasets, min.ctrls, min.cpm, 
                                with.clin.vars, paired, map2entrez)
    else
    {
        if(mode == "ehub") dat <- .retrieveGSE62944FromEHub(with.clin.vars)
        else dat <- .downloadGSE62944(data.dir, with.clin.vars)
        el <- .splitByCancerType(dat$tum, dat$norm, nr.datasets, 
                                    min.ctrls, min.cpm, paired, map2entrez)
    } 
    if(interactive()) cacheResource(el, "tcga")
    return(el)
}

#
# preferred way of getting GSE62944 (fast)
#
.retrieveGSE62944FromEHub <- function(with.clin.vars = FALSE)
{
    suppressMessages({
        hub <- ExperimentHub::ExperimentHub()
        gse <- AnnotationHub::query(hub, "GSE62944")
        tum <- gse[["EH1043"]]
   }) 

    ind <- grep("CancerType", colnames(colData(tum)))
    colnames(colData(tum))[ind] <- "type"

    # should clinical variables be dropped?
    if(!with.clin.vars)
    {
        colData(tum) <- colData(tum)[,c(1,ind)]
        colData(tum)[,1] <- colnames(tum)
        colnames(colData(tum))[1] <- "sample"
    }

    suppressMessages( norm <- gse[["EH1044"]] )
    names(assays(norm)) <- names(assays(tum)) <- "exprs" 
    
    list(tum = tum, norm = norm)
}

#
# alternative way of getting GSE62944 in case ExperimentHub fails
#
.downloadGSE62944 <- function(data.dir, with.clin.vars = FALSE)
{
    # output data directory
    if(is.null(data.dir))
        data.dir <- rappdirs::user_data_dir("GSEABenchmarkeR")

    data.dir <- file.path(data.dir, "GSE62944")
    tum.file <- "GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FeatureCounts.txt.gz"
    norm.file <- "GSM1697009_06_01_15_TCGA_24.normal_Rsubread_FeatureCounts.txt.gz"
    tum.cl.file <- "GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt.gz"
    norm.cl.file <- "GSE62944_06_01_15_TCGA_24_Normal_CancerType_Samples.txt.gz"
    clin.var.file <- "GSE62944_06_01_15_TCGA_24_548_Clinical_Variables_9264_Samples.txt.gz"
    
    rel.files <- c(tum.file, norm.file, tum.cl.file, norm.cl.file, clin.var.file)
    rel.files <- file.path(data.dir, rel.files)

    if(!all(file.exists(rel.files)))
    {
        message("Downloading GSE62944 from GEO ...")
        getGEOSuppFiles <- NULL
        EnrichmentBrowser::isAvailable("GEOquery", type = "software")
        getGEOSuppFiles("GSE62944", baseDir = dirname(data.dir))
        message(paste("Data files are stored under:", data.dir))
        tar.file <- file.path(data.dir, "GSE62944_RAW.tar") 
        untar(tar.file, basename(rel.files[1:2]), exdir = data.dir)
    }

    message("Reading feature counts ...")
    tum.cont <- read.delim(rel.files[1], row.names=1L, as.is=TRUE)
    norm.cont <- read.delim(rel.files[2], row.names=1L, as.is=TRUE)
    colnames(norm.cont) <- gsub("\\.", "-", colnames(norm.cont))
    colnames(tum.cont) <- gsub("\\.", "-", colnames(tum.cont))

    # cancer types
    message("Reading sample annotation (cancer/normal) ...")
    tum.cl <- read.delim(rel.files[3], header=FALSE, as.is=TRUE)
    norm.cl <- read.delim(rel.files[4], header=FALSE, as.is=TRUE)
    colnames(tum.cl) <- colnames(norm.cl) <- c("sample", "type")

    # clinical variables
    if(with.clin.vars)
    {
        clin.var <- read.delim(rel.files[5], as.is=TRUE) 
        n <- clin.var[,1] 
        clin.var <- clin.var[,-c(1,2,3)]
        colnames(clin.var) <- gsub("\\.", "-", colnames(clin.var))
        clin.var <- t(clin.var)
        colnames(clin.var) <- n
        tum.cl <- cbind(tum.cl, clin.var)
    }

    # clean up
    message("Cleaning up ...")
    rm.files <- list.files(data.dir, full.names=TRUE) 
    rm.files <- rm.files[!(rm.files %in% rel.files)]
    file.remove(rm.files)

    # create SummarizedExperiment for tumor and normal 
    message("Creating a SummarizedExperiment for tumor and normal samples ...")
    
    tum.cont <- tum.cont[tum.cl$sample]
    tum <- SummarizedExperiment(assays = list(exprs = as.matrix(tum.cont)),
                                colData = DataFrame(tum.cl))    
 
    norm.cont <- norm.cont[norm.cl$sample]
    norm <- SummarizedExperiment(assays = list(exprs = as.matrix(norm.cont)),
                                 colData = DataFrame(norm.cl))

    list(tum = tum, norm = norm)
}

.retrieveFromCTD <- function(nr.datasets = NULL, min.ctrls = 9, min.cpm = 2, 
    with.clin.vars = FALSE, paired = TRUE, map2entrez = TRUE)
{
    message("Obtaining datasets through curatedTCGAData ...")
    curatedTCGAData <- experiments <- NULL
    EnrichmentBrowser::isAvailable("curatedTCGAData", type = "data")
  
    rnacomp <- .getResourceFromCache(rname = "cTD_rnacomp")
    if(is.null(rnacomp))
    {
        rnacomp <- curatedTCGAData("*", "RNASeq2*", FALSE)
        cacheResource(rnacomp, "cTD_rnacomp")
    }

    # extract the list of SEs
    exp.list <- experiments(rnacomp)
    .first <- function(n) unlist(strsplit(n, "_"))[1]    
    names(exp.list) <- vapply(names(exp.list), .first, character(1))

    message("Cancer types with tumor samples:")
    message(paste(names(exp.list), collapse=", "))

    # extract samples types (01 = tumor, 11 = adj. normal)
    stypes <- lapply(exp.list, .getSampleType) 
    stab <- lapply(stypes, table)

    .hasAControl <- function(x) "11" %in% names(x)
    ind <- vapply(stab, .hasAControl, logical(1))
    message("Cancer types with adj. normal samples:")
    message(paste(names(exp.list[ind]), collapse=", "))

    # has sufficient controls?
    .hasControls <- function(x) 
        all(c("01", "11") %in% names(x)) && x["11"] >= min.ctrls
    ind <- vapply(stab, .hasControls, logical(1))
    exp.list <- exp.list[ind]
    
    # specified number of datsets?
    if(!is.null(nr.datasets)) 
        exp.list <- exp.list[seq_len(min(nr.datasets, length(exp.list)))]

    message("Cancer types with sufficient tumor and adj. normal samples:")
    message(paste(names(exp.list), collapse=", "))
    message("Creating a SummarizedExperiment for each of them ...")
    
    exp.list <- lapply(exp.list, .groupSE, report = !paired)
    for(i in seq_along(exp.list)) 
        metadata(exp.list[[i]])$dataId <- names(exp.list)[i]   
    
    # paired?
    if(paired)
    { 
        exp.list <- lapply(exp.list, .pairSE)
        nr.samples <- vapply(exp.list, ncol, integer(1)) 
        exp.list <- exp.list[nr.samples >= min.ctrls * 2]
    }

    # exclude genes not sufficiently expressed
    exp.list <- lapply(exp.list, .filterByCPM, min.cpm = min.cpm, is.raw = FALSE)

    # map gene symbols to Entrez Gene IDs?
    if(map2entrez) 
        suppressMessages(exp.list <- lapply(exp.list, EnrichmentBrowser::idMap, 
                                org = "hsa", from = "SYMBOL", to = "ENTREZID"))
 
    # annotate clinical variables?
    if(with.clin.vars)
        exp.list <- lapply(exp.list, .extractColData, cdat = colData(rnacomp))

    return(exp.list) 
}

.filterByCPM <- function(se, min.cpm, is.raw = TRUE)
{
    mat <- assay(se)
    if(is.raw) mat <- edgeR::cpm(mat)
    rs <- rowSums(mat > min.cpm)
    keep <- rs >= ncol(se) / 2
    se <- se[keep,]
    if(!is.raw)
    { 
        assay(se) <- log2(assay(se) + 2)
        metadata(se)$dataType <- "ma"
    }
    return(se)
}

.getSampleType <- function(se) substring(colnames(se), 14, 15)

.groupSE <- function(se, report = TRUE)
{
    GRP.COL <- EnrichmentBrowser::configEBrowser("GRP.COL")
    stypes <- .getSampleType(se)
    ind <- stypes %in% c("01", "11")
    se <- se[,ind]
    stypes <- stypes[ind]
    se[[GRP.COL]] <- ifelse(stypes == "01", 1, 0) 
    metadata(se)$dataType <- "rseq"
    
    if(report) 
        message(metadata(se)$dataId, 
                " tumor: ", sum(se[[GRP.COL]] == 1), 
                " adj.normal: ", sum(se[[GRP.COL]] == 0))

    return(se)
}

.pairSE <- function(se)
{
    GRP.COL <- EnrichmentBrowser::configEBrowser("GRP.COL")
    BLK.COL <- EnrichmentBrowser::configEBrowser("BLK.COL")

    patient.ids <- substring(colnames(se), 1, 12)
    ctrls <- patient.ids[se[[GRP.COL]] == 0] 
    cases <- patient.ids[se[[GRP.COL]] == 1] 
    isect <- intersect(ctrls, cases)    
    ind <- patient.ids %in% isect
    se <- se[,ind]
    se[[BLK.COL]] <- patient.ids[ind] 

    message(metadata(se)$dataId, 
            " tumor: ", sum(se[[GRP.COL]] == 1), 
            " adj.normal: ", sum(se[[GRP.COL]] == 0))
    return(se)
}

.extractColData <- function(se, cdat)
{
    snames <- substring(colnames(se), 1, 12)    
    colData(se) <- cbind(colData(se), cdat[snames,])
    return(se)    
}

.splitByCancerType <- function(tum, norm, nr.datasets, 
    min.ctrls = 9, min.cpm = 2, paired = FALSE, map2entrez = TRUE)
{
    GRP.COL <- EnrichmentBrowser::configEBrowser("GRP.COL")
    BLK.COL <- EnrichmentBrowser::configEBrowser("BLK.COL")

    if(map2entrez)
    {
        # map IDs hgnc -> entrez
        suppressMessages({
            tum <- EnrichmentBrowser::idMap(tum, org = "hsa", 
                                            from = "SYMBOL", to = "ENTREZID")
            norm <- EnrichmentBrowser::idMap(norm, org = "hsa", 
                                                from = "SYMBOL", to = "ENTREZID")
        })
    }

    # creating SEs
    tum.cts <- table(colData(tum)[,"type"])
    message("Cancer types with tumor samples:")
    message(paste(names(tum.cts), collapse=", "))
    norm.cts <- table(colData(norm)[,"type"])
    message("Cancer types with adj. normal samples:")
    message(paste(names(norm.cts), collapse=", "))

    norm.cts.min.ctrls <- norm.cts[norm.cts >= min.ctrls]
    rel.cts <- intersect(names(norm.cts.min.ctrls), names(tum.cts))

    if(!is.null(nr.datasets)) 
    {
        nr.datasets <- min(nr.datasets, length(rel.cts)) 
        rel.cts <- rel.cts[seq_len(nr.datasets)] 
    }

    message("Cancer types with sufficient tumor and adj. normal samples:")
    message(paste(rel.cts, collapse=", "))
    message("Creating a SummarizedExperiment for each of them ...")

    diff.cols <- setdiff(colnames(colData(tum)), colnames(colData(norm)))
    if(length(diff.cols))
    {
        mat <- matrix(NA, nrow = ncol(norm), ncol = length(diff.cols))
        colnames(mat) <- diff.cols
        colData(norm) <- cbind(colData(norm), DataFrame(mat))
        colData(norm) <- colData(norm)[,colnames(colData(tum))]
    }

    exp.list <- lapply(rel.cts,
        function(ct) 
        {

            ctrls <- norm[, norm$type == ct]
            cases <- tum[, tum$type == ct]

            if(paired)
            {
                pre.norm <- substring(colnames(ctrls), 1, 12)
                pre.tum <- substring(colnames(cases), 1, 12)
                isect <- intersect(pre.tum, pre.norm)    
                cases <- cases[,match(isect, pre.tum)]
                ctrls <- ctrls[,match(isect, pre.norm)]    
            }

            message(paste(ct, "tumor:", ncol(cases), 
                            "adj.normal:", ncol(ctrls)))
            
            se <- cbind(cases, ctrls) 
            ind <- unique(names(metadata(se)))
            metadata(se) <- metadata(se)[ind]
            se[[GRP.COL]] <- c( rep(1, ncol(cases)), rep(0, ncol(ctrls)) )
            if(paired) se[[BLK.COL]] <- rep(isect, 2)
            metadata(se)$dataId <- ct 
            metadata(se)$dataType <- "rseq"

            # preprocessing
            # rm genes with all 0 and less than min.reads
            .filterByCPM(se, min.cpm)
        })
    names(exp.list) <- rel.cts
    
    nr.samples <- vapply(exp.list, ncol, integer(1)) 
    exp.list[nr.samples >= min.ctrls * 2]
}
