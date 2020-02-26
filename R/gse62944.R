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
.loadTCGA <- function(nr.datasets=NULL, cache=TRUE, mode=c("ehub", "geo"), 
    data.dir=NULL, min.ctrls=9, min.cpm=2, with.clin.vars=FALSE, paired=TRUE)
{
    message("Loading TCGA data compendium ...")
    # should a cached version be used?
    if(cache)
    {  
        # ignore all other arguments 
        el <- .getResourceFromCache(rname="tcga", update.value=NA)
        if(!is.null(el)) return(el[seq_len(min(nr.datasets, length(el)))])
    }      
 
    # otherwise download it
    mode <- match.arg(mode)
    
    if(mode == "ehub") dat <- .retrieveGSE62944FromEHub(with.clin.vars)
    else dat <- .downloadGSE62944(data.dir, with.clin.vars)
    el <- .splitByCancerType(dat$tum, dat$norm, 
        nr.datasets, min.ctrls=min.ctrls, min.cpm=min.cpm, paired=paired)
    
    if(interactive()) cacheResource(el, "tcga")
    return(el)
}

#
# preferred way of getting GSE62944 (fast)
#
.retrieveGSE62944FromEHub <- function(with.clin.vars=FALSE)
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
.downloadGSE62944 <- function(data.dir, with.clin.vars=FALSE)
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
        EnrichmentBrowser::isAvailable("GEOquery", type = "software")
        GEOquery::getGEOSuppFiles("GSE62944", baseDir = dirname(data.dir))
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

.splitByCancerType <- function(tum, norm, 
    nr.datasets, min.ctrls=9, min.cpm=2, paired=FALSE)
{
    # map IDs hgnc -> entrez
    suppressMessages({
        suppressPackageStartupMessages(
            tum <- EnrichmentBrowser::idMap(tum, 
                    org="hsa", from="SYMBOL", to="ENTREZID")
        )
        norm <- EnrichmentBrowser::idMap(norm, 
                    org="hsa", from="SYMBOL", to="ENTREZID")
    })

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
            se$GROUP <- c( rep(1, ncol(cases)), rep(0, ncol(ctrls)) )
            if(paired) se$BLOCK <- rep(isect, 2)
            metadata(se)$dataId <- ct 
            metadata(se)$dataType <- "rseq"

            # preprocessing
            # rm genes with all 0 and less than min.reads
            rs <- rowSums(edgeR::cpm(assay(se)) > min.cpm)
            keep <-  rs >= ncol(se) / 2
            se[keep,]
        })
    names(exp.list) <- rel.cts
    
    nr.samples <- vapply(exp.list, ncol, integer(1)) 
    exp.list[nr.samples >= min.ctrls * 2]
}
