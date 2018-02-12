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
.loadTCGA <- function(nr.datasets=NULL,
                        mode=c("ehub", "geo"), data.dir=NULL,
                        min.ctrls=9, min.cpm=2, with.clin.vars=FALSE)
{
    # output data directory
    if(is.null(data.dir))
        data.dir <- rappdirs::user_data_dir("GSEABenchmarkeR")
   
    # check whether GSE62944 is already there
    gse.dir <- file.path(data.dir, "GSE62944")
    if(file.exists(gse.dir))
    {
        message("Found GSE62944 on disk. Loading from file ...")
        exp.list <- .loadEDataFromFile(gse.dir, nr.datasets=nr.datasets)
    }
    # otherwise download it
    else
    {
        mode <- match.arg(mode)
        if(mode == "ehub") dat <- .retrieveGSE62944FromEHub(with.clin.vars)
        else dat <- .downloadGSE62944(data.dir, with.clin.vars)
        exp.list <- .splitByCancerType(dat$tum, dat$norm, nr.datasets,
            data.dir=gse.dir, min.ctrls=min.ctrls, min.cpm=min.cpm)
    }
     
    return(exp.list)
}

#
# preferred way of getting GSE62944 (fast)
#
.retrieveGSE62944FromEHub <- function(with.clin.vars=FALSE)
{
    hub <- ExperimentHub::ExperimentHub()
    gse <- AnnotationHub::query(hub, "GSE62944")
    tum <- gse[["EH1043"]]
    
    ind <- grep("CancerType", colnames(colData(tum)))
    colnames(colData(tum))[ind] <- "type"

    # should clinical variables be dropped?
    if(!with.clin.vars)
    {
        colData(tum) <- colData(tum)[,c(1,ind)]
        colData(tum)[,1] <- colnames(tum)
        colnames(colData(tum))[1] <- "sample"
    }

    norm <- gse[["EH1044"]]
    names(assays(norm)) <- names(assays(tum)) <- "exprs" 
    
    res <- list(tum=tum, norm=norm)
    return(res)
}

#
# alternative way of getting GSE62944 in case ExperimentHub fails
#
.downloadGSE62944 <- function(data.dir, with.clin.vars=FALSE)
{
    message("Downloading GSE62944 from GEO ...")
    GEOquery::getGEOSuppFiles("GSE62944", baseDir=data.dir)
    data.dir <- file.path(data.dir, "GSE62944")
    message(paste("Data files are stored under:", data.dir))

    # feature counts
    tum.file <- "GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FeatureCounts.txt.gz"
    norm.file <- "GSM1697009_06_01_15_TCGA_24.normal_Rsubread_FeatureCounts.txt.gz"
   
    tar.file <- file.path(data.dir, "GSE62944_RAW.tar") 
    untar(tar.file, c(norm.file,  tum.file), exdir=data.dir)
    
    message("Reading feature counts ...")
    norm.file <- file.path(data.dir, norm.file)
    norm.cont <- read.delim(norm.file, row.names=1L, as.is=TRUE)
    tum.file <- file.path(data.dir, tum.file)
    tum.cont <- read.delim(tum.file, row.names=1L)
   
    colnames(norm.cont) <- gsub("\\.", "-", colnames(norm.cont))
    colnames(tum.cont) <- gsub("\\.", "-", colnames(tum.cont))

    # cancer types
    message("Reading sample annotation (cancer/normal) ...")
    tum.cl.file <- "GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt.gz"
    tum.cl.file <- file.path(data.dir, tum.cl.file)
    norm.cl.file <- "GSE62944_06_01_15_TCGA_24_Normal_CancerType_Samples.txt.gz"
    norm.cl.file <- file.path(data.dir, norm.cl.file)
    tum.cl <- read.delim(tum.cl.file, header=FALSE, as.is=TRUE)
    norm.cl <- read.delim(norm.cl.file, header=FALSE, as.is=TRUE)
    colnames(tum.cl) <- colnames(norm.cl) <- c("sample", "type")

    # clinical variables
    if(with.clin.vars)
    {
        clin.var.file <- "GSE62944_06_01_15_TCGA_24_548_Clinical_Variables_9264_Samples.txt.gz"
        clin.var.file <- file.path(data.dir, clin.var.file)
        clin.var <- read.delim(clin.var.file, as.is=TRUE) 
        n <- clin.var[,1] 
        clin.var <- clin.var[,-c(1,2,3)]
        colnames(clin.var) <- gsub("\\.", "-", colnames(clin.var))
        clin.var <- t(clin.var)
        colnames(clin.var) <- n
        tum.cl <- cbind(clin.var, tum.cl[,"type"])
        colnames(tum.cl)[ncol(tum.cl)] <- "type" 
    }

    # clean up
    message("Cleaning up ...")
    rel.files <- c(tum.file, tum.cl.file, norm.file, norm.cl.file)
    if(with.clin.vars) rel.files <- c(rel.files, clin.var.file)
    rm.files <- list.files(data.dir, full.names=TRUE) 
    rm.files <- rm.files[!(rm.files %in% rel.files)]
    file.remove(rm.files)

    # create SummarizedExperiment for tumor and normal 
    message("Creating a SummarizedExperiment for tumor and normal samples ...")
    tum <- SummarizedExperiment(assays=list(exprs=as.matrix(tum.cont)))
    colData(tum) <- DataFrame(tum.cl)    
 
    norm <- SummarizedExperiment(assays=list(exprs=as.matrix(norm.cont)))
    colData(norm) <- DataFrame(norm.cl)    

    res <- list(tum=tum, norm=norm)
    return(res)
}

.splitByCancerType <- function(tum, norm, 
    nr.datasets, data.dir, min.ctrls=9, min.cpm=2)
{
    if(!file.exists(data.dir)) dir.create(data.dir, recursive=TRUE)
    # map IDs hgnc -> entrez
    tum <- EnrichmentBrowser::map.ids(tum, org="hsa", from="SYMBOL", to="ENTREZID")
    norm <- EnrichmentBrowser::map.ids(norm, org="hsa", from="SYMBOL", to="ENTREZID")

    # creating SEs
    tum.cts <- table(colData(tum)[,"type"])
    message("Cancer types with tumor samples:")
    message(paste(names(tum.cts), collapse=", "))
    norm.cts <- table(colData(norm)[,"type"])
    message("Cancer types with adj. normal samples:")
    message(paste(names(norm.cts), collapse=", "))

    norm.cts.min.ctrls <- norm.cts[norm.cts >= min.ctrls]
    rel.cts <- names(norm.cts.min.ctrls)[names(norm.cts.min.ctrls) %in% names(tum.cts)] 

    if(!is.null(nr.datasets)) rel.cts <- rel.cts[seq_len(nr.datasets)] 
    message("Cancer types with sufficient tumor and adj. normal samples:")
    message(paste(rel.cts, collapse=", "))
    message("Creating a SummarizedExperiment for each of them ...")

    exp.list <- lapply(rel.cts,
        function(ct) 
        {
            message(paste(ct, "tumor:", tum.cts[ct], "adj.normal:", norm.cts[ct]))

            ctrls <- norm[, norm$type == ct]
            cases <- tum[, tum$type == ct]

            se <- cbind(cases, ctrls) 
            colData(se)$GROUP <- c( rep(1, ncol(cases)), rep(0, ncol(ctrls)) )
            metadata(se)$dataId <- ct 
            metadata(se)$annotation <- "hsa"
            metadata(se)$dataType <- "rseq"

            # preprocessing
            # rm genes with all 0 and less than min.reads
            rs <- rowSums(edgeR::cpm(assay(se)) > min.cpm)
            keep <-  rs >= ncol(se) / 2
            se <- se[keep,]
   
            out.file <- paste(ct, "rds", sep=".")
            out.file <- file.path(data.dir, out.file) 
            saveRDS(se, file=out.file)
            return(se)
        })
    names(exp.list) <- rel.cts
    return(exp.list)
}
