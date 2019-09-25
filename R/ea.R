#' Differential expression analysis for datasets of a compendium
#' 
#' This function applies selected methods for differential expression (DE)
#' analysis to selected datasets of an expression data compendium.
#' 
#' DE studies typically report a gene as differentially expressed if the
#' corresponding DE p-value, corrected for multiple testing, satisfies the
#' chosen significance level.  Enrichment methods that work directly on the
#' list of DE genes are then substantially influenced by the multiple testing
#' correction.
#' 
#' An example is the frequently used over-representation analysis (ORA), which
#' assesses the overlap between the DE genes and a gene set under study based
#' on the hypergeometric distribution (see Appendix A of the
#' \code{EnrichmentBrowser} vignette for an introduction).
#' 
#' ORA is inapplicable if there are few genes satisfying the significance
#' threshold, or if almost all genes are DE.
#' 
#' Using \code{padj.method="flexible"} accounts for these cases by applying
#' multiple testing correction in dependence on the degree of differential
#' expression:
#' 
#' \itemize{ \item the correction method from Benjamini and Hochberg (BH) is
#' applied if it renders >= 1\% and <= 25\% of all measured genes as DE, \item
#' the p-values are left unadjusted, if the BH correction results in < 1\% DE
#' genes, and \item the more stringent Bonferroni correction is applied, if the
#' BH correction results in > 25\% DE genes.  }
#' 
#' Note that resulting p-values should not be used for assessing the
#' statistical significance of DE genes within or between datasets.  They are
#' solely used to determine which genes are included in the analysis with ORA -
#' where the flexible correction ensures that the fraction of included genes is
#' roughly in the same order of magnitude across datasets.
#' 
#' Alternative stratgies could also be applied - such as taking a constant
#' number of genes for each dataset or excluding ORA methods in general from
#' the assessment.
#' 
#' @param exp.list Experiment list.  A \code{list} of datasets, each being of
#' class \code{\linkS4class{SummarizedExperiment}}.
#' @param de.method Differential expression method.  See documentation of
#' \code{\link{deAna}}.
#' @param padj.method Method for adjusting p-values to multiple testing.  For
#' available methods see the man page of the stats function \code{p.adjust}.
#' Defaults to 'flexible', which applies a dataset-specific correction
#' strategy. See details.
#' @param parallel Parallel computation mode.  An instance of class
#' \code{\linkS4class{BiocParallelParam}}.  See the vignette of the
#' \code{BiocParallel} package for switching between serial, multi-core, and
#' grid execution.  Defaults to \code{NULL}, which then uses the first element
#' of \code{BiocParallel::registered()} for execution.  If not changed by the
#' user, this accordingly defaults to multi-core execution on the local host.
#' @param alpha Statistical significance level. Defaults to 0.05.
#' @param beta Absolute log2 fold change cut-off. Defaults to 1 (2-fold).
#' @param out.dir Character.  Determines the output directory where DE results
#' for each dataset are written to.  Defaults to \code{NULL}, which then writes
#' to a subdir named 'de' in \code{rappdirs::user_data_dir("GSEABenchmarkeR")}.
#' @param max.na Integer. Determines for which genes a meta fold change is 
#' computed. Per default, excludes genes for which the fold change is not 
#' annotated in >= 1/3 of the datasets in \code{exp.list}.
#' @param ...  Additional arguments passed to \code{EnrichmentBrowser::deAna}.
#' @return \code{runDE} returns \code{exp.list} with DE measures annotated to
#' the \code{\link{rowData}} slot of each dataset, \code{writeDE} writes to file,
#' and \code{plotDEDistribution} plots to a graphics device.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{loadEData} to load a specified expression data compendium.
#' @examples
#' 
#'     # reading user-defined expression data from file
#'     data.dir <- system.file("extdata/myEData", package="GSEABenchmarkeR")
#'     edat <- loadEData(data.dir)
#' 
#'     # differential expression analysis
#'     edat <- runDE(edat)
#' 
#'     # visualization of per-dataset DE distribution
#'     plotDEDistribution(edat)
#' 
#'     # calculating meta fold changes across datasets 
#'     mfcs <- metaFC(edat, max.na=0) 
#' 
#'     # writing DE results to file
#'     out.dir <- tempdir()
#'     out.dir <- file.path(out.dir, "de")
#'     if(!file.exists(out.dir)) dir.create(out.dir)
#'  
#'     writeDE(edat, out.dir)    
#' 
#' @export runDE
runDE <- function(exp.list, 
    de.method=c("limma", "edgeR", "DESeq2"), 
    padj.method="flexible", parallel=NULL, ...)
{
    flex <- padj.method == "flexible"
    if(flex) padj.method <- "none"
    
    PCOL <- EnrichmentBrowser::configEBrowser("PVAL.COL")
    ADJP.COL <- EnrichmentBrowser::configEBrowser("ADJP.COL")
    .de <- function(i, ...)
    { 
        se <- EnrichmentBrowser::deAna(i, de.method=de.method, 
                                            padj.method=padj.method, ...)

        if(flex)
        {
            padj <- p.adjust(rowData(se)[,PCOL], method="BH")
            fracSigP <- mean(padj < 0.05, na.rm=TRUE)
            if(fracSigP > 0.25) 
                padj <- p.adjust(rowData(se)[,PCOL], method="bonf") 
            if(fracSigP > 0.01) rowData(se)[,ADJP.COL] <- padj
            else rowData(se)[,ADJP.COL] <- rowData(se)[,PCOL]
        }
        return(se)
    }
   
    exp.list <- .iter(exp.list, .de, ..., parallel=parallel)
    return(exp.list)
}

#' @rdname runDE
#' @export
metaFC <- function(exp.list, max.na=round(length(exp.list) / 3))
{
    all.names <- unique(unlist(lapply(exp.list, names)))
    glen <- length(all.names)

    # fold changes
    FC.COL <- EnrichmentBrowser::configEBrowser("FC.COL")
    fcs <- vapply(exp.list, function(d) 
                    rowData(d,use.names=TRUE)[all.names,FC.COL], 
                    numeric(glen))
    rownames(fcs) <- all.names

    # filter out genes not measured in >= 1/3 of datasets
    few.nas <- rowSums(t(apply(fcs, 1, is.na))) <= max.na
    fcs <- fcs[few.nas,]

    # weights ~ sample size
    weights <- vapply(exp.list, ncol, integer(1))
    weights <- sqrt(weights) 

    # weighted mean
    wfcs <- apply(fcs, 1, weighted.mean, w=weights, na.rm=TRUE)  
    wfcs <- wfcs[order(abs(wfcs), decreasing=TRUE)]
    return(wfcs)
}

#
# fraction of DE genes according to p-value and fold change threshold
#
.fractDE <- function(se, alpha=0.05, beta=1, freq=c("rel", "abs"))
{
    freq <- match.arg(freq)
    FC.COL <- EnrichmentBrowser::configEBrowser("FC.COL")
    ADJP.COL <- EnrichmentBrowser::configEBrowser("ADJP.COL")

    fract.p <- sum(rowData(se)[,ADJP.COL] < alpha)
    fract.fc <- sum(abs(rowData(se)[,FC.COL]) > beta)
    
    if(freq == "rel")
    { 
        fract.p <- fract.p / nrow(se) * 100
        fract.fc <- fract.fc / nrow(se) * 100
    }
    res <- c(fract.p, fract.fc)
    names(res) <- c("p", "fc")
    return(res)
}


#' Application of enrichment methods to multiple datasets
#' 
#' This function applies selected methods for enrichment analysis to selected
#' datasets of a compendium.
#' 
#' 
#' @param exp.list Experiment list.  A \code{list} of datasets, each being of
#' class \code{\linkS4class{SummarizedExperiment}}. In case of just one dataset
#' a single \code{\linkS4class{SummarizedExperiment}} is also allowed. 
#' See the documentation of \code{\link{sbea}} for required minimal annotations.
#' @param methods Methods for enrichment analysis.  A character vector with
#' method names chosen from \code{\link{sbeaMethods}} and
#' \code{\link{nbeaMethods}}, or user-defined functions
#' implementing methods for enrichment analysis.
#' @param gs Gene sets, i.e. a list of character vectors of gene IDs.
#' @param perm Number of permutations of the sample group assignments. 
#' Defaults to 1000. Can also be an integer vector matching
#' the length of \code{methods} to assign different numbers of permutations for
#' different methods.
#' @param parallel Parallel computation mode.  An instance of class
#' \code{\linkS4class{BiocParallelParam}}.  See the vignette of the
#' \code{BiocParallel} package for switching between serial, multi-core, and
#' grid execution.  Defaults to \code{NULL}, which then uses the first element
#' of \code{BiocParallel::registered()} for execution.  If not changed by the
#' user, this accordingly defaults to multi-core execution on the local host.
#' @param save2file Logical. Should results be saved to file for subsequent
#' benchmarking?  Defaults to \code{FALSE}.
#' @param out.dir Character.  Determines the output directory where results are
#' saved to.  Defaults to \code{NULL}, which then writes to
#' \code{rappdirs::user_data_dir("GSEABenchmarkeR")} in case \code{save2file}
#' is set to \code{TRUE}.
#' @param ...  Additional arguments passed to the selected enrichment methods.
#' @return A list with an entry for each method applied.  Each method entry is
#' a list with an entry for each dataset analyzed.  Each dataset entry is a
#' list of length 2, with the first element being the runtime and the second
#' element being the gene set ranking, as obtained from applying the respective
#' method to the respective dataset.
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
#'     gs.file <- system.file("extdata/hsa_kegg_gs.gmt", package="EnrichmentBrowser")
#'     kegg.gs <- EnrichmentBrowser::getGenesets(gs.file)
#' 
#'     # applying two methods to two datasets 
#'     res <- runEA(geo2kegg, methods=c("ora", "camera"), gs=kegg.gs, perm=0)
#'     
#' 
#' @export runEA
runEA <- function(exp.list, methods, gs, perm=1000,
    parallel=NULL, save2file=FALSE, out.dir=NULL, ...)
{
    # singleton call?
    if(is(exp.list, "SummarizedExperiment"))
    {
        if(length(methods) == 1)
        {
            res <- .ea(exp.list, methods, gs, perm[1], 
                        save2file=save2file, out.dir=out.dir, ...)
            return(res)
        }
        exp.list <- list(exp.list)
    }

    # setup
    if(is.function(methods)) methods <- list(method=methods) 
    else
    {   
        .eaPkgs(methods)
        names(methods) <- methods
    }
    nr.meth <- length(methods)

    show.progress <- interactive() && nr.meth > 2
    if(show.progress) pb <- txtProgressBar(0, nr.meth, style=3)

    if(length(perm) != nr.meth) perm <- rep(perm[1], nr.meth)
    names(perm) <- names(methods)
  
    for(i in seq_along(exp.list)) 
        metadata(exp.list[[i]])$dataId <- names(exp.list)[i]

    res <- lapply(names(methods), 
        function(m, ...)
        {
            if(show.progress) setTxtProgressBar(pb, match(m, methods))
            r <- .iter(exp.list, .ea, 
                        method=methods[[m]], perm=perm[m], ..., parallel=parallel)
            return(r)
        },
        gs=gs, save2file=save2file, out.dir=out.dir, ...
    )
   
    if(show.progress) close(pb) 
    names(res) <- names(methods)
    return(res)
}

# check on availability of ea packages
.eaPkgs <- function(ea.methods)
{
    if(is.list(ea.methods)) ea.methods <- names(ea.methods)
    sbea.pkgs <- EnrichmentBrowser::configEBrowser("SBEA.PKGS")
    nbea.pkgs <- EnrichmentBrowser::configEBrowser("NBEA.PKGS")
    ea.pkgs <- c(sbea.pkgs, nbea.pkgs)
    for(m in ea.methods)
        if(m %in% names(ea.pkgs)) 
            EnrichmentBrowser::isAvailable(ea.pkgs[m], type="software")
}

# execute indiviual methods
.ea <- function(se, method, gs, perm=1000, save2file=FALSE, out.dir=NULL, ...)
{
    id <- metadata(se)$dataId
    res <- .execEA(method, se, gs, perm, ...)
    ti <- res$ti
    res <- res$res

    if(is.function(method)) method <- "method"
    if(is(res, "try-error"))
    {
        message(paste(method, "could not be evaluated on", id))
        message(res)
        message("Returning NULL")
        return(NULL)
    }

    # output
    res <- EnrichmentBrowser::gsRanking(res, signif.only=FALSE)

    if(save2file) .save2file(res, out.dir, method, id, ti[3])
    return(list(runtime=ti[3], ranking=res))
}

.execEA <- function(method, se, gs, perm, ...)
{
    f <- file()
    sink(file=f) ## silence

    # execute
    suppressMessages(
    if(is.function(method) || method %in% EnrichmentBrowser::sbeaMethods())
    {
        ti <- system.time(
            res <- try(
                EnrichmentBrowser::sbea(method, se, gs, perm=perm, ...), 
            silent=TRUE)
        )
    }
    else
    { 
        ti <- system.time(
            res <- try(
                EnrichmentBrowser::nbea(method, se, gs, perm=perm, ...),
            silent=TRUE)
        )
    })
    sink() ## undo silencing
    close(f)
    return(list(res=res,ti=ti))
}

.save2file <- function(res, out.dir, method, id, ti=NULL)
{
    if(is.null(out.dir)) 
            out.dir <- rappdirs::user_data_dir("GSEABenchmarkeR")
    out.dir <- file.path(out.dir, method)
    if(!file.exists(out.dir)) dir.create(out.dir, recursive=TRUE)

    out.file <- file.path(out.dir, paste0(id, ".rds"))
    saveRDS(res, file=out.file)
    if(!is.null(ti)) cat(ti, file=sub("rds$", "txt", out.file))
}


