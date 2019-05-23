############################################################
# 
# author: Ludwig Geistlinger
# date: 2016-06-29 12:07:36
# 
# descr: plotting of eval results
# 
############################################################

#' Customized boxplot visualization of benchmark results
#' 
#' This is a convenience function to create customized boxplots for specific
#' benchmark criteria such as runtime, statistical significance and phenotype
#' relevance.
#' 
#' 
#' @param data Numeric matrix or list of numeric vectors.  In case of a matrix,
#' column names are assumed to be method names and rownames are assumed to be
#' dataset IDs.  In case of a list, names are assumed to be methods names and
#' each element corresponds to a numeric vector with names assumed to be
#' dataset IDs.
#' @param what Character.  Determines how the plot is customized.  One of
#' \itemize{ \item runtime: displays runtime of methods across datasets, \item
#' sig.sets: displays percentage of significant gene sets, \item rel.sets:
#' displays phenotype relevance scores.  }
#' @return None. Plots to a graphics device.
#' @author Ludwig Geistlinger <Ludwig.Geistlinger@@sph.cuny.edu>
#' @seealso \code{\link{evalNrSigSets}} to evaluate fractions of significant 
#' gene sets; \code{\link{evalRelevance}} to evaluate phenotype relevance of
#' gene set rankings.
#' @examples
#' 
#'     # simulated setup:
#'     # 3 methods & 5 datasets
#'     methods <- paste0("m", 1:3)
#'     data.ids <- paste0("d", 1:5)
#' 
#'     # runtime data
#'     rt <- vapply(1:3, function(m) runif(5, min=m, max=m+1), numeric(5))
#'     rownames(rt) <- data.ids 
#'     colnames(rt) <- methods
#' 
#'     # plot
#'     bpPlot(rt, what="runtime")
#' 
#' @export bpPlot
bpPlot <- function(data, what=c("runtime", "sig.sets", "rel.sets", "typeI")) 
{
    if(is.matrix(data)) 
        data <- sapply(colnames(data), function(i) data[,i], simplify=FALSE) 
    names(data) <- substring(names(data), 1, 7)
    data <- data[order(vapply(data, median, numeric(1), na.rm=TRUE))]

    what <- what[1]
    if(what == "runtime") data <- lapply(data, log, base=10)
    
    ylab <- switch(what,
                    runtime = "log10 runtime [sec]",
                    sig.sets = "% significant sets",
                    rel.sets = "%opt",
                    typeI = "type I error rate",
                    what)

    par(las=2)
    boxplot(data, col=rainbow(length(data)), ylab=ylab)
    if(what == "typeI") abline(h=0.05, col="red", lty=2)
}

plotTypeIError <- function(data)
{
    names(data) <- substring(names(data), 1, 7)
    data <- data[,order(data["Max.",] - data["Mean",])]
    par(las=2)
    par(pch=20)
    boxplot(data, ylab="type I error rate")
    points(x=seq_len(ncol(data)), y=data["Mean",], col="darkviolet")
    abline(h=0.05, col="red", lty=2)
}

plotGSSizeRobustness <- function(mat, what=c("sig.sets", "typeI"))
{
    what <- what[1]
    ylab <- switch(what,
                    sig.sets = "% significant sets",
                    typeI = "type I error rate",
                    what)


    xmax <- nrow(mat)
    ymax <- max(mat, na.rm=TRUE)
    plot(NA, axes=FALSE, 
            xlim=c(1, xmax), ylim=c(0, ymax),
            xlab="GS size", ylab=ylab)

    grid <- sub("gs", "", rownames(mat))
    axis(1, at=seq_len(xmax), labels=grid)
    axis(2)
    box()
    matlines(mat)
}   


#' @rdname runDE
#' @export
plotDEDistribution <- function(exp.list, alpha=0.05, beta=1)
{
    x <- lapply(exp.list, function(se) .fractDE(se, alpha=alpha, beta=beta))
    x <- do.call(cbind, x)
    par(pch=20)
    plot(x=x["fc",],y=x["p",], 
        xlab="%[abs(log2FC) > 1]", ylab="%[adjp < 0.05]", col="red")
    text(x=x["fc",], y=x["p",], names(exp.list), cex=0.8, pos=4)
}

#' @rdname runDE
#' @export
plotNrSamples <- function(exp.list)
{
    GRP.COL <- EnrichmentBrowser::configEBrowser("GRP.COL")
    nr.ctrls <- vapply(exp.list, function(se) sum(se[[GRP.COL]] == 0), numeric(1))
    nr.cases <- vapply(exp.list, function(se) sum(se[[GRP.COL]] == 1), numeric(1))
    par(pch=20)
    plot(x=nr.ctrls, y=nr.cases, xlab="#controls", ylab="#cases", 
        col="red", xlim=c(0, max(nr.ctrls)), ylim=c(0, max(nr.cases)))
    text(x=nr.ctrls, y=nr.cases, names(exp.list), cex=0.8, pos=4)
}

# plots aggregated relscoresum distribution of enrichment methods over all datasets
.plotOverallRankDistrib <- function(k, best="max", ylab="%opt")
{
    k <- k[,order(apply(k,2,median, na.rm=TRUE), decreasing=(best=="max"))]
    par(las=2)
    colnames(k) <- substring(colnames(k), 1, 7)
    boxplot(k, col=rainbow(11), ylab=ylab)
    par(las=1)
    f <- ifelse(best == "max", which.max, which.min)
    sel <- apply(k, 1, f)
    best <- colnames(k)[sel]
    axis(3, at=seq_len(ncol(k)), 
        labels=vapply(colnames(k), function(n) sum(best==n), integer(1)))
    mtext("#best", side=3, line=3, las=1, cex=1)
}

# barplot relscores of enrichment methods per dataset
.barplotRelScores <- function(x, file)
{
    opt <- x[,"opt"]
    x <- x[,-ncol(x)]
    pdf(file)
    par(las=2)
    for(n in rownames(x))
    {
        dat <- sort(x[n, ], decreasing=TRUE)
        barplot(dat, col="lightblue", ylab="weighted relScore sum", 
            main=paste(n, "opt:", round(opt[n], digits=1)))
    }
    dev.off()
}

.compRankDistrib <- function(x)
{
    res <- apply(x, 1, 
        function(l) 
        {
            o <- order(l, decreasing=TRUE)
            r <- match(colnames(x), colnames(x)[o])
            return(r)
        })
    rownames(res) <- colnames(x)
    res <- t(res)
    res <- res[, order(apply(res, 2, median))] 
    return(res)
}

