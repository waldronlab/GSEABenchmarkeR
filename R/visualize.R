############################################################
# 
# author: Ludwig Geistlinger
# date: 2016-06-29 12:07:36
# 
# descr: plotting of eval results
# 
############################################################

# general boxplot functionality for enrichment method data
bpPlot <- function(data, what=c("runtime", "sig.sets", "rel.sets")) 
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
                    what)

    par(las=2)
    boxplot(data, col=rainbow(length(data)), ylab=ylab)
}

# plots aggregated relscoresum distribution of enrichment methods over all datasets
plotOverallRankDistrib <- function(k, best="max", ylab="%opt")
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
barplotRelScores <- function(x, file)
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

compRankDistrib <- function(x)
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

