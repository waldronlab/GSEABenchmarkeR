\name{deApply}
\alias{deApply}
\alias{writeDE}
\title{ 
    Differential expression analysis for datasets of a compendium
}
\description{
    This function applies selected methods for differential expression (DE)
    analysis to selected datasets of an expression data compendium.
}
\usage{
    deApply( data.ids, 
                de.method=c("limma", "edgeR", "DESeq"), 
                padj.method="flexible", parallel=NULL, ... )

    writeDE( data.ids, out.dir=NULL )
}
\arguments{
    \item{data.ids}{
        A character vector of dataset IDs.
    }
    \item{de.method}{
        Differential expression method.
        See documentation of \code{EnrichmentBrowser::de.ana}.
    }
    \item{padj.method}{
        Method for adjusting p-values to multiple testing. 
        For available methods see the man page of the stats
        function \code{p.adjust}.  
        Defaults to 'flexible', which applies a dataset-specific
        correction strategy. See details.   
    }
    \item{parallel}{
        Parallel computation mode.
        An instance of class \code{\linkS4class{BiocParallelParam}}.
        See the vignette of the \code{BiocParallel} package for switching 
        between serial, multi-core, and grid execution.
        Defaults to \code{NULL}, which then uses the first element of 
        \code{BiocParallel::registered()} for execution.
        If not changed by the user, this accordingly defaults to multi-core
        execution on the local host. 
    }
    \item{out.dir}{
        Character.
        Determines the output directory where DE results for each
        dataset are written to. 
        Defaults to \code{NULL}, which then writes to a subdir named
        'de' in \code{rappdirs::user_data_dir("GSEABenchmarkeR")}.
    }
    \item{...}{
        Additional arguments passed to \code{EnrichmentBrowser::de.ana}.
    }
}
\details{
    Large-scale DE studies typically report a gene as differentially expressed
    if the corresponding DE p-value, corrected for multiple testing, satisfies 
    the chosen significance level.
    In this case, the method used for multiple testing correction can substantially 
    influence the performance of enrichment methods that work directly on the list of
    DE genes.

    An example is the frequently used over-representation analysis (ORA), which 
    assesses the overlap between the DE genes and a gene set under study based
    on the hypergeometric distribution (see Appendix A of the 
    \code{EnrichmentBrowser} vignette for an introduction).

    In the extreme case, where there are no, or only a few, genes satisfying the 
    significance threshold, this can even render ORA inapplicable.
    On the other hand, application of ORA is useless if almost all genes turn out to  
    be DE. 

    Using \code{padj.method="flexible"} accounts for these cases by applying 
    multiple testing correction in dependence on the degree of differential 
    expression: 

    \itemize{
        \item the correction method from Benjamini and Hochberg (BH) is applied
                if it renders >= 1\% and <= 25\% of all measured genes as DE, 
        \item the p-values are left unadjusted, if the BH correction results 
                in < 1\% DE genes, and 
        \item the more stringent Bonferroni correction is applied, if the BH 
                correction results in > 25\% DE genes.
    } 

    Note that resulting p-values should not be used for assessing the 
    statistical significance of DE genes within or between datasets. 
    The intended use is to determine which genes to include in the analysis
    with ORA - and thus to ensure that the fraction of included genes is 
    roughly in the same order of magnitude across datasets. 
    
    Alternative stratgies could also be applied - such as taking a constant 
    number of genes for each dataset or excluding ORA methods in general from 
    the assessment. 
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
\seealso{
    \code{loadEData} to load a specified expression data compendium.
}
\examples{
    # reading user-defined expression data from file
    data.dir <- system.file("extdata/myEData", package="GSEABenchmarkeR")
    data.ids <- loadEData(data.dir)

    # differential expression analysis
    deApply(data.ids)
    se <- getDataset(data.ids[1])

    # writing DE results to file
    out.dir <- tempdir()
    out.dir <- file.path(out.dir, "de")
    if(!file.exists(out.dir)) dir.create(out.dir)
 
    writeDE(data.ids, out.dir)    

}   

