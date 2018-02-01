\name{loadEData}
\alias{loadEData}
\title{ 
    Loading pre-defined and user-defined expression data
}
\description{
    This function implements a general interface for loading the pre-defined 
    GEO2KEGG microarray compendium and the TCGA RNA-seq compendium.
    It also allows loading of user-defined data from file. 
}
\usage{

    loadEData( edata, reload=FALSE, nr.datasets=NULL, ... )
}
\arguments{
    \item{edata}{
        Expression data compendium.
        A character vector of length 1 that must be either
        \itemize{
            \item 'geo2kegg': to load the GEO2KEGG microarray compendium,
            \item 'tcga': to load the TCGA RNA-seq compendium, or
            \item an absolute file path pointing to a directory, in which a 
            user-defined compendium has been saved in RDS files. 
        }   
        See details.
    }
    \item{reload}{
        Logical.
        If the chosen compendium has already been loaded into the workspace, 
        should it be reloaded, thereby resetting all changes performed to the 
        compendium during the current R session?
        Defaults to FALSE.
    }
    \item{nr.datasets}{
        Integer.
        Number of datasets that should be loaded from the compendium.
        This is mainly for demonstration purposes.
    }
    \item{...}{
        Additional arguments passed to the internal loading routines of the 
        GEO2KEGG and TCGA compendia.
        This currently includes for loading of the GEO2KEGG compendium
        \itemize{
            \item preproc: logical. Should probe level data automatically be 
                summarized to gene level data? Defaults to FALSE.
            \item de.only: logical. Include only datasets in which 
                differentially expressed genes have been found? Defaults to FALSE.
            \item excl.metac: logical. Exclude datasets for which MetaCore rather
                than KEGG pathways have been assigned as target pathways?
                Defaults to FALSE.
        }
        And for loading of the TCGA compendium
        \itemize{
            \item mode: character, determines how GSE62944 is obtained. Either 
                'ehub' (default, via ExperimentHub) or 'geo' (direct download 
                from GEO, slow).
            \item data.dir: character. Absolute file path indicating where 
                processed RDS files for each dataset are written to. 
                Defaults to NULL, which will then write to 
                \code{rappdirs::user_data_dir("GSEABenchmarkeR")}.
            \item min.ctrls: integer. Minimum number of controls, i.e. adjacent 
                normal samples, for a cancer type to be included. Defaults to 9.
            \item min.cpm: integer. Minimum counts-per-million reads mapped. 
                    See the edgeR vignette for details. The default filter is to
                    exclude genes with cpm < 2 in more than half of the samples.
            \item with.clin.vars: logical. Should clinical variables (>500) be 
                kept to allow for more advanced sample groupings in addition to
                the default binary grouping (tumor vs. normal)?
        }
    }
}
\value{
    A character vector of dataset IDs.
}
\details{
    The pre-defined GEO2KEGG microarray compendium consists of 42 datasets 
    investigating a total of 19 different human diseases as collected by  
    Tarca et al. (2012 and 2013).

    The pre-defined TCGA RNA-seq compendium consists of datasets from 
    The Cancer Genome Atlas (TCGA, 2013) investigating a total of 34 different 
    cancer types.
 
    User-defined data can also be loaded, given that datasets, preferably of 
    class \code{SummarizedExperiment}, have been saved as \code{RDS} files.
}
\references{
    Tarca et al. (2012)
    Down-weighting overlapping genes improves gene set analysis.
    BMC Bioinformatics, 13:136.    

    Tarca et al. (2013)
    A comparison of gene set analysis methods in terms of sensitivity,
    prioritization and specificity.
    PLoS One, 8(11):e79217.

    The Cancer Genome Atlas Research Network (2013)
    The Cancer Genome Atlas Pan-Cancer analysis project.
    Nat Genet, 45(10):1113-20.

    Rahman et al. (2015)
    Alternative preprocessing of RNA-Sequencing data in The Cancer Genome Atlas
    leads to improved analysis results.
    Bioinformatics, 31(22):3666-72.
}
\author{
    Ludwig Geistlinger <Ludwig.Geistlinger@sph.cuny.edu>
}
\seealso{
    \code{\linkS4class{SummarizedExperiment}},

    \code{\link{getDataset}} to obtain an expression dataset from a loaded 
    compendium.
}
\examples{
    # (1) Loading the GEO2KEGG microarray compendium
    ma.ids <- loadEData("geo2kegg", nr.datasets=2)

    # (2) Loading the TCGA RNA-seq compendium
    rseq.ids <- loadEData("tcga", nr.datasets=2)

    # (3) reading user-defined expression data from file
    data.dir <- system.file("extdata/myEData", package="GSEABenchmarkeR")
    data.ids <- loadEData(data.dir)
}   

