\name{HuntingtonsDiseaseWGCNACollection}
\alias{HuntingtonsDiseaseWGCNACollection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Returns the collection of WGCNA modules from Huntington's Disease (HD)-related data.
}
\description{
  Returns the collection of WGCNA modules from Huntington's Disease (HD)-related data.
}
\usage{
HuntingtonsDiseaseWGCNACollection(
        organism = "human", 
        useHomology = TRUE,
        addOldOrganismToSetNames = FALSE,
        namePattern = ".convertedFrom.\%o",
        addOldOrganismToSetDescriptions = FALSE,
        descriptionPattern = " (Converted from \%o.)")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{organism}{
   Character string specifying the organism. Can be common name (such as "human"), scientific name ("Homo
sapiens"), or scientific shorthand ("Hs"). To see all recognized organisms, use function
\code{\link{organismLabels}}. If \code{NULL} or \code{NA}, gene sets in the collection will retain their
original content and organism they correspond to.
}
  \item{useHomology}{logical: should genes in non-human organisms by matched by homology
information (where available)? If this is \code{FALSE} or the homology information is not available, genes
will be matched by name.}

  \item{addOldOrganismToSetNames}{Logical: should an indication of the old organism be appended to the set
name? This may be useful to indicate that the gene sets was originally constructed for a different organism
and its content may not accurately reflect  the new organism.}

  \item{namePattern}{Character string giving the pattern that is to be added to
the set name. The characters "\%o" will be substituted by the actual organism name.}

  \item{addOldOrganismToSetDescriptions}{Logical: should an  indication of the old organism be appended to the
set description?}

  \item{descriptionPattern}{Character string giving the pattern that is to be added to
the set description. The characters "\%o" will be substituted by the actual organism name.}
}
\details{
 The gene sets in the WGCNA HD collection consist of modules determined by WGCNA analyses of various
HD-related data, some from human and some from mouse tissues. This function converts all gene sets to the
specified \code{organism} using \code{convertCollectionToOrganism}.
Note that each such conversion may introduce some differences (for
example, drop genes that cannot be mapped), and the gene set interpretation in the new organism may not be
as accurate as in the original organism (human). 

For detailed description and references for each individual gene set or gene property in the WGCNA.HD.
collection, the user may consider using the function \code{\link{collection2dataFrames}} to create a data
frame that contains meta-information about gene sets (gene set names, descriptions, sources (references)
etc.) 
}
\value{
  Collection of gene lists and properties.
}
\author{
  WGCNA analyses were performed by Peter Langfelder based on publicly available data and data proprietary
to CHDI, Inc.
}
\examples{
  \dontrun{
  collection = HuntingtonsDiseaseWGCNACollection(organism = "human");
  knownGroups(collection)
  }

}
\keyword{misc}% __ONLY ONE__ keyword per line
