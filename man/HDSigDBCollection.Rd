\name{HDSigDBCollection}
\alias{HDSigDBCollection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Returns a collection of expression pattern gene sets curated from public data related to Huntington's
Disease
}
\description{
  Returns a collection of expression pattern gene sets curated from public data related to Huntington's
Disease.
}
\usage{
HDSigDBCollection(
        organism = "human", 
        useHomology = TRUE,
        addOldOrganismToSetNames = FALSE,
        namePattern = ".convertedFrom.\%o",
        addOldOrganismToSetDescriptions = FALSE,
        descriptionPattern = " (Converted from \%o.)",
        ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{organism}{
   Character string specifying the organism. Can be common name (such as "human"), scientific name ("Homo
sapiens"), or scientific shorthand ("Hs"). To see all recognized organisms, use function
\code{\link{organismLabels}}.  If \code{NULL} or \code{NA}, gene sets in the collection will retain their
original content and organism they correspond to.
}
  \item{useHomology}{logical: should genes in organisms other than human and mice by matched by homology
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

  \item{...}{Unused.}
}
\details{
  Gene sets in this collection are based on expression patterns (differential expression, region or cell type
specific expression etc) in public data. Most of the studies compiled in this collection 
relate in some shape or form to Huntington's Disease. 

The genes were originally compiled for the human and
mouse genomes. If another
organism is specified, the gene sets and gene properties will be converted to the new organism (see
\code{convertCollectionToOrganism}). Note that each such conversion may introduce some differences (for
example, drop genes that cannot be mapped), and the gene set interpretation in the new organism may not be
as accurate as in the original organism. 

For detailed description and references for each individual gene set or gene property in the internal
collection, the user may consider using the function \code{\link{collection2dataFrames}} to create a data
frame that contains meta-information about gene sets (gene set names, descriptions, sources (references)
etc.) 
}
\value{
  Collection of gene lists and properties.
}
\author{
  Data were analyzed and gene sets compiled by Rancho Biosciences under contract from CHDI, Inc. 
  Code by Peter Langfelder.
}
\examples{
  \dontrun{
  collection = HDSigDBCollection(organism = "human");
  knownGroups(collection)
  }

}
\keyword{misc}% __ONLY ONE__ keyword per line
