\name{BrainDiseaseCollection}
\alias{BrainDiseaseCollection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Returns a collection containing gene sets implicated in neurological diseases.
}
\description{
  Returns a collection containing gene sets implicated in neurological diseases.
}
\usage{
BrainDiseaseCollection(
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
original content and organism  (in this case mouse).
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
SCSBrainCellType gene sets were downloaded on December 14, 2015 and correspond to version 5.0. 
}
\value{
  Collection of gene sets implicated in neurological diseases.
}
\author{
Code by Peter Langfelder, gene sets compiled by Jason Chen from on-line sources.
}
\examples{
  \dontrun{
  collection = BrainDiseaseCollection(organism = "mouse");
  knownGroups(collection)
  }

}
\keyword{misc}% __ONLY ONE__ keyword per line
