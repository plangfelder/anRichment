\name{PhenopediaCollection}
\alias{PhenopediaCollection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Returns the collection of disease-associated genes collected in the CDC Phenopedia database
}
\description{
  Returns the collection of disease-associated genes collected in the CDC Phenopedia database.
}
\usage{
PhenopediaCollection(
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

  \item{...}{Other arguments, currently ignored.}
}
\details{

 Phenopedia database is part of the Public Health Genomics and Precision Health Knowledge Base compiled and maintained by
Centers for Disease Control and Prevention. The database collects human gene-disease and -trait associations by curating
published literature. 

The database version 5.6 was downloaded May 21, 2019.

This function converts all gene sets to the
specified \code{organism} using \code{convertCollectionToOrganism}.
Note that each such conversion may introduce some differences (for
example, drop genes that cannot be mapped), and the gene set interpretation in the new organism may not be
as accurate as in the original organism (human). 

For detailed description and references for each individual gene set or gene property in the 
collection, the user may consider using the function \code{\link{collection2dataFrames}} to create a data
frame that contains meta-information about gene sets (gene set names, descriptions, sources (references)
etc.). 
}
\value{
  Collection of gene lists and properties.
}
\references{
  CDC Phenopedia web page, https://phgkb.cdc.gov/PHGKB/startPagePhenoPedia.action

   Wei Yu, Melinda Clyne, Muin J. Khoury and Marta Gwinn. Phenopedia and Genopedia: Disease-centered and Gene-centered
Views of the Evolving Knowledge of Human Genetic Associations. Bioinformatics Vol. 26 no. 1 2010, pages 145-146; doi:10.1093/bioinformatics/btp618. 
}
\author{
  Data from CDC, code by Peter Langfelder  
}
\examples{
  \dontrun{
  collection = PhenopediaCollection(organism = "human");
  as.character(dataSetNames(collection)[1:100]);
  knownGroups(collection)
  }

}
\keyword{misc}% __ONLY ONE__ keyword per line
