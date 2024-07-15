\name{convertCollectionToOrganism}
\alias{convertCollectionToOrganism}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert a collection to another organism.
}
\description{
Converts Entrez IDs in gene sets and gene properties between organisms.
}
\usage{
convertCollectionToOrganism(
   collection, 
   organism, 
   useHomology = TRUE,
   addOldOrganismToSetNames = FALSE,
   namePattern = ".convertedFrom.\%o",
   addOldOrganismToSetDescriptions = FALSE,
   descriptionPattern = " (Converted from \%o.)",
   dropEmptySets = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{collection}{
  Collection to be converted.
}
  \item{organism}{
  The target organism to which to convert the collection. Organism can be specified
by a common name, scientific name, or shorthand. See \code{\link{organismLabels}} for known organisms.
}
  \item{useHomology}{
  Logical: should homology information be used? See details.
}
  \item{addOldOrganismToSetNames}{Logical: should an indication of the old organism be appended to the set
name? This may be useful to indicate that the gene sets was originally constructed for a different organism
and its content may not accurately reflect  the new organism.}

  \item{namePattern}{Character string giving the pattern that is to be added to
the set name. The characters "\%o" will be substituted by the actual organism name.}

  \item{addOldOrganismToSetDescriptions}{Logical: should an  indication of the old organism be appended to the
set description?}

  \item{descriptionPattern}{Character string giving the pattern that is to be added to
the set description. The characters "\%o" will be substituted by the actual organism name.}

  \item{dropEmptySets}{Logical: should empty sets be dropped from the returned collection?}
}
\details{
Entrez IDs can be converted using homology information, if such information is available. At present
homology information from Mouse Genome Informatics (http://informatics.jax.org/) is used for homology
mapping. If homology mapping is not available or not selected, genes are matched by first converting Entrez
IDs to gene symols, then looking up Entrez IDs for the same gene symbols (disregarding upper/lowecase) for
the target organism. It is recommended that homology be used whenever available.

If a gene cannot be mapped to the new organism, it is removed from the resulting collection; thus, the
number of genes in each gene set or in the gene properties may be lower.
}
\value{
A collection containing converted gene Entrez IDs.
}
\references{
 MGI homology mapping data are available at ftp://ftp.informatics.jax.org/pub/reports/index.html#homology 
}
\author{
Peter Langfelder
}
\seealso{
 \code{\link{convertGeneSetToOrganism}} for converting individual gene sets to another organism;

 \code{\link{mapEntrez}} for the workhorse Entrez ID mapping function;

 \code{\link{organismLabels}} for a list of known organisms and the various ways of identifying them.
}
\keyword{misc}% __ONLY ONE__ keyword per line
