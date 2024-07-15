\name{convertGeneSetToOrganism}
\alias{convertGeneSetToOrganism}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert gene Entrez IDs of a gene set to a new organism. 
}
\description{
This function converts the gene Entrez IDs of a gene set to a new organism.
}
\usage{
convertGeneSetToOrganism(
        geneSet, 
        organism, 
        entrezMap = NULL,
        useHomology = TRUE,
        addOldOrganismToName = FALSE,
        namePattern = ".convertedFrom.\%o",
        addOldOrganismToDescription = FALSE,
        descriptionPattern = " (Converted from \%o.)")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneSet}{
 The gene set to be converted.
}
  \item{organism}{
 Character string giving the organism that the gene set is to be converted to. Organism can be specified by
a common name, scientific name, or shorthand. See \code{\link{organismLabels}} for known organisms.
}

  \item{entrezMap}{An optional 2-column matrix or data frame giving the mapping between identifiers for old
and new organism in columns 1 and 2, respectively. If not given, the map will be generated as described in
Details. }

  \item{useHomology}{
   Logical: should homology information be used? See details.
}
  \item{addOldOrganismToName}{Logical: should an indication of the old organism be appended to the set
name? This may be useful to indicate that the gene sets was originally constructed for a different organism
and its content may not accurately reflect  the new organism.}

  \item{namePattern}{Character string giving the pattern that is to be added to
the set name. The characters "\%o" will be substituted by the actual organism name.}

  \item{addOldOrganismToDescription}{Logical: should an  indication of the old organism be appended to the
set description?}

  \item{descriptionPattern}{Character string giving the pattern that is to be added to
the set description. The characters "\%o" will be substituted by the actual organism name.}
}
\details{
Entrez IDs can be converted using homology information, if such information is available. At present
homology information from Mouse Genome Informatics (http://informatics.jax.org/) is used for homology
mapping. If homology mapping is not available or not selected, genes are matched by first converting Entrez
IDs to gene symols, then looking up Entrez IDs for the same gene symbols (disregarding upper/lowecase) for
the target organism. It is recommended that homology be used whenever available.

If a gene cannot be mapped to the new organism, it is removed from the resulting gene set; thus, the
number of genes in each gene set or in the gene properties may be lower.
}
\value{
 A gene set with gene Entrez IDs converted to the new ogranism, and the \code{organism} component set
appropriately.
}
\references{
 MGI homology mapping data are available at ftp://ftp.informatics.jax.org/pub/reports/index.html#homology
}
\author{
Peter Langfelder
}
\seealso{
 \code{\link{convertCollectionToOrganism}} for converting entire collections to another organism;

 \code{\link{mapEntrez}} for the workhorse Entrez ID mapping function;

 \code{\link{organismLabels}} for a list of known organisms and the various ways of identifying them.
}
\keyword{ misc }% __ONLY ONE__ keyword per line
