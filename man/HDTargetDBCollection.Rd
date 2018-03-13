\name{HDTargetDBCollection}
\alias{HDTargetDBCollection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Returns the collection of Huntington's Disease literature gene sets
}
\description{
  Returns the collection of Huntington's Disease literature gene sets.
}
\usage{
HDTargetDBCollection(
        organism = "human", 
        useHomology = TRUE,
        addOldOrganismToSetNames = FALSE,
        namePattern = ".convertedFrom.\%o",
        addOldOrganismToSetDescriptions = FALSE,
        descriptionPattern = " (Converted from \%o.)",
        traceable = TRUE,
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

  \item{traceable}{Logical: should only a subset of the gene sets that are traceable to a published
article or website be included?}

  \item{...}{Other arguments, currently ignored.}
}
\details{
 The gene sets in the HD literature collection are compiled from various literature sources. They represent
functional categories, sets of genes from perturbation studies and others. Although the genes may have been
studied in various model organisms, they are all stored as human orthologs.

This function converts all gene sets to the
specified \code{organism} using \code{convertCollectionToOrganism}.
Note that each such conversion may introduce some differences (for
example, drop genes that cannot be mapped), and the gene set interpretation in the new organism may not be
as accurate as in the original organism (human). 

For detailed description and references for each individual gene set or gene property in the 
collection, the user may consider using the function \code{\link{collection2dataFrames}} to create a data
frame that contains meta-information about gene sets (gene set names, descriptions, sources (references)
etc.) 
}
\value{
  Collection of gene lists and properties.
}
\references{
  Many of the gene sets, as well as their example use, are described in

  James K. T. Wang, Peter Langfelder, Steve Horvath and Michael J. Palazzolo. Exosomes and Homeostatic Synaptic Plasticity Are Linked to Each other and to Huntington's, Parkinson's, and Other Neurodegenerative Diseases by Database-Enabled Analyses of Comprehensively Curated Datasets. Front. Neurosci. 11 149, 2017; online 31 March 2017. https://doi.org/10.3389/fnins.2017.00149 
}
\author{
  Gene sets were originally collected by Michael Palazzolo and Jim Wang under contract with CHDI
Management/CHDI Foundation. The collection was created by Peter Langfelder.
}
\examples{
  \dontrun{
  collection = HDTargetDBCollection(organism = "human");
  knownGroups(collection)
  }

}
\keyword{misc}% __ONLY ONE__ keyword per line
