\name{buildMSigDBCollection}
\alias{buildMSigDBCollection}
\alias{MSigDBCollection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Returns a collection containing Molecular Signatura Database (MSigDB) gene sets.
}
\description{
  Starting from an XML file, these functions return collections containing Molecular Signatura Database (MSigDB) gene sets.
}
\usage{
buildMSigDBCollection(
  file, 
  MSDBVersion = "5.0", 
  excludeCategories = c("C2", "C5"))

MSigDBCollection(
  file, 
  MSDBVersion = "5.0", 
  organism = "human", 
  useHomology = TRUE,
  addOldOrganismToSetNames = FALSE, 
  namePattern = ".convertedFrom.%o",
  addOldOrganismToSetDescriptions = FALSE,
  descriptionPattern = " (Converted from %o.)",
  verbose = 1, indent = 0, 
  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{XML file (including path) such as the one distributed by Broad Institute.}

  \item{MSDBVersion}{Character string containing the MSigDB version. }

  \item{excludeCategories}{Character vector containing the MSigDB category labels to be excluded. The default
corresponds to "curated gene sets  from online pathway databases, publications in PubMed, and knowledge of
domain experts" (C2) and GO (C5). }
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

  \item{verbose}{Level of verbosity. Zero means silent, higher numbers mean more verbose output.}

  \item{ident}{Indentation of the diagnostic messages. Each unit adds 2 spaces.}

  \item{...}{Other arguments, currently ignored.}
}
\details{
These functions convert the XML file representation of MSigDB gene sets in a collection usable within the
anRichment package. The XML file can be downloaded from http://software.broadinstitute.org/gsea/msigdb/.
}
\value{
  Collection of MSigDB gene lists. \code{buildMSigDBCollection} returns gene sets as they appear in MSigDB;
function \code{MSigDBCollection} also converts the gene sets to a target organism.
}
\author{
 MSigDB is maintained by the GSEA team at Broad Institute, http://software.broadinstitute.org/gsea/msigdb. R
code by Peter Langfelder.
}
\keyword{misc}% __ONLY ONE__ keyword per line
