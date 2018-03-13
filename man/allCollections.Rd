\name{allCollections}
\alias{allCollections}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
allCollections(
  merge = FALSE,
  organism = "human", 
  useHomology = TRUE,
  addOldOrganismToSetNames = FALSE, 
  namePattern = ".convertedFrom.%o",
  addOldOrganismToSetDescriptions = FALSE,
  descriptionPattern = " (Converted from %o.)",

  trimInternal = TRUE,
  traceableHDTargetDB = TRUE,

  buildExternal = FALSE,
  genomicSpacings = 5e6,

  MSigDBxml)
  
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{merge}{Logical: should the various collections be merged into a single collection (\code{TRUE}), 
or returned as a list with each component corresponding to a separate collection (\code{FALSE})? }
  \item{organism}{
   Character string specifying the organism. Can be common name (such as "human"), scientific name ("Homo
sapiens"), or scientific shorthand ("Hs"). To see all recognized organisms, use function
\code{\link{organismLabels}}.  If \code{NULL} or \code{NA}, gene sets in the collection will retain their
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

  \item{trimInternal}{Logical: should sets that are also contained in \code{\link{HDSigDBCollection}} be removed?}

  \item{traceableHDTargetDB}{Logical: should only a subset of the gene sets that are traceable to a published
article or website be included?}

  \item{buildExternal}{Logical: should external collections (namely GO, MSigDB and genomic position
collection) be built? Building of MSigDB can be suppressed by setting \code{MSigDBxml = NULL}, and building
of genomic position collection can be suppressed by setting \code{genomicSpacings = NULL}.) }

  \item{genomicSpacings}{A vector of spacings around which to build the genomic position collection. Each
spacing is used as a sliding window, moved by half the spacing between two sets.}

  \item{MSigDBxml}{Path to an MSigDB file in XML format. If \code{NULL}, building of MSigDB is skipped.}
}
\details{
The function loads and returns all internally stored collections: (1) the "Yang literature collection" (see
\code{\link{YangLiteratureCollection}}), (2) HDSigDB collection (\code{\link{HDSigDBCollection}}),
(3) HD Target DB (\code{\link{HDTargetDBCollection}}), (4) "internal" collection largely inherited from older
versions and originally compiled by Jeremy Miller for the \code{userListEnrichment} function in WGCNA
(\code{\link{internalCollection}}), (5)
collection of gene sets collected by Jeremy Miller at AIBS (\code{\link{MillerAIBSCollection}}), and (6)
collection of pathways from NCBI BioSystems (\code{\link{BioSystemsCollection}}).

Additionally, the folliwng external collections are optionally built and returned as well: (1) Gene Ontology
(GO) collection, based on annotation packages from the Bioconductor project (see
\code{\link{buildGOcollection}}); (2) Molecular Signatures Database (MSigDB, see
\code{\link{buildMSigDBcollection}}); and (3) genomic position collection
(\code{\link{genomicPositionCollection}}).

}
\value{
Either a merged collection, or a list in which each component is a collection.
}
\author{
Code by Peter Langfelder.
}
\seealso{
Functions that return individual collections: \code{\link{YangLiteratureCollection}},
\code{\link{HDSigDBCollection}}, \code{\link{HDTargetDBCollection}}, \code{\link{internalCollection}},
\code{\link{MillerAIBSCollection}}, \code{\link{BioSystemsCollection}}, \code{\link{buildGOcollection}}, 
\code{\link{buildMSigDBcollection}}, \code{\link{genomicPositionCollection}}.
}
\keyword{ ~kwd2 }
