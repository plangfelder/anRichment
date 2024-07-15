\name{newGeneSet}
\alias{newGeneSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create a new gene set
}
\description{
  This function creates a new gene set.
}
\usage{
newGeneSet(
   # Gene information
   geneEntrez, geneEvidence, geneSource, 
   # Set information
   ID, name, shortName = name, description, 
   source, organism,
   internalClassification, groups,
   lastModified = Sys.Date(),
   format = "\%Y-\%m-\%d",
   alternateNames = character(0),
   externalDB = "",
   externalAccession = "",
   webLink = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneEntrez}{
     A vector of Entrez identifiers of the genes to be included in this gene set. A gene can be included
multiple times (for example, with different evidence codes, see below).
}
  \item{geneEvidence}{
     Character vector giving the evidence codes of the genes. Must have the same length as
\code{geneEntrez}. For available evidence codes, see \code{\link{knownEvidenceCodes}}.
}
  \item{geneSource}{
     Character vector giving the source (reference, web site, etc) for each gene/evidence. Must have the
same length as \code{geneEntrez}.  
}
  \item{ID}{
     Identifier for the gene set. This should be a sequence of characters (or numbers) that is unique for
this gene set.
}
  \item{name}{
     Descriptive short name for the gene set.
}
  \item{shortName}{
     An optional short version of the gene set name, used for printing when space is at premium.
}
  \item{description}{
     Description of the gene set.
}
  \item{source}{
     Literature or other reference for the gene set. This can be the same as (one of) the \code{geneSource}
references above, but could also be different.
}
  \item{organism}{
     Character string identifying the organism for which the gene set is defined. Can be a
common name (e.g., "human"), scientific name
(e.g., "Homo sapiens"), or the scientific name shorthand (e.g., "Hs"). To see available organisms, use
\code{\link{organismLabels}} with the default arguments (i.e., 'organismLabels()').

}
  \item{internalClassification}{
     A character vector of arbitrary length that the user can use to categorize the gene set.
}
  \item{groups}{
     Character vecotor giving the groups the gene set belongs to.
}
  \item{lastModified}{Date of last modification of the gene set, given 
either as a \code{\link{Date}} object, or as a character string specifying date. If a character string, its
format can be specified using the argument \code{format} below.
}
   \item{format}{Specification of format of character dates. See \code{\link{strptime}} for the
description of the format string; the default means <year with century>-<month>-<day>. }

  \item{alternateNames}{Optional alternate names for the gene set. }
  \item{externalDB}{Character string giving the name of the external database from which the gene set was sourced. }
  \item{externalAccession}{Character string giving the accession identifier in the external database 
above that corresponds to the gene set.}
  \item{webLink}{Character string giving the web link (URL) to the external database (or directly to the gene
set in the external database).}

}

\details{
  The function does not try to ensure validity of entries in the \code{groups} vector by comparing
\code{groups} to some known list of groups. The user should beware of misspellings.
}
\value{
  A new gene set, which is a list with components encoding the input information. The list can be accessed
and modified by the user, but the user is explicitly discouraged from writing code that relies on the
structure on the list since it may change in the future.
}
\author{
  Peter Langfelder
}
\seealso{
  \code{\link{newGroup}} for creating new groups,

  \code{\link{newCollection}} and \code{\link{addToCollection}} for creating collections and adding gene
sets to them.
}
\keyword{ misc }% __ONLY ONE__ keyword per line
