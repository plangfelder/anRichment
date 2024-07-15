\name{collectionFromGeneLists}
\alias{collectionFromGeneLists}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a rudimentary gene set collection
}
\description{
This function creates a rudimentary reference collection from given gene Entrez identifiers and the
corresponding gene set names. This can be useful for quick tests of enrichment; 
for long-term use we strongly recommend
using the full functionality for creating annotated collections, for example
\code{\link{collectionFromDataFrames}}.
}
\usage{
collectionFromGeneLists(
   identifiers, 
   setNames = NULL, 
   organism, 
   userName = "user",
   sortSets = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{identifiers}{
Gene identifiers. This can be either a named list in which each component corresponds to a gene set,
or a simple vector, in which case gene set membership can be either given in \code{setNames} below or by
the \code{names} attribute.
}
  \item{setNames}{
A vector specifying the set name to which each gene in \code{entrez} belongs. Not used if \code{entrez} is
a list or named vector.
}
  \item{organism}{
Character string identifying the organism. This can be a common name, scientific name, or
scientific shorthand. See \code{\link{organismLabels}} or simply type \code{organismLabels()} to see all
available organism labels.}
  \item{userName}{
A character string that will be used as a fill for the various information fields in gene sets, groups and
collection. 
}
  \item{sortSets}{
Logical: should the sets be sorted (by name, in increasing order) in the final collection?
}
}
\details{
The function creates a collection (see \code{\link{newCollection}}) containig the gene sets specified in
\code{entrez} and \code{setNames}, and a single group whose name equals \code{userName}. 
}
\value{
A collection containing specified gene sets.
}
\author{
Peter Langfelder
}
\seealso{
\code{\link{newGeneSet}}, \code{\link{newGroup}}, \code{\link{newCollection}} for creating new gene sets,
groups and collections programmatically;

\code{\link{collectionFromDataFrames}} for creating a collection from appropriately formatted data frames
or text tables. 
}
\keyword{misc}% __ONLY ONE__ keyword per line
