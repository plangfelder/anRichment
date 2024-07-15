\name{collectionFromClasses}
\alias{collectionFromClasses}
\title{
Create a collection in with gene sets that correspond to given classes (e.g., gene modules).
}
\description{
Create a collection in with gene sets that correspond to given classes (e.g., gene modules).
}
\usage{
collectionFromClasses(
  identifiers, 
  labels,
  organism,
  shortAnalysisName,

  classDescription,
  classCol = "class",
  nameCol,
  shortNameCol,
  descriptionCol,
  ignoreClasses = 0,
  evidence = "other",
  source = shortAnalysisName,
  IDBase, IDBaseSep = ".",
  groups,
  internalClassification,
  lastModified = Sys.Date(),
  format = "\%Y-\%m-\%d")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{identifiers}{
Vector of indetifiers (usually Entrez) for the genes.
}
  \item{labels}{
Class labels in the same order as \code{identifiers} above.
}
  \item{organism}{
Organism to which the gene IDs correspond to.
}
  \item{shortAnalysisName}{
A short name for the analysis in which the classes were determined.
}
  \item{classDescription}{
A data frame containing information about the classes present in \code{labels}. 
}
  \item{classCol}{
Name of the column that contains the class ID, that is the value of \code{labels} that corresponds to the
class.
}
  \item{nameCol}{
Name of the column in \code{classDescription} that contains human-readable names for the classes.
}
  \item{shortNameCol}{
Name of the column in \code{classDescription} that contains a short version of the human-readable names for
the classes.
}
  \item{descriptionCol}{
Name of the column in \code{classDescription} that contains a more detailed class description.
}
  \item{ignoreClasses}{
Optional vector of labels of "uninteresting" classes (e.g., labeling genes that do not belong to any proper
class).
}
  \item{evidence}{
A character string giving the evindence code to be used for the gene content. The same code will be used for
all genes. For all recognized evidence codes, see \code{\link{knownEvidenceCodes}}.
}
  \item{source}{
Character string to be used as the "source" field of the generated gene sets.
}
  \item{IDBase}{
Character string that will serve as the base of the gene set IDs. \code{IDBaseSep} and 
numbers starting from 1 will be added to the base to form the actual set IDs.
}
  \item{IDBaseSep}{
Character string to separate \code{IDBase} and numbers in the gene set IDs.
}
  \item{groups}{
A list in which every component is a "group" (see \code{\link{newGroup}}) that the gene sets will belog to.
All gene sets will belong to the same set of groups.
}
  \item{internalClassification}{
A character vector giving the "internal classification" for the gene sets.
}
  \item{lastModified}{
Date of last modification to be recorded with the gene sets. 
}
  \item{format}{
Date format.
}
}
\value{
A collection of gene sets (see \code{\link{newCollection}}) containing the classes as gene sets.
}
\author{
Peter Langfelder
}

\seealso{
\code{\link{newGroup}}, \code{\link{newGeneSet}} and \code{\link{newCollection}} for descrition of gene
sets, groups and collections;
\code{\link{collectionFromGeneLists}} for a simpler way to create a collection from gene lists.
}
\keyword{misc}% __ONLY ONE__ keyword per line
