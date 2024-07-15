\name{dataSetMatrix}
\alias{dataSetMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a matrix representation of the collection
}
\description{
Create a matrix representation of the collection, i.e., a matrix in which every column corresponds to a data
set, and each row to a gene that is present the background (which could be all genes present in the
collection). For gene sets, the matrix contains a one
if the gene belongs to the gene set, and zero otherwise.
}
\usage{
dataSetMatrix(
   collection, 
   background = NULL, 
   includeProperties = TRUE,
   includeGeneSets = TRUE,
   dropZeroOverlapSets = FALSE,
   namesFrom = c("ID", "name"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{collection}{
A collection of gene sets and/or gene properties.
}
  \item{background}{
Optional specification of background, i.e., the genes to be represented in the matrix. If not given, all
genes in the collection will be used as background.
}
  \item{includeProperties}{
Logical: should gene properties be included in the matrix?
}
  \item{includeGeneSets}{
Logical: should gene sets be included in the matrix?
}
  \item{dropZeroOverlapSets}{
Logical: should gene sets with zero overlap with the specified background be dropped from the matrix?
}
  \item{namesFrom}{
Determines the component of the gene set or property meta data that will be used as column names of the
result.
}
}
\value{
A matrix in which rows correspond to genes and columns to data sets (gene sets or properties). For gene sets,
membership is encoded by 1 (and non-membership by 0); for properties, the property values are copied.
}
\author{
Peter Langfelder
}
\keyword{misc}% __ONLY ONE__ keyword per line
