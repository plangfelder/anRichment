\name{geneSetsToHumanReadable}
\alias{geneSetsToHumanReadable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert gene sets in a collection to a "human readable" format
}
\description{
This function creates a data frame that represents the gene sets in a collection, with 4 components: set
name, set description, source, and gene symbols of the genes in the collection concatenated into a single
string.
}
\usage{
geneSetsToHumanReadable(collection, sep = ", ")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{collection}{
A collection of gene sets.
}
  \item{sep}{
Separator for gene symbols.
}
}
\value{
Data frame with one row per gene set and the following columns:
\item{setName}{The gene set name.}
\item{setDescription}{The description of the gene set.}
\item{source}{The source of the gene set (the \code{source} component of the gene set meta-data).}
\item{Gene symbols}{Gene symbols of the genes in the gene set, concatenated with \code{sep} as the
separator.}
}
\author{
Peter Langfelder
}
\keyword{misc}% __ONLY ONE__ keyword per line
