\name{mergeCollections}
\alias{mergeCollections}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Merge collections into one.
}
\description{
  Given multiple collections, this function merges them into a single structure. If groups with the same name
exist in several collections, they will be merged. Data sets with duplicate names will be removed;
optionally, the function can throw an error on duplicate data set names.
}
\usage{
mergeCollections(..., stopOnDuplicates = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
     Two or more collections.
}
  \item{stopOnDuplicates}{Logical: should presence of two data sets with same name cause an error to be
thrown?}
}
\value{
  The merged collection.
}
\author{
  Peter Langfelder
}

\seealso{
 \code{\link{newCollection}} to create a new collection.
}
\keyword{ misc }% __ONLY ONE__ keyword per line
