\name{subCollectionIndex}
\alias{subCollectionIndex}
\title{
Find the index of sub-collection data sets in a (larger) collection.
}
\description{
This function finds the index of sub-collection data sets in a (larger) collection.
}
\usage{
subCollectionIndex(
  subCollection, 
  combinedCollection, 
  combinedIDNames = NULL, 
  stopOnUnmatched = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subCollection}{
The sub-collection.
}
  \item{combinedCollection}{
The larger collection.
}
  \item{combinedIDNames}{
Optional; if given, must be a vector with one component per data set in the combined collection. Each
component must be a concatenation of ID and name of the data set. If this function is called multiple times,
pre-calculating this argument can save some time.
}
  \item{stopOnUnmatched}{
Logical: should the function throw an error if it cannot match a data set in the sub-collection to any of the
data sets in the combined collection?
}
}
\details{
Data sets are matched based on the concatenation of their ID and name.
}
\value{
A numeric vector giving the indices of the data sets in the sub-collection in the combined collection.
}
\author{
Peter Langfelder
}

\seealso{
\code{\link{dataSetIDs}} and \code{\link{dataSetNames}} for getting data set IDs and names, respectively.
}
\keyword{misc}% __ONLY ONE__ keyword per line
