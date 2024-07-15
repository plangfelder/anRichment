\name{renameMultipleGroups}
\alias{renameMultipleGroups}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rename multiple groups in one call.
}
\description{
Rename multiple groups in one call.
}
\usage{
renameMultipleGroups(collection, old2newTable)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{collection}{
A collection of data sets.
}
  \item{old2newTable}{
A matrix or data frame with two columns giving the old and new group names, respectively.
}
}
\value{
The input collection with groups renamed according to the translation table.
}
\author{
Peter Langfelder
}
\keyword{misc}
