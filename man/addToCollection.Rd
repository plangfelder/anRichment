\name{addToCollection}
\alias{addToCollection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add gene sets or groups to an existing collection.
}
\description{
Add gene sets or groups to an existing collection.
}
\usage{
addToCollection(collection, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{collection}{
Collection to add to.
}
  \item{\dots}{
One or more gene sets or groups. The type of each object is checked and determined automatically.
}
}
\value{
The input collection with the gene sets and groups added.
}
\author{
Peter Langfelder
}

\seealso{
\code{\link{newGeneSet}}, \code{\link{newCollection}}
}
\keyword{ misc}% __ONLY ONE__ keyword per line
