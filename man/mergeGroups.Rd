\name{mergeGroups}
\alias{mergeGroups}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Merge groups into one
}
\description{
Merge given groups into one.
}
\usage{
mergeGroups(collection, groupNames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{collection}{
A collection of gene sets and gene properties.
}
  \item{groupNames}{
A character vector with names of the groups to be merged into one. All groups in this vector will be merged
into a single group which will carry the name, description and other information of the first group in
\code{groupNames}.
}
}
\value{
A collection with the indicated groups merged.
}
\author{
Peter Langfelder
}
\keyword{misc}% __ONLY ONE__ keyword per line
