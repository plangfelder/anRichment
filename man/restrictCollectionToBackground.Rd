\name{restrictCollectionToBackground}
\alias{restrictCollectionToBackground}
\title{
Restrict A Collection To A Given Background Set Of Genes
}
\description{
Restrict genes in each data set to only those that are present in the given background.
}
\usage{
restrictCollectionToBackground(collection, background)
}
\arguments{
  \item{collection}{
A collection.
}
  \item{background}{
A vector of identifiers.
}
}
\value{
  The same collection with gene sets restricted to the background. This function does not check whether the restricted
sets actually contain any genes.
}
\author{
Peter Langfelder
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{misc}
