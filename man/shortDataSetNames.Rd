\name{shortDataSetNames}
\alias{shortDataSetNames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Short names of gene sets in a collection
}
\description{
Returns the short gene set names in the given collection, optionally restricted to only
the given groups or gene sets. Supplying gene set IDs in argument \code{dataSets} will convert them to short gene set names.
}
\usage{
shortDataSetNames(collection, groups = NULL, dataSets = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{collection}{
   A collection of gene sets.
}
  \item{groups}{
     Names of the groups to be considered. The default \code{NULL} is equivalent to considering all groups.
}
  \item{dataSets}{
     Names or IDs of gene sets to be considered. The default \code{NULL} is equivalent to considering all
gene sets. 
}

}
\value{
  A vector of short gene set names
}
\author{
  Peter Langfelder
}

\seealso{

  \code{\link{dataSetIDs}} for obtaining gene set IDs;
  \code{\link{dataSetNames}} for obtaining (full) gene set names;

  \code{\link{newGeneSet}} for creating new gene sets;
  \code{\link{newGeneProperty}} for creating new gene properties (quantitative or ordinal characteristics
for genes);
}
\keyword{ misc }% __ONLY ONE__ keyword per line
