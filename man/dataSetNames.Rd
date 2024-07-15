\name{dataSetNames}
\alias{dataSetNames}
\alias{dataSetIDs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Names and IDs of gene sets in a collection
}
\description{
These functions return the gene set names and IDs in the given collection, optionally restricted to only
the given groups or gene sets. Supplying gene set names to \code{dataSetIDs} will effectively convert the
names to IDs, while supplying IDs to dataSetNames will convert them to gene set names.
}
\usage{
dataSetNames(collection, groups = NULL, dataSets = NULL)
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
  A vector of gene set names or IDs.
}
\author{
  Peter Langfelder
}

\seealso{

  \code{\link{dataSetIDs}} for obtaining gene set IDs;

  \code{\link{newGeneSet}} for creating new gene sets;
  \code{\link{newGeneProperty}} for creating new gene properties (quantitative or ordinal characteristics
for genes);
}
\keyword{ misc }% __ONLY ONE__ keyword per line
