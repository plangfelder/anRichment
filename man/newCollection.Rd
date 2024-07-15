\name{newCollection}
\alias{newCollection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create a new collection.
}
\description{
  A collection is a structure that contains gene sets and some associated information. This function sets up
the required structure and fills it with supplied information.
}
\usage{
newCollection(dataSets = list(), groups = list(), 
              identifiers = numeric(0),
              weights = matrix(nrow = 0, ncol = 1),
              ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSets}{
     Gene sets to be stored in the collection. A list where each component corresponds to a gene set.
}
  \item{groups}{
     Information about groups the gene sets belong to. A list where each component corresponds to a
group.
}
  \item{identifiers}{A vector of gene identifiers (usually Entrez IDs). Users should normally not have to
use this argument.}

  \item{weights}{A matrix of weights corresponding to the genes in \code{identifiers}. As with
\code{identifiers}, users should normally leave this argument at the default value.}

  \item{\dots}{
     Additional gene sets or groups to add to the collection.
}
}
\value{
  An collection containing the gene sets and the groups specified.
}
\author{
  Peter Langfelder
}
\seealso{
 \code{\link{newGeneSet}} to create new gene sets;

  \code{\link{newGeneProperty}} for creating new gene properties (quantitative or ordinal characteristics 
for genes);

 \code{\link{newGroup}} to create new groups;

 \code{\link{addToCollection}} to add additional gene set(s) and/or group(s) to an existing collection.

}
\keyword{ misc}% __ONLY ONE__ keyword per line
