\name{mapCollectionIDs}
\alias{mapCollectionIDs}
\title{
Map or translate gene identifiers in a collection.
}
\description{
This function maps gene identifiers in a collection using a user-supplied map. The map must map each new identifier to an
old identifier; several new identifiers may be mapped to one old identifier and they will appear in each data set that the
old identifier appeared in, with the same evidence codes.
}
\usage{
mapCollectionIDs(collection, new2old, keepOrder = TRUE)
}
\arguments{
  \item{collection}{
  A collection (see \code{link{newCollection}}).
}
  \item{new2old}{
  A two-column matrix or data frame in which column 1 are new IDs and column 2 the corresponding old IDs. 
Each new identifier may appear only once
but several new identifiers may map to a single old identifiers. Old identifiers that do not appear in the second column
are discarded.
}
  \item{keepOrder}{
  Logical: should the new identifiers be ordered in the order of the corresponding old identifiers? 
}
}
\value{
  A collection with mapped identifiers.
}
\author{
  Peter Langfelder
}

\seealso{
\code{\link{newCollection}} for creating collections;

\code{\link{convertCollectionToOrganism}} for mapping collections between organisms.
}
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{misc}
% \concept{enrichment analysis}
% Use only one concept per line.
