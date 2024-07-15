\name{knownGroups}
\alias{knownGroups}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Groups present in given collection.
}
\description{
  The function returns a character vector of names of all groups present in the given collection. The
groups are determined from the actual gene sets.
}
\usage{
knownGroups(collection,
            sortBy = c("alphabet", "size"), reverseSort = FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{collection}{
     The collection object in which to look for groups.
}
  \item{sortBy}{
     Character string spcifying how to sort the groups. Recognized choices are (unique abbreviations of) 
\code{"alphabet"} (for alphabetical sort) and \code{"size"} (for sorting by the number of gene sets that
belong to the group in decreasing order).
}
  \item{reverseSort}{
     Logical indicating whether the sort should be reversed from its default direction (which is increasing
for the alphabetical sort and decreasing for sorting by size).
}
}
\value{
  A character vector of names of all groups present in the given \code{collection}.
}
\author{
  Peter Langfelder
}

\seealso{
  \code{\link{newGroup}} for creating new groups, 

  \code{\link{addToCollection}} to 
}
\keyword{ misc }% __ONLY ONE__ keyword per line
