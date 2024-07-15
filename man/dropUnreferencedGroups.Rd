\name{dropUnreferencedGroups}
\alias{dropUnreferencedGroups}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Drop unreferenced groups from a collection.
}
\description{
This function drops those groups from the collections that do not tag any data sets.
}
\usage{
dropUnreferencedGroups(collection, ignore.case = TRUE, verbose = 1, indent = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{collection}{
A collection of data sets.
}
  \item{ignore.case}{Logical: should case (uppercase vs. lowercase) be ignored when matching patterns? }

  \item{verbose}{ integer level of verbosity. Zero means silent, higher values make the output
progressively more and more verbose. }

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }

}
\details{
Unreferened groups are groups that exist in the \code{group} component of a collection but none of the gene
sets are tagged by these groups or by their "children" or sub-groups.
}
\value{
The input collection with the unreferenced groups removed.
}
\author{
Peter Langfelder
}
\keyword{misc}
