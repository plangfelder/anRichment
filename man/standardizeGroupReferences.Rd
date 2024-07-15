\name{standardizeGroupReferences}
\alias{standardizeGroupReferences}
\title{
Change group references from alternate to standard names
}
\description{
The function changes the \code{groups} component of the data sets in the given collection from alternate to standard
names.
}
\usage{
standardizeGroupReferences(collection)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{collection}{
An anRichment collection.
}
}
\value{
The input collection, with the group references in the \code{group} components of the data sets changed from alternate to
standard names.
}
\author{
Peter Langfelder
}
\seealso{
\code{\link{newGroup}} for creating new groups with standard and alternate names.
}
\keyword{misc}% __ONLY ONE__ keyword per line
