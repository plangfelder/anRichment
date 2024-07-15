\name{sameOrganism}
\alias{sameOrganism}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine whether two strings describe the same organism
}
\description{
Determine whether two strings describe the same organism. For example, \code{"human"}, \code{"Hs"} and
\code{"Homo Sapiens"} all refer to the same organism. Organisms can be specified
by a common name, scientific name, or shorthand. See \code{\link{organismLabels}} for known organisms.
}
\usage{
sameOrganism(org1, org2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{org1}{
 Character string describing an organism.
}
  \item{org2}{
 Character string describing an organism.
}
}
\value{
  Logical indicating whether the two input strings denote the same organism.
}
\author{
  Peter Langfelder
}
\seealso{
  \code{\link{organismLabels}} for a list of known organisms and the various ways of identifying them.
}
\examples{
  sameOrganism("Hs", "Human");
  sameOrganism("Mm", "Human")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }% __ONLY ONE__ keyword per line
