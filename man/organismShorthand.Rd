\name{organismShorthand}
\alias{organismShorthand}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Organism shorthand
}
\description{
  This function returns the scientific shorthand of the given organism.
}
\usage{
organismShorthand(organism)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{organism}{
      Character string giving one of the common name, scientific name, or scientific shorthand.
}
}
\details{
  The recognized organisms are those for which there exists a Bioconductor organism annotation package
(usually named org.Xx.eg.db, where Xx is the scientific shorthand.)

See \code{\link{organismLabels}} for all recognized organisms.

}
\value{
  Character string giving the organism shorthands.
}
\author{
  Peter Langfelder
}
\seealso{
  \code{\link{organismLabels}}
}
\examples{
  organismShorthand("human")
  organismShorthand("mouse")
}
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
