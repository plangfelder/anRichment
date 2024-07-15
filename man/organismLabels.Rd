\name{organismLabels}
\alias{organismLabels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Organism labels.
}
\description{
  Returns 3 organism labels (common name, scientific name, and scientific shorthand) for the given
organism.
}
\usage{
organismLabels(organism = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{organism}{
    Character string giving one of the common name, scientific name, or scientific shorthand. 
}
}
\value{
  The recognized organisms are those for which there exists a Bioconductor organism annotation package
(usually named org.Xx.eg.db, where Xx is the scientific shorthand.)

  If the input \code{organism} is non-\code{NULL} and recognized, its 3 labels (common name, scientific
name, and scientific shorthand) are returned. If the input \code{organism} is \code{NULL}, the function
returns a data frame with columns giving the common names, scientific names, and scientific shorthands for
all recognized organisms.
}
\author{
  Peter Langfelder
}
\examples{
organismLabels()

organismLabels("human")
}
\keyword{ misc }% __ONLY ONE__ keyword per line
