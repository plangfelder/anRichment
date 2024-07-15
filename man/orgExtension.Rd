\name{orgExtension}
\alias{orgExtension}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Package extension for Bioconductor's org annotation packages.
}
\description{
  Utility function returning the package extension for Bioconductor's org annotation packages.
}
\usage{
orgExtension(organism)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{organism}{
    Character string identifying the organism. Can be a common name (e.g., "human"), scientific name
(e.g., "Homo sapiens"), or the scientific name shorthand (e.g., "Hs"). To see available organisms, use
\code{\link{organismLabels}} with the default arguments (i.e., 'organismLabels()').
}
}
\details{
  Bioconductor's annotation package names currently use the form org.Xx.eg.db where Xx is the shorthand of
the scientific name of the organism, except of yeast, for which the "eg" is replaced by "sgd".
}
\value{
  Character string, \code{"eg"} or \code{"sgd"}, appropriate for the given \code{organism}.
}
\author{
  Peter Langfelder
}
\note{This function relies on the current name convention of Bioconductor's packages which may change in
the future.}
\examples{
  orgExtension("human")
  orgExtension("yeast")
}
\keyword{ misc }
