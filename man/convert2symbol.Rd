\name{convert2symbol}
\alias{convert2symbol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Convert Entrez identifiers to gene symbols.
}
\description{
  Convert Entrez identifiers to gene symbols using Bioconductor annotation packages.
}
\usage{
convert2symbol(entrez, organism)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{entrez}{
     A (numeric) vector of Entrez identifiers, or a list of such vectors.
}
  \item{organism}{
     Character string identifying the organism. Can be a common name (e.g., "human"), scientific name
(e.g., "Homo sapiens"), or the scientific name shorthand (e.g., "Hs"). To see available organisms, use
\code{\link{organismLabels}} with the default arguments (i.e., 'organismLabels()').
}
}
\value{
  A vector of gene symbols, or a lits of such vectors (depending on whether the input 'entrez' was a vector or a list). 
  Unrecognized Entrez identifiers results in a NA (missing) output.
}
\author{
  Peter Langfelder
}
\seealso{
  \code{\link{convert2entrez}} for the reverse mapping.
}
\examples{
  entrez = c(348, 54674, 123456789);
  convert2symbol(entrez = entrez, organism = "human")
}
\keyword{ misc }% __ONLY ONE__ keyword per line
