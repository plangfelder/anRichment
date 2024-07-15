\name{allOrgGenes}
\alias{allOrgGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Entrez identifiers of all genes for an organism.
}
\description{
  This function utilizes the Bioconductor's annotation packages org.*.eg.db to extract a list of all entrez
identifiers known for the given organism. 
}
\usage{
allOrgGenes(organism)
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
This function should be used with caution since the returned list of all Entrez identifiers may not quite
be complete.  
}
\value{
  A vector of Entrez identifiers.
}
\author{
  Peter Langfelder
}
\seealso{
  Bioconductor annotation packages org.*.eg.db
}
\examples{
  entrez = allOrgGenes("human")
  length(entrez)
  entrez[1:10]
}
\keyword{misc }
