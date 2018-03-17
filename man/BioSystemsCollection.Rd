\name{BioSystemsCollection}
\alias{BioSystemsCollection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Returns a collection containing NCBI Biosystems gene sets.
}
\description{
  Returns a collection containing NCBI Biosystems gene sets.
}
\usage{
BioSystemsCollection( organism = "human", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{organism}{
   Character string specifying the organism. Can be common name (such as "human"), scientific name ("Homo
sapiens"), or scientific shorthand ("Hs"). To see all recognized organisms, use function
\code{\link{organismLabels}}. 
}
  \item{...}{Ignored.}
}
\details{
BioSystems data were downloaded on 2017-10-20.
}
\value{
  Collection of gene sets corresponding to NCBI BioSystems entities.
}
\author{
Code by Peter Langfelder, gene sets from NCBI BioSystems.
}
\examples{
  collection = BioSystemsCollection(organism = "mouse");
  knownGroups(collection)
}
\keyword{misc}% __ONLY ONE__ keyword per line
