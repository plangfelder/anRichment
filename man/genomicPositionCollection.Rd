\name{genomicPositionCollection}
\alias{genomicPositionCollection}

\title{
Build genomic position collection for anRichment
}
\description{
Using annotation data from Bioconductor, this function builds a collection whose gene sets "tile" the entire
genome, i.e., each gene set contains genes within a window around a specific position.
}
\usage{
genomicPositionCollection(
  organism, 
  spacings,
  namePrefix = "Genomic position", 
  dataSource = "",
  overlapFactor = 2,
  membershipBy= c("start", "end", "overlap"),
  useUnits = "Mb",
  unit = 1e6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{organism}{
Organism for which the collection should be creataed. At present only human and mouse are supported.
}
  \item{spacings}{A numeric vector giving the window sizes.
}
  \item{namePrefix}{Character string giving the prefix for gene set IDs.  
}
  \item{dataSource}{Character string giving the source of gene sets.
}
  \item{overlapFactor}{A single number. Windows will be moved by spacing divided by \code{overlapFactor}.
}
  \item{membershipBy}{One of "start", "end", "overlap". Determines which exact point on the gene will be used
to determine which gene set the gene belongs to.
}
  \item{useUnits}{Character string representing the unit of genomic distance to use in gene set names.
}
  \item{unit}{Length of the unit.
}
}
\details{
Each chromosome is processed individually. 
}
\value{
A collection of gene sets.
}
\author{
Peter Langfelder
}

\seealso{
  \code{\link{allCollections}} for loading all collections available in anRichment.
}
\keyword{misc}
