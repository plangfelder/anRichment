\name{mapEntrez}
\alias{mapEntrez}

\title{
Map entrez identifiers between organisms
}
\description{
Attempts to map Entrez gene identifiers between organisms, either using a suplied map, homology information (if
available), or by mapping gene symbols.
}
\usage{
mapEntrez(entrez, 
          orgFrom, orgTo, 
          entrezMap = NULL,
          useHomology = TRUE, version = NULL)
}
\arguments{
  \item{entrez}{
A vector of Entrez gene identifiers.
}
  \item{orgFrom}{
Organism to which the input Entrez identifiers correspond. This can be a scientific name (such as Homo
Sapiens), common name (human), or scientific shorthand (Hs). See \code{\link{organismLabels}} for
recognized organisms.
}
  \item{orgTo}{
Organism to which the input Entrez identifiers should be translated. 

}
    \item{entrezMap}{An optional 2-column matrix or data frame giving the mapping between identifiers for
\code{orgFrom} and \code{orgTo} in columns 1 and 2, respectively. If not given, the map will be generated as
described in Details. } 

  \item{useHomology}{
Logical: should homology information, if available, be used? If this is \code{FALSE} or homology is not
available, genes will be matched by gene symbols.
}
  \item{version}{Version to be used. Experimental.}
}
\details{
The homology files were downloaded from the MGI web site

http://www.informatics.jax.org/homology.shtml

in March 2019.

}
\value{
A vector of translated Entrez identifiers. Identifiers that could not be translated are turned into
\code{NA}.}

\author{
Peter Langfelder
}

\seealso{
\code{\link{entrezHomologs}} for creating suitable Entrez mappings;

\code{\link{organismLabels}}
}
\keyword{misc}% __ONLY ONE__ keyword per line
