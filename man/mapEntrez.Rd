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
  \item{version}{Caution, experimental. Version to be used; a real number of the form year.month. Currently, 3 thresholds
are recognized. Versions below 2014.04 will use a human-mouse homology; versions between 2014.04 and 2019.03 will use
multi-organism homology downloaded in April 2014, versions between 2019.03 and 2025.12 will use multi-organism homology
downloaded in March 2019, and versions higher or equal to 2025.12 as well as the default \code{NULL} will use
the multi-organism homology information downloaded on December 2, 2025.}
}
\details{
The current homology file at was downloaded from the MGI web site

https://www.informatics.jax.org/homology.shtml

on December 2, 2025. The file download link is currently (December 2025)

https://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt

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
