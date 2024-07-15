\name{entrezHomologs}
\alias{entrezHomologs}

\title{
Creates a translation table of Entrez IDs between two organisms.
}
\description{
This function creates a translation table of Entrez IDs between the given organisms, either using external homology
information or by matching gene symbols (a crude approximation used when homology information is not
available).
}
\usage{
entrezHomologs(orgFrom, orgTo, useHomology = TRUE, version = NULL)
}
\arguments{
  \item{orgFrom}{
Organism to which the input Entrez identifiers correspond. This can be a scientific name (such as Homo
Sapiens), common name (human), or scientific shorthand (Hs). See \code{\link{organismLabels}} for
recognized organisms.
}
  \item{orgTo}{
Organism to which the input Entrez identifiers should be translated.
}

  \item{useHomology}{
Logical: should homology information, if available, be used? If this is \code{FALSE} or homology is not
available, genes will be matched by gene symbols.
}
  \item{version}{Caution, experimental. Version to be used; a real number of the form year.month. Currently, 2 thresholds
are recognized. Versions below 2014.04 will use a human-mouse homology; versions between 2014.04 and 2019.03 will use
multi-organism homology downloaded in April 2014, and versions above 2019.03 (as well as the default \code{NULL}) will use
the multi-organism homology information download in March 2019.}
}
\details{
The homology file at was downloaded from the MGI web site

http://www.informatics.jax.org/homology.shtml

on March 27, 2019. The file download link is currently (March 2019) 

http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt

}
\value{
A 2-column matrix giving pairs of corresponding Entrez IDs in \code{orgFrom} (column 1) and \code{orgTo}
(column 2).
}
\author{
Peter Langfelder
}
\seealso{
\code{\link{mapEntrez}}
}
\keyword{misc}
