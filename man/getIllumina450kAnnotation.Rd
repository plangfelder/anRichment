\name{getIllumina450kAnnotation}
\alias{getIllumina450kAnnotation}
\title{
Retrieve a shortened form of the Illumina 450k manifest (probe annotation table).
}
\description{
This function returs a shortened form of the Illumina 450k manifest (probe annotation table). The current
manifest is version 1.2, downloaded 2016/02/20. The shortened version omits the sequence information columns.
}
\usage{
getIllumina450kAnnotation(addEntrez = TRUE)
}
\arguments{
  \item{addEntrez}{
Logical: add a column with Entrez IDs?
}
}
\value{
A data frame containing the 450k probe information.
}
\references{
Data provided by Illumina at
https://support.illumina.com/downloads/infinium_humanmethylation450_product_files.html .
}
\author{
Data are provided by Illumina; R code by Peter Langfelder.
}
\keyword{misc}% __ONLY ONE__ keyword per line
