\name{orgObject}
\alias{orgObject}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Specified object from Bioconductor's organism packages
}
\description{
  This function returns the specifed object from the appropriate Bioconductor org package.
}
\usage{
orgObject(organism, 
          objectName, 
          convertToList = TRUE, 
          uniqueOnly = FALSE, 
          warnNonUnique = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{organism}{
      Character string specifying the organism for which to build the collection. Can be a common name
(e.g.,
"human"), scientific name
(e.g., "Homo sapiens"), or the scientific name shorthand (e.g., "Hs"). To see available organisms, use
\code{\link{organismLabels}} with the default arguments (i.e., 'organismLabels()').
     
}
  \item{objectName}{
    Name of the object to return.
}
  \item{convertToList}{
    Logical: should the (database) object be converted to list?
}
  \item{uniqueOnly}{
     Logical: should only the list elements with length 1 be kept?
}
  \item{warnNonUnique}{
     Logical: should the function issue a warning when there are list elements with length 0 or more than
1?
}
}
\value{
  The org object, possibly converted to a list.
}
\author{
  Peter Langfelder
}
\keyword{misc }% __ONLY ONE__ keyword per line
