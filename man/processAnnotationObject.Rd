\name{processAnnotationObject}
\alias{processAnnotationObject}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simple processing of annotation objects.
}
\description{
  This function is intended to simplify processing annotation objects that follow the standard Bioconductor
structure of \code{AnnotationDBI} objects. It returns a named list (or vector, if possible) containing the
data in the annotation object. If a vector of reference identifiers is supplied, the output will match the
components specified by the reference identifiers.
}
\usage{
processAnnotationObject(object, uniqueOnly = TRUE, warnNonUnique = FALSE, reference = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  An annotation object. 
}
  \item{uniqueOnly}{
  Logical: should only components with a unique value be returned? If all components are unique (i.e., all
identifiers are mapped to a single value), the output will be simplified to a vector.
}
  \item{warnNonUnique}{
  Logical: should a warning be issued when some components are non-unique?
}
  \item{reference}{
  An optional vector of reference identifiers. If given, the components of the output will match the input
identifiers, with missing values for those reference identifiers that are not present in the annotation
object.
}
}
\value{
  If all identifiers map to a unique value, a named vector, otherwise a list.
}
\author{
  Peter Langfelder
}
\examples{
  require(org.Hs.eg.db);

  # Convert Entrez IDs to symbols (the last number is an invalid Entrez ID):
  entrez = c(1, 10, 284748374);
  symbols = processAnnotationObject(org.Hs.egSYMBOL, reference = entrez);

  symbols
  #        1        10 UnknownID 
  #   "A1BG"    "NAT2"        NA 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
