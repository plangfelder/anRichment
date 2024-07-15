\name{geneAnnotationFromEntrez}
\alias{geneAnnotationFromEntrez}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create a rudimentary gene annotation table from Entrez IDs
}
\description{
 This function is basically a wrapper for Bioconductor's org.Xx.eg.db packages and creates a simple
3-column annotation data frame containing the Entrez ID, gene symbol, and gene name.
}
\usage{
geneAnnotationFromEntrez(entrez, organism, includePosition = FALSE)
}
\arguments{
  \item{entrez}{
 A vector of Entrez identifiers.
}
  \item{organism}{
 Character string specifying the organism to which the Entrez identifiers belong. Organism can be specified
by a common name, scientific name, or shorthand. See \code{\link{organismLabels}} for known organisms.
}
  \item{includePosition}{Logical: should chromosome and chromosomal position be included in the output?}
}
\value{
  A data frame with columns giving the gene Entrez ID, official symbol and name, and optionally chromosome
and base pair location (position).
}
\author{
  Peter Langfelder
}
\note{
  The requisite org.Xx.eg.db package must be installed; if it is not, the function will throw an error.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

 \code{org.Xx.eg.db} packages from Bioconductor for the actual data;

 \code{\link{organismLabels}} for a list of known organisms and the various ways of identifying them.

}
\keyword{misc}% __ONLY ONE__ keyword per line
