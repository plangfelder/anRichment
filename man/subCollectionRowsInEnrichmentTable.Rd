\name{subCollectionRowsInEnrichmentTable}
\alias{subCollectionRowsInEnrichmentTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finds index of rows in an enrichment table that correspond to data sets from a given sub-collection.
}
\description{
Finds index of rows in an enrichment table (typically returned by \code{enrichmentAnalysis})
that correspond to data sets from a given sub-collection of the collection that was used in the enrichment
analysis calculation.

}
\usage{
subCollectionRowsInEnrichmentTable(enrichmentTable, subCollection)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{enrichmentTable}{
Enrichment table, typically the \code{enrichmentTable} component of the output of
\code{\link{enrichmentAnalysis}}.
}
  \item{subCollection}{
A collection tht is a sub-collection of the collection used to calculate the enrichment table.
}
}
\details{
Data sets are matched based on names and IDs (if present in the enrichment table).
}
\value{
Numeric verctor of indices of rows of the enrichment table in which the data sets are present in the given
\code{subCollection}.}
\author{
Peter Langfelder
}
\keyword{misc}% __ONLY ONE__ keyword per line
