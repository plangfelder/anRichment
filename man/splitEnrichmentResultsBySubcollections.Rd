\name{splitEnrichmentResultsBySubcollections}
\alias{splitEnrichmentResultsBySubcollections}
\alias{splitEnrichmentTableBySubcollections}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Split results of enrichment analysis into parts corresponding to subcollections
}
\description{
Split results of enrichment analysis into parts corresponding to subcollections of the reference collection
used for the enrichment analysis.
}
\usage{
splitEnrichmentResultsBySubcollections(enr, combinedCollection, collections, dropColumns = character(0),
        includeComponents = c("countsInDataSet", "pValues", "enrichmentRatio"))
splitEnrichmentTableBySubcollections(enrTab, combinedCollection, collections, dropColumns = character(0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{enr}{
A list with results of enrichment analysis. Typically the value of \code{\link{enrichmentAnalysis}}.
}
  \item{enrTab}{
Enrichment table, typically the \code{enrichmentTable} component of the output of
\code{\link{enrichmentAnalysis}}.
}
  \item{combinedCollection}{
The combined collection against which the enrichment calculation was performed. Usually the argument
\code{refCollection} to \code{\link{enrichmentAnalysis}}.
}
  \item{collections}{
A list in which each component is a collection. Each collection should be a subcollection of
\code{combinedCollection} in the sense that the gene sets contained in the subcollection are also present in
\code{combinedCollection}.
}
  \item{dropColumns}{
Optional character vector of column names of \code{enrTab} or the \code{enrichmentTable} component of
\code{enr} that should be dropped in output. 
}
  \item{includeComponents}{
  Character vector specifying component tables in argument \code{enr} that are to be split. The components must be
matrices with rows corresponding to gene sets and columns to classes. Within the object returned by 
function \code{\link{enrichmentAnalysis}}, the relevant components are \code{pValues}, \code{countsInDataSet},
\code{enrichmentRatio} and optionally \code{Bonferroni} and \code{FDR}.}
}
\details{
Although the \code{\link{enrichmentAnalysis}} returns a list with multiple components, the function
\code{splitEnrichmentResultsBySubcollections} only splits and returns the components
\code{enrichmentTable}, \code{countsInDataSets} and \code{pValues}.
}
\value{
Both functions return a list with each component corresponding to onbe of the collections in
\code{collections}. 

In the value of \code{splitEnrichmentResultsBySubcollections}, each component is a list with three components,
\code{enrichmentTable}, \code{countsInDataSets} and \code{pValues} in which the rows (corresponding to
reference data sets) are restricted to those in the relevant subcollection.

In the value of \code{splitEnrichmentTableBySubcollections}, each component is the input \code{enrTab} in
which the rows (corresponding to 
reference data sets) are restricted to those in the relevant subcollection.

}
\author{
Peter Langfelder
}
\seealso{
\code{\link{enrichmentAnalysis}}
}
\keyword{misc}
