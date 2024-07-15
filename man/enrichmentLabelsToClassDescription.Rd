\name{enrichmentLabelsToClassDescription}
\alias{enrichmentLabelsToClassDescription}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert module enrichment labels to suitable class descriptions
}
\description{
Converts enrichment labels (such as those returned by \code{\link{enrichmentLabels}}) to class descriptions
that are needed by \code{\link{collectionFromClasses}}
}
\usage{
enrichmentLabelsToClassDescription(
  eLabels, 
  shortAnalysisName, 
  fullAnalysisName,
  shortNameFirst = TRUE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eLabels}{
A data frame with information about enrichment analysis of modules. Typically the output of
\code{\link{enrichmentLabels}}.
}
  \item{shortAnalysisName}{
A short version of the analysis name, to be added to the \code{enrichmentLabel} column in
\code{eLabels} to form the class name.
}
  \item{fullAnalysisName}{
A long analysis name, added to the \code{enrichmentLabel} column in
\code{eLabels} to form class description.
}
  \item{shortNameFirst}{Logical: Should the short analysis name precede the enrichment label in the final
class description?}
}
\value{
A data frame with the columns of the input, except \code{enrichmentLabel} is prepended by the
\code{shortAnalysisName}, and a new column \code{description} is added where the enrichment label is
explained in a bit more detail.
}
\author{
Peter Langfelder
}
\section{Warning}{
  The actual output enrichment labels and descriptions should be considered experimental and subject to
future changes.
}
\seealso{
\code{\link{enrichmentLabels}} to get the input \code{enrichmentLabels}, and
\code{\link{collectionFromClasses}} that typically uses the output of this function.
}
\keyword{misc}% __ONLY ONE__ keyword per line
