\name{knownEvidenceCodes}
\alias{knownEvidenceCodes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Known evidence codes.
}
\description{
  Returns a data frame containing known evidence codes and their description.
}
\usage{
knownEvidenceCodes()
}
\details{
  Most of the known evidence codes come from the GO consortium.
}
\value{
  A data frame containing the following columns
  \item{evidenceCode}{Short evidence codes.}
  \item{evidenceDescription}{Description of the evidence code.}
  \item{evidenceType}{Type of the evidence code.}
}
\references{
 GO consortium web site, especially http://geneontology.org/page/guide-go-evidence-codes
}
\author{
  Peter Langfelder
}
\seealso{
 \code{\link{matchEvidenceCode}}
}
\examples{
  knownEvidenceCodes()
}
\keyword{ misc}% __ONLY ONE__ keyword per line
