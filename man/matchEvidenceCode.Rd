\name{matchEvidenceCode}
\alias{matchEvidenceCode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Match given evidence to all known evidence codes.
}
\description{
  Returns the numeric code of the given evidence character string. 

}
\usage{
matchEvidenceCode(evidence)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{evidence}{
     Character vector that specifies the evidence(s) to be matched. See Details.
}
}
\details{
The translation between characters and
numeric codes is encoded in the output of \code{knownEvidenceCodes()}. This function will return the numeric
codes of all rows in which the input \code{evidence} is found. Note that while individual evidence codes
and descriptions are unique, a single evidence type typically  encompasses several evidence codes, and is
also a valid input for this function.

Duplicates are removed from the returned codes and the codes are sorted, i.e., they do not in general
correspond to the order of the entries in \code{evidence}. 

If the input \code{evidence} is \code{"all"}, all knwon numeric codes will be returned. The query
\code{"all"} must be the only element of \code{evidence}.
}
\value{
  A vector of numeric evidence codes that match the given \code{evidence}.
}
\author{
  Peter Langfelder
}
\seealso{
  \code{knownEvidenceCodes} for all recognized evidence codes and types.
}
\examples{

# Display all evidence codes:
knownEvidenceCodes()

# Numeric code for evidence codes "IEP" and "IKR":
matchEvidenceCode(c("IEP", "IKR"))

# Numeric codes for all "Experimental" evidence codes:
matchEvidenceCode("Experimental")

# All numeric codes:

matchEvidenceCode("all")

# 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
