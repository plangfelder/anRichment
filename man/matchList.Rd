\name{matchList}
\alias{matchList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A faster way to match a list of vectors to common reference
}
\description{
This function implements a faster alternative to \code{lapply(lst, match, ref, incomparables = incomparables)}.
}
\usage{
matchList(lst, ref, nomatch = NA_integer_, incomparables = NA, setAttributes = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lst}{
  A list of vectors to be \code{match}ed against \code{ref}.
}
  \item{ref}{
  A vector to be matched against.
}
  \item{setAttributes}{
  Logical: should attributes of \code{lst} (such as names, potential dimensions and dimnames) be copied into the output?
}
  \item{nomatch}{Value that signifies a component without a match in \code{ref}.}

  \item{incomparables}{Vector of values that cannot be matched.}
}
\value{
  A list of match indices.
}
\author{
Peter Langfelder
}
\seealso{
  The workhorse function \code{\link{match}}.
}
% \keyword{misc}
% \concept{match}
