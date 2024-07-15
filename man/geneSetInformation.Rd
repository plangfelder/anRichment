\name{geneSetInformation}
\alias{geneSetInformation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return gene set meta-information 
}
\description{
This function returns a data frame with meta-information about selected gene sets in a collection.
}
\usage{
geneSetInformation(
   collection, 
   tags = NULL,
   matchComponents = c("ID", "name", "groups", "alternateNames", "source",
                       "groupAlternateNames", "nameAndAlternates", "groupsAndAlternates"),
   searchType = c("any", "all"), invertSearch = FALSE, 
   exactMatch = TRUE, fixed = TRUE, ignore.case = TRUE,
   firstDate = NULL, lastDate = NULL,
   dateFormat = "\%Y-\%m-\%d", invertDateSearch = FALSE,
   sep = "|", organismName = c("common", "scientific", "short")) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{collection}{
A collection.
}
  \item{tags}{
     A character vector that gives search terms that the selected data sets must match. The tags can be
data set IDs, names, and group names.
}
  \item{matchComponents}{Character vector denoting which of the identifying information should be searched for
the \code{tags} given above. The default is to search all of the searchable fields. Values
"nameAndAlternates" and "groupAndAlternates" are shorthands for "name", "alternateNames", and "groups", 
"alternateGroupNames", respectively.}

\item{searchType}{One of \code{"any"} or \code{"all"}. For \code{"any"}, data sets that match at least
one of the supplied tags will be returned. For \code{"all"}, data sets that match all of the supplied tags
will be returned.
}
  \item{invertSearch}{logical: should data sets that do NOT match the \code{dataSetTags} be added to the
group?
}
  \item{exactMatch}{Logical: should only gene sets matching the tag exactly be returned? If \code{FALSE},
partial matches (anywhere within the target string) will be allowed; the matching is then done using
\code{\link{grep}}.
}
  \item{fixed}{Logical: if \code{exactMatch} above is \code{FALSE}, should the supplied \code{tags} be
interpreted as literal search terms (\code{fixed=TRUE}) or as regular expressions (\code{FALSE})? This
argument is passed on to \code{\link{grep}}.
}

  \item{ignore.case}{Logical: should case (uppercase vs. lowercase) be ignored when matching patterns? }

  \item{firstDate}{Optional specification of the earliest date for the data sets to be selected. Should be
either a \code{\link{Date}} object, or a character string specifying date. If a character string, its 
format can be specified using the argument \code{dateFormat} below. If \code{NULL}, dates will not be 
restricted from below.}

  \item{lastDate}{Optional specification of the latest date for the data sets to be selected. Should be
either a \code{\link{Date}} object, or a character string specifying date. If a character string, its
format can be specified using the argument \code{dateFormat} below. If \code{NULL}, dates will not be
restricted from above.}

  \item{dateFormat}{Specification of format of character dates. See \code{\link{strptime}} for the
description of the format string; the default means <year with century>-<month>-<day>. }


  \item{invertDateSearch}{Logical: should the date restriction be inverted, that is, data sets between
\code{firstDate} and \code{lastDate} be excluded from the operation?}

  \item{sep}{
For components that have several values (such as groups), this argument will be used as separator.
}
  \item{organismName}{A character string specifying which of the three versions of each organism name to be
included in the generated data frames. Must be one of "common", "scientific", or "short".}

}
\value{
A data frame with one row per matching gene set, giving the gene set ID, name, short name, description,
source, internal classification, groups, and number of genes. If no gene sets match the query, the function
returns \code{NULL}.
}
\author{
Peter Langfelder
}

\seealso{
\code{\link{collection2dataFrames}} for turning an entire collection to a series of data frames.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}% __ONLY ONE__ keyword per line
