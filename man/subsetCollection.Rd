\name{subsetCollection}
\alias{subsetCollection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Subset a collection by groups and/or gene sets
}
\description{
  This function subsets a collection to given groups and/or gene sets.
}
\usage{
subsetCollection(
   collection, 
   tags = NULL, 
   matchComponents = c("ID", "name", "groups", "alternateNames", "source",
                       "groupAlternateNames", "nameAndAlternates", "groupsAndAlternates"),
   searchType = c("any", "all"), invertSearch = FALSE,
   exactMatch = TRUE, fixed = TRUE, ignore.case = TRUE,
   firstDate = NULL, lastDate = NULL,
   dateFormat = "\%Y-\%m-\%d", invertDateSearch = FALSE)
}
\arguments{
  \item{collection}{
     A collection to be subset.
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

  \item{invertSearch}{Logical: should the search be inverted, that is the result will contain all data sets
except ones that match the sarch pattern?}

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

}
\value{
  The subset collection.
}
\author{
  Peter Langfelder
}
\seealso{
  \code{\link{newCollection}} for creating new collections.
}
\keyword{misc }% __ONLY ONE__ keyword per line
