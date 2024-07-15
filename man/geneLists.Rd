\name{geneLists}
\alias{geneLists}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extract Entrez identifiers from the given collection.
}
\description{
  Extracts Entrez identifiers of genes with given evidence from gene lists in the given \code{collection}. 
Only one copy of each Entrez is retained for each gene list.
}
\usage{
geneLists(
  collection, 
  evidence = "all",
  tags = NULL, 
  matchComponents = c("ID", "name", "groups", "alternateNames", "source",
                      "groupAlternateNames", "nameAndAlternates", "groupsAndAlternates"),
  searchType = c("any", "all"), invertSearch = FALSE,
  exactMatch = TRUE, fixed = TRUE, ignore.case = TRUE,
  firstDate = NULL, lastDate = NULL,
  dateFormat = "\%Y-\%m-\%d", invertDateSearch = FALSE,
  namesFrom = c("ID", "name"),
  simplify = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{collection}{
     A collection object.
}
  \item{evidence}{
     Specifies evidence codes or genes that are to be included in the returned lists. The special value
\code{"all"} will return all genes. Available evidence codes can be obtained using
\code{knownEvidenceCodes()}.
}
  \item{tags}{Optional specification of search terms (group names, gene set names, or gene set IDs). If
given, gene lists will only be returned for gene sets that match the tags according to further search
specifications.}

  \item{matchComponents}{Character vector denoting which of the identifying information should be searched for
the \code{tags} given above. The default is to search all of the searchable fields. Values
"nameAndAlternates" and "groupAndAlternates" are shorthands for "name", "alternateNames", and "groups", 
"alternateGroupNames", respectively.}

  \item{searchType}{One of \code{"any"} or \code{"all"}. For \code{"any"}, data sets that match at
least
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


  \item{namesFrom}{One of \code{"ID", "name"}. Specifies where the \code{names} attribute of the returned
list will be taken from: the identifier (\code{ID}) or the name (\code{name}) of each gene set.}

  \item{simplify}{Logical: if only one set matches, should the resulting one-component list be simplified
by simply returning the first component?}

}
\value{
  A list in which every component contains a vector containing the Entrez identifiers of the genes from the
corresponding gene set in \code{collection}. Each component is named by either the identifier or the name
of the corresponding gene set.
}
\author{
  Peter Langfelder
}

\seealso{
  \code{\link{matchEvidenceCode}} for details on specifying evidence codes;
  \code{\link{knownEvidenceCodes}} for all known evidence codes;
  \code{\link{grep}} for partial and regular expression matching.
}
\keyword{misc}% __ONLY ONE__ keyword per line
