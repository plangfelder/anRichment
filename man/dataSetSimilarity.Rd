\name{dataSetSimilarity}
\alias{dataSetSimilarity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute a similarity matrix between data sets
}
\description{
This function implements a calculation of data similarity matrix for all pairs of data sets in a collection.
The input can be a collection, or a pre-calculated matrix representation of the collection.
}
\usage{
dataSetSimilarity(
    collection,

    # background against which to calculate obverlaps. defaults to all genes in the collection.
    background = NULL,

    # Subsetting options
    tags = NULL,
    matchComponents = c("ID", "name", "groups", "alternateNames", "source",
                        "groupAlternateNames", "nameAndAlternates", "groupsAndAlternates"),
    searchType = c("any", "all"), invertSearch = FALSE,
    exactMatch = TRUE, fixed = TRUE, ignore.case = TRUE,

    # Special subsetting options
    includeProperties = TRUE,
    includeGeneSets = TRUE,
    dropZeroOverlapSets = FALSE,

    namesFrom = c("ID", "name"),

    # Alternative to all of the above: pre-calculated data set matrix
    dataSetMat = NULL,

    # Distance options
    type = c("product", "cosine", "cov", "cor"),
    productNormalization = c("min", "mean", "max"),
    TOMType = c("none", "min", "mean"),
    verbose = 1,
    indent = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{collection}{
A collection of "data sets" (gene sets or gene properties)
}
  \item{background}{
Optional specification of background. If not given, all genes appearing in the collection (after optional
subsetting below) will be used.
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

    \item{includeProperties}{
Logical: should gene properties be included in the matrix?
}
  \item{includeGeneSets}{
Logical: should gene sets be included in the matrix?
}
  \item{dropZeroOverlapSets}{
Logical: should gene sets with zero overlap with the specified background be dropped from the matrix?
}
  \item{namesFrom}{
Determines the component of the gene set or property meta data that will be used as column names of the
result.
}
\item{dataSetMat}{
Optional matrix representation of the collection. If given, all of the above arguments are disregarded; if
not given, the above arguments are used to calculate it.
}
  \item{type}{
Similarity type.
}
  \item{productNormalization}{
Normalization for the product similarity. Only used when \code{type="product"}.
}
  \item{TOMType}{
One of \code{"none"}, \code{"min"}, \code{"mean"}, indicating whether the similarity should include a TOM
transformation, and which denominator should be used in the calculation of TOM.
}
  \item{verbose}{
     Controls the level of verbosity. Higher numbers mean more diagnostic output.
}
  \item{indent}{
     Controls the indentation of the diagnostic messages. Each unit adds two spaces.
}
}
\details{
If the matrix representation of the collection is not supplied in argument \code{dataSetMat}, the function
first calls \code{\link{dataSetMatrix}} to create the matrix representation, possibly after first subsetting
the collection according to the supplied \code{tags} and other subsetting arguments.

The set-setsimilarity is then calculated from the matrix representation in which each columns corresponds to
a data set. The product similarity equals the scalar product of the two vectors, divided (normalized) by
1 (no normalization), minimum or mean of the sums of each column as indicated in argument
\code{productNormalization}. Other available similarities include the cosine (product divided by the product
of norms of the vector), covariance and correlation of the columns.

Finally, the similarity matrix is optionally transformed using the topological overlap transformation.

}
\value{
A similarity matrix.
}
\author{
Peter Langfelder
}
\seealso{
\code{dataSetMatrix} to create the matrix representation of a collection.
}
\keyword{misc}% __ONLY ONE__ keyword per line
