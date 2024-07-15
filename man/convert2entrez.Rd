\name{convert2entrez}
\alias{convert2entrez}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Convert gene symbols or RefSeq indentifiers to Entrez.
}
\description{
  Convert gene symbols or RefSeq indentifiers to Entrez using Bioconductor's annotation packages.
}
\usage{
convert2entrez(
         organism, 
         symbol = NULL, refSeq = NULL, 
         ignoreCase = FALSE,
         dropNonConverted = FALSE,
         useAlias = TRUE,
         noMatch = NA, 
         multipleMatches = NA,
         annotation = NULL,
         symbolCol = "symbol",
         refSeqCol = "refSeq",
         aliasCol = "synonyms",
         aliasSep = "; ",
         aliasSplit.fixed = TRUE,
         entrezCol = "Entrez", 
         ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{organism}{
     Character string identifying the organism. Can be a common name (e.g., "human"), scientific name
(e.g., "Homo sapiens"), or the scientific name shorthand (e.g., "Hs"). To see available organisms, use
\code{\link{organismLabels}} with the default arguments (i.e., 'organismLabels()').
}
  \item{symbol}{
     Character vector of gene symbols or a list of such character vectors.
}
  \item{refSeq}{
     Character vector of RefSeq identifiers, or a list of such vectors.
}
  \item{ignoreCase}{logical: should upper- and lower-case characters be treated as the same?}

  \item{dropNonConverted}{logical: should entries whose Entrez cannot be found be dropped from the output?
If such entries are retained, their Entrez ID is set to \code{NA}.}

  \item{useAlias}{logical: should gene alias symbols be searched for a unique match as well?}

  \item{noMatch}{Value to return for unmatched input elements.}

  \item{multipleMatches}{Value to return for input elements that have multiple matches.}

  \item{annotation}{Optional annotation data frame that, when given, replaces the annotation database specified by
\code{organism}. }

  \item{symbolCol}{Column name in \code{annotation} that contains symbols.}

  \item{refSeqCol}{Column name in \code{annotation} that contains refSeq identifiers.}

  \item{aliasCol}{Column name in \code{annotation} that contains aliases (alternative gene symbols) separated by
\code{aliasSep}. }

  \item{aliasSep}{Character string that separates aliases.}

  \item{aliasSplit.fixed}{Logical; this argument is passed on to \code{\link{strsplit}} (as \code{fixed}) when splitting
values in \code{aliasCol} into individual aliases.}

  \item{entrezCol}{Column name in \code{annotation} that contains Entrez identifiers.}

  \item{...}{Other arguments, currently undocumented.}
}
\details{
If \code{annotation} is given, it is used as a lookup/translation table between symbols or refSeq identifiers and Entrez
identifiers. Otherwise, the appropriate Bioconductor organism database is used.

The function first searches all official gene symbols for
matches of the given \code{symbol}. If the input \code{useAlias} is \code{TRUE}, 
\code{symbol} entries with no match among official symbols will be matched against alias symbols and those for
which a unique match is found are returned.

}

\value{
  A vector of Entrez identifiers if the input identifiers formed a single atomic vector, or a list of such
vectors if the input identifiers were a list.
}
\author{
  Peter Langfelder
}
\note{The mapping from gene symbol to Entrez is often not unique, and the function may print a warning to
that effect.}

\seealso{
 \code{\link{convert2symbol}} for the Entrez to gene symbol conversion.
}
\examples{
  require(org.Hs.eg.db);
  symbols = c("APOE", "LRRN3", "Wrong");
  convert2entrez(organism = "human", symbol = symbols)
}
\keyword{misc}% __ONLY ONE__ keyword per line
