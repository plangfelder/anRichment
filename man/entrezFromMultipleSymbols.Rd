\name{entrezFromMultipleSymbols}
\alias{entrezFromMultipleSymbols}
\title{
Convert gene symbols to Entrez IDs, with possibly multiple symbols in each entry
}
\description{
Convert gene symbols to Entrez IDs, with possibly multiple symbols in each entry. 

}
\usage{
entrezFromMultipleSymbols(symbol, split = ";", organism = "human", ignoreCase = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{symbol}{
A character vector of gene symbols. Each entry may contain multiple symbols separated by \code{split}. 
}
  \item{split}{
A character or string giving the separator of when there are multiple gene symbols in one element of
\code{symbol}. 
}
  \item{organism}{
Organism for which the Entrez IDs are desired.
}
  \item{ignoreCase}{
Logical: should case (upper- vs. lowe-case) be ignored when matching symbols?
}
}
\details{
The given character strings are first split into symbols using the argument \code{split} as the divider. Each
fraction is then interpreted as a single symbol. Symbols are matched against official gene symbols; ones that
cannot be matched are then matched against aliases. 
When multiple symbols are
present, an Entrez ID is assigned only if all of the gene symbols in that element map to the same Entrez ID
(or to no Entrez ID at all)

}
\value{
A vector of Entrez IDs. Elements that cannot be mapped to a unique Entrez have NA in the output.
}
\author{
Peter Langfelder
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
The workhorse function \code{\link{convert2entrez}}.
}
\keyword{misc}% __ONLY ONE__ keyword per line
