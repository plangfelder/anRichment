\name{buildGOcollection}
\alias{buildGOcollection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Build a collection object that contains GO terms as its gene lists.
}
\description{
  This function creates a collection object with GO terms for the given organism as its gene sets. It
defines 4 groups: GO (for all terms), GO.BP (for GO Biological Process terms), GO.MF (for GO Molecular
Function terms), and GO.CC (for GO Cellular Component terms).
}
\usage{
buildGOcollection(organism, 
                 termNames = NULL, 
                 includeOffspring = TRUE, 
                 strict = TRUE, 
                 verbose = 2, indent = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{organism}{
     Character string specifying the organism for which to build the collection. Can be a common name (e.g.,
"human"), scientific name
(e.g., "Homo sapiens"), or the scientific name shorthand (e.g., "Hs"). To see available organisms, use
\code{\link{organismLabels}} with the default arguments (i.e., 'organismLabels()').
}
  \item{termNames}{
     Character vector giving the term names (or GO identifiers) of terms to be included. If \code{NULL},
all terms will be included.
}
  \item{includeOffspring}{
     Logical: should genes in offspring terms of a term X be included in the gene list for term X?
}
  \item{strict}{
     Logical: should the function throw and error (\code{TRUE}) or a warning (\code{FALSE}) if a supplied
\code{termName} cannot be found among the existing GO terms?
}
  \item{verbose}{
     Controls the level of verbosity. Higher numbers mean more diagnostic output.
}
  \item{indent}{
     Controls the indentation of the diagnostic messages. Each unit adds two spaces.
}
}
\details{ This function requires the GO.db package, as well as the appropriate org.Xx.eg.db package for 
the input organism.
 }

\value{
  A collection object with GO terms as its gene lists.
}
\references{
 The Gene Ontology Consortium. Gene ontology: tool for the unification of biology. Nat. Genet.. May
2000;25(1):25-9


}
\author{
  Peter Langfelder
}

\seealso{

\code{\link{newCollection}} for creating new empty collections,

\code{\link{internalCollection}} for a collection of gene sets from the literature.
 
}
\examples{
GOcollection = buildGOcollection(organism = "human") # Can also use synonyms "Homo sapiens" and "Hs"
# GOcollection contains these groups:
knownGroups(GOcollection)

}
\keyword{misc}% __ONLY ONE__ keyword per line
