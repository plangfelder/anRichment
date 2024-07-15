\name{collection2dataFrames}
\alias{collection2dataFrames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Turn a collection into data frames
}
\description{
Turn a collection into data frames that can be saved as spreadsheets.
}
\usage{
collection2dataFrames(
    collection, 
    sep = "|", 
    propertyNamesFrom = "ID",
    organismName = c("common", "scientific", "short"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{collection}{
A collection.
}
  \item{sep}{
For components that have several values (such as groups), this argument will be used as separator.
}
  \item{propertyNamesFrom}{
Component whose value will be used for column names of numerical properties. Recommended values are
\code{"ID"} or \code{"name"}.
}
  \item{organismName}{A character string specifying which of the three versions of each organism name to be
included in the generated data frames. Must be one of "common", "scientific", or "short".}

}
\details{
A collection is a structure that contains gene sets, numerical or categorical gene properties, and some
associated information. This function turns the entire information stored in a collection into data frames
that can be easily saved as plain text tables or spreadsheets.
}
\value{
  A list with the following components:
  \item{geneSetInfo}{Attributes of gene sets, one gene set per row.
  Columns include ID, name, short name, description, source, organism, internal classification, and groups.
}
  \item{dataPropertyInfo}{Attributes of gene properties, one gene set per row.
  Columns include ID, name, short name, description, source, organism, internal classification, groups, and
weight index. }

  \item{geneSetContent}{A data frame containing the gene set content. Each row corresponds to one gene in a
gene set. Columns give gene set ID, gene set
name, gene Entrez ID, evidence code, and source. For storage efficiency, columns ID, name and source are
returned as factors since it is expected that they will contain multiple repeats of the same character
strings.}

  \item{genePropertyContent}{A data frame giving the values of gene properties. Each row corresponds to a
gene and each column to a gene property (i.e., a row in \code{dataSetInfo.geneProperties} above). }

  \item{genePropertyWeights}{A data frame giving the weights of each gene in gene properties. Each row
corresponds to a gene. The weights of genes in each gene property are encoded in the column given by the
\code{weightIndex} entry of the corresponding row in \code{dataSetInfo.geneProperties}. The first column of
this data frame is typically unit (i.e., equal) weights; other columns may or may not be present.}

  \item{evidenceCodes}{A data frame that summarizes known evidence codes. This can be used to interpret the
numeric evidence codes used in \code{geneSetContent}.}
}
\author{
Peter Langfelder
}

\seealso{
\code{\link{newGeneSet}}, \code{\link{newGeneProperty}}, \code{\link{newGroup}},
\code{\link{newCollection}} for creating new gene sets, gene properties, groups and collections.
}
\keyword{misc}% __ONLY ONE__ keyword per line
