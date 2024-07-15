\name{representativeProbes}
\alias{representativeProbes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select a representative probe
}
\description{
Select a representative probe (variable) based on given data.
}
\usage{
representativeProbes(
   mtmat,

   # Mapping of variables to genes: either Entrez...
   entrez = NULL,

   #... or an annotation table with other necessary information
   annot = NULL,
   annotIDCol = "Name",
   entrezCol = NULL,
   symbolCol = "UCSCRefGeneName",
   symbolSep = ";",
   organism = "human",

   # Optional filtering based on annotation
   keepValues = NULL,
   keepValueColumns = NULL,
   keepPatterns = NULL,
   keepPatternColumns = NULL,
   ignoreCase = TRUE,

   # Options governing selection of representatives
   selectionMethod = "maxVariance",
   selectionStatisticFnc = NULL,
   useGroupHubs = TRUE,
   connectivityPower = 6,
   minProportionPresent = 1,
   statisticFncArguments = list(),
   adjacencyArguments = list(),

   # Consensus options
   consensusQuantile = 0,
   calibration = "full quantile",

   getRepresentativeData = TRUE,
   dataToReduce = NULL,
   verbose = 2, indent = 0)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mtmat}{
A matrix, or a \code{\link{multiData}} structure with marices, on which to base the representative probe
selection.
}

  \item{entrez}{
A vector of Entrez identifiers for each variable (column) in \code{mtmat}.
}

  \item{annot}{
An optional data frame or matrix of annotations for the variables (columns) of \code{mtmat}. The row order
does not have to correspond to the rows in \code{mtmat}; see \code{annotIDCol} below for data to annotation
matching.
}

  \item{annotIDCol}{
Name or number of the column in \code{annot} that contains variable names, to be matched against the column names
in \code{mtmat}.
}

  \item{entrezCol}{
Name or number of the column in \code{annot} that contains Entrez identifiers for each variable/probe. If not
given, symbols (specified by \code{symbolCol} below) will be converted to Entrez IDs.
}

  \item{symbolCol}{
Name or number of the column in \code{annot} that contains gene symbols. This is used only if \code{entrezCol}
is missing.
}

  \item{symbolSep}{
Separator that separates the (potentially multiple) gene symbols for each variable/probe.
}

  \item{organism}{
Organism for which the analysis is carried out.
}

  \item{keepValues}{
Optional filter of probes in \code{annot}. Can be a vector or a list of values that must be present in 
\code{annot} columns specified in \code{keepValueColumns}. Values must be matched exactly.
}

  \item{keepValueColumns}{
A vector of column names that correspond to the filter values specified in \code{keepValues}.
}

  \item{keepPatterns}{
A second optional filter of probes in \code{annot}. Can be a vector or a list of regular expressions 
that match values in \code{annot} columns specified in \code{keepPatternColumns}. Regular expressions are
matched against candidate values using \code{\link{grep}}.
}

  \item{keepPatternColumns}{
A vector of column names that correspond to the \code{grep} filter values specified in \code{keepValues}.
}

  \item{ignoreCase}{
Logical: should case be ignored when filtering variables?
}

  \item{selectionMethod}{character string for determining which method is used to choose the representative
(when \code{useGroupHubs} is \code{TRUE}, this method is only used for groups with 2
variables).
The following values can be used:
"MaxMean" (default) or "MinMean" return the variable with the highest or lowest mean value, respectively;
"maxRowVariance" return the variable with the highest variance;
"absMaxMean" or "absMinMean" return the variable with the highest or lowest mean absolute value; and
"function" will call a user-input function (see the description of the argument
\code{selectionStatisticFnc}). The built-in functions can be instructed to use robust analogs (median and
median absolute deviation) by also specifying \code{statisticFncArguments=list(robust = TRUE)}.
}

  \item{selectionStatisticFnc}{User-supplied function used to calculate the selection statistic when
\code{method} above equals \code{"function"}.  The function must take argumens \code{x} (a matrix) and
possibly other arguments that can be specified using \code{statisticFncArguments} below. The return value
must be a vector with one component per column of \code{x} giving the selection statistic for each column.
}

  \item{useGroupHubs}{Logical: if \code{TRUE}, groups with 3 or more variables will be
represented by the variable with the highest
connectivity according to a signed weighted correlation network adjacency matrix among the corresponding
rows. The connectivity is defined as the row sum of the adjacency matrix. The signed weighted
adjacency matrix is defined as A=(0.5+0.5*COR)^power where power is determined by the argument
\code{connectivityPower} and COR denotes the matrix of pairwise correlation coefficients among the
corresponding rows. Additional arguments to the underlying function \code{\link{adjacency}} can be specified
using the argument \code{adjacencyArguments} below.
}

  \item{connectivityPower}{
Positive number (typically integer) for specifying the soft-thresholding power used
to construct the signed weighted adjacency matrix, see the description of \code{useGroupHubs}.
This option is only used if \code{useGroupHubs} is \code{TRUE}.
}

  \item{minProportionPresent}{
A number between 0 and 1 specifying a filter of candidate probes. Specifically, for each group, the variable
with the maximum consensus proportion of present data is found. Only variables whose consensus proportion of
present data is at least \code{minProportionPresent} times the  maximum consensus proportion are retained as
candidates for being a representative.
}

  \item{consensusQuantile}{
A number between 0 and 1 giving the quantile probability for consensus calculation.
0 means the minimum value (true consensus) will be used.
}

  \item{calibration}{
Character string describing the method of calibration of the selection statistic among
the data sets. Recognized values are \code{"none"} (no calibration) and \code{"full quantile"} (quantile
normalization). 
}

  \item{statisticFncArguments}{
 A list giving further arguments to the selection statistic function. Can be
used to supply additional arguments to the user-specified \code{selectionStatisticFnc}; the value
\code{list(robust = TRUE)} can be used with the built-in functions to use their robust variants.
}

  \item{adjacencyArguments}{
Further arguments to the function \code{adjacency}, e.g.
\code{adjacencyArguments=list(corFnc = "bicor", corOptions = "use = 'p', maxPOutliers = 0.05")} will select
the robust correlation \code{\link{bicor}} with a good set of options. Note that the \code{\link{adjacency}}
arguments \code{type} and \code{power} cannot be changed.
}

  \item{getRepresentativeData}{
Logical: should the representative data, i.e., \code{mtmat} restricted to
the representative variables, be returned?
}

  \item{dataToReduce}{
Optional additional data (e.g., module labels) that correspond to columns in \code{mtmat} and should be
reduced to representatives. This can be a vector or a matrix in which each _row_ corresponds to a variable in
\code{mtmat}.
}

  \item{verbose}{
Level of verbosity; 0 means silent, larger values will cause progress messages to be printed.
}

  \item{indent}{
Indent for the diagnostic messages; each unit equals two spaces.
}
}
\details{
This function is mostly a wrapper for \code{\link[WGCNA]{consensusRepresentatives}}, with extra code intended
to make the use of \code{\link[WGCNA]{consensusRepresentatives}} easier for microarray and especially
methylation data. 
}
\value{
\item{representatives}{A named vector giving, for each group, the selected representative (input \code{rowID}
or the variable (column) name in \code{mdx}). Names correspond to groups.}
\item{varSelected}{A logical vector with one entry per variable (column) in input \code{mdx} (possibly
after restriction to variables occurring in \code{colID}), \code{TRUE} if the column was selected as a
representative.}
\item{representativeData}{Only present if \code{getRepresentativeData} is \code{TRUE};
the input \code{mdx} restricted to the representative variables, with column
names changed to the corresponding groups.}
\item{reducedData}{The input \code{dataToReduce} reduced to representatives; only present if
\code{dataToReduce} was non-NULL.}
}
\author{
Peter Langfelder
}
\seealso{
\code{\link[WGCNA]{consensusRepresentatives}} for the actual representative selection function.
}
\keyword{misc}
