\name{enrichmentLabels}
\alias{enrichmentLabels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize an enrichment table
}
\description{
This function creates a data fram with a few most highly enriched terms for each class, and creates an
"enrichment label" for each class that summarizes the enrichement on a single line.
}
\usage{
enrichmentLabels(
    # Enrichment table and the relevant columns
    enrichmentTable,
    classCol = "class",
    pValueCol = "Bonferroni",
    setNameCol = "dataSetName",
    groupCol = "inGroups",
    fractionCol = "fracOfEffectiveClassSize",
    classSizeCol = "classSize",

    classLabels = NULL,

    # Restrictions on terms to be included
    pValueThreshold = 0.05,
    minSize = 0.1,

    # Which groups to include as separate entries
    focusOnGroups = c("all", "GO.BP", "GO.MF", "GO.CC", "Brain region markers", "Cell type markers"),
    excludeGroups = NULL,
    groupShortNames = NULL,
    numericClassLabels = TRUE,

    entrySeparator = "|",

    # Optional: collection or at least group information (component of collection) if implied groups are to be used
    collection = NULL,
    groupList = NULL,

    # What to include in the one-line summary
    includeClassSize = TRUE,
    outputClassPrefix = if (includeClassSize) "M." else "",
    orderBySignificance = TRUE,
    includeFraction = TRUE,
    fractionDigits = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{enrichmentTable}{
  Enrichment table: a data frame where each row corresponds to a class-term combination and gives the
overlap and enrichment statistics for the class-term overlap
}
  \item{classCol}{
 Name of the column in \code{enrichmentTable} that contains the class labels.
}
  \item{pValueCol}{
 Name of the column in \code{enrichmentTable} that contains the enrichment p-value.
}
  \item{setNameCol}{
 Name of the column in \code{enrichmentTable} that contains the gene set name.
}
  \item{groupCol}{
 Name of the column in \code{enrichmentTable} that contains the groups to which the gene set belongs.
}
  \item{fractionCol}{
 Name of the column in \code{enrichmentTable} that contains the overlap as a fraction of the class size.
}
  \item{classSizeCol}{
 Name of the column in \code{enrichmentTable} that contains the class size.
}
  \item{classLabels}{Optional vector of class labels that were used to carry out the enrichment analysis (function
\code{enrichmentAnalysis}). If given, class
sizes are calculated from this vector; otherwise, they are taken from the \code{classSizeCol}. Note that
\code{enrichmentAnalysis} restricts class labels to those elements that have a valid identifier and hence the class sizes
reported in \code{enrichmentAnalysis} may be smaller.}

  \item{pValueThreshold}{
 Threshold of the values contained in \code{pValueCol} for the term to be included in the labels.
}
  \item{minSize}{
 Minimum fraction of the class size that the overlap should be to get priority.
}
  \item{focusOnGroups}{
 A charcter vector giving the groups whose strongest enriched gene sets should be included in the
enrichment summary.
}
  \item{excludeGroups}{
 Optional character vector giving the groups whose data sets should be excluded when constructing the enrichment label.
This argument takes precedence over \code{focusOnGroups}, that is, if a data set belongs to one of the
\code{focusOnGroups} and \code{excludeGroups} groups, it will be excluded from the enrichment label.}

  \item{groupShortNames}{
 Optional specifications of shortened names for the groups to appear in the one-line summary.
}
  \item{numericClassLabels}{
 Logical: should class labels be assumed to be numeric? This only affects the final ordering of the output
table.
}
  \item{collection}{Optional collection that was used as reference collection in \code{enrichmentAnalysis}. When given,
group information in the collection is used to determine groups implied by \code{focusOnGroups}; these are effectively added to
\code{focuesOnGroups} beforeenrichment terms are selected. }

  \item{groupList}{Optional list of groups (see \code{newGroup}). When given, the parent-child relationships are used to determine groups implied by \code{focusOnGroups}; these are effectively added to
\code{focuesOnGroups} beforeenrichment terms are selected. }

  \item{entrySeparator}{Character (or character string) that will be used to separate multiple entries
within a single character string.}

  \item{includeClassSize}{
 Logical: should class sizes be included in the one-line summary?
}
  \item{outputClassPrefix}{
 A character string giving the prefix to be prepended in front of the class name for the one-line summaries.
}
  \item{orderBySignificance}{Logical: should the enrichment terms in each resulting label be ordered by significance? }

  \item{includeFraction}{
 Logical: should the overlap fraction be included in the one-line summaries?
}
  \item{fractionDigits}{
 Integer giving the number of digits the fraction should be rounded to.
}
}
\details{
This function is intended to provide a short summary of highest-enriched terms for each class (see below
how the highest-enriched set is determined). It outputs a
data frame with one row for each class (as opposed to one row per class/gene set combination). It
identifies the highest enriched gene set overall, and highest-enriched gene set in each of the
\code{focusOnGroups} groups, and puts the information together. Finally, it also puts together a single
character string that summarizes the information; this can be used as an enrichment-derived label for each
class (hence the function name). An example output enrichment label follows:

\code{M.1 (1380): intracellular part (0.74; All); RedBloodCell (0.09; CT)}

Here M. is the output class prefix \code{outputClassPrefix}, 1 is class label, 1380 is class size,
\code{intracellular part} is the highest enriched gene set (in this case a GO term), 0.74 is the fraction
of class 1 genes that are also part of the gene set \code{intracellular part}, and All indicates this is
the highest-enriched gene set among all sets. The \code{RedBloodCell} gene set is the highest enriched gene
set from among the cell type (CT) markers, and 9\% (=0.09) of the class 1 genes are also in the gene set
\code{RedBloodCell}. Gene sets from other groups are not listed because they did not pass the p-value
threshold of 0.05.

The highest-enriched gene set in each group is the set with the most significant p-value from among those
that contain at least \code{minSize} fraction of the class genes and pass the p-value threshold. 
If none of the gene sets pass the overlap
size and p-value criteria, the gene set with the most significant p-value is selected if it passes the
p-value threshold. If none of the gene sets pass the p-value threshold, the corresponding entries in the
output are empty.

}
\value{
A data frame containing the following columns:
\item{enrichmentLabel}{A character string summary of highest-enriched terms (see details).}
\item{class}{class label}
\item{classSize}{class size}
Further columns give information about highest enriched gene sets overall, and in each group specified in
\code{focusOnGroups}.

}
\author{
Peter Langfelder
}
\seealso{
  \code{\link{enrichmentAnalysis}} for enrichment calculations
}
\keyword{misc}% __ONLY ONE__ keyword per line
