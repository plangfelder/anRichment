\name{enrichmentAnalysis.general}
\alias{enrichmentAnalysis.general}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Enrichment analysis
}
\description{
  This function calculates enrichment statistics of given classes of genes in a reference collection.
}
\usage{
enrichmentAnalysis.general(
   classLabels = NULL,
   identifiers = NULL,
   classSuffixes = NULL,
   # Alternative input specification
   active = NULL,
   inactive = NULL,
   activeNames = names(active),
   # Optional map from active to inactive sets, may be useful to save memory
   active2inactive = NULL,

   # Reference collection: should be appropriate to the organism
   refCollection,

   # If data IDs and collections IDs aren't the same, a map from class/active to collection. One class ID may
   # map to at most one valid collection ID. 
   classID2setID = NULL,

   # optional sub-selection
   useGroups = NULL,

   # Specify what to use as the background list
   useBackground = c("intersection", "given", "reference", "custom"),
   customBackground = NULL,
   removeMissing = TRUE,

   # Included evidence levels
   useEvidence = "all",

   # Various calculation options
   removeDuplicatesInDifferentClasses = TRUE,
   ignoreLabels = NULL,
   alternative = "greater",

   # How much output to generate
   nBestDataSets = 10,
   threshold = 0.05,
   thresholdType = c("Bonferroni", "FDR", "nominal"),

   # Calculate multiple testing corrections? 
   getBonferroniCorrection = TRUE,
   getFDR = TRUE,

   # Return overlap features in the enrichment table?
   getOverlapIDs = TRUE,
   getOverlapSymbols = FALSE,
   ID2symbol = NULL,
   maxReportedOverlapGenes = 50,

   # Output options
   entrySeparator = "|",
   groupSeparator = entrySeparator,
   geneSeparator = entrySeparator,
   classColName = "class",

   # Return gene set details?
   getDataSetDetails = TRUE,

   # Diagnostic message options
   verbose = 1, indent = 0,

   # Undocumented arguments whose effect may change in the future
   ... )

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{classLabels}{
     A vector or a matrix of labels (such as class or cluster labels) 
     corresponding to the Entrez identifiers below. If a matrix, each column is taken to be a separate set
     of class labels. If given, overrides \code{active} and \code{inactive} below.
}
  \item{identifiers}{
     A vector of identifiers, one for each entry of \code{labels} if \code{labels} is a vector or
one for each row if \code{labels} is a matrix.
}
  \item{classSuffixes}{
     Optional character vector of length equal to the number of columns in \code{classLabels}. When putting together
enrichment results, \code{classSuffixes} will be added to levels of \code{classLabels} in each column to make them unique.
If not given, defaults to  " (label set 1)", " (label set 2)" etc.}

  \item{active}{
     Alternative specification of the input genes: a vector of identifiers considered 'active', or a
list of such vectors. See details.
}
  \item{inactive}{
     Alternative specification of the input: a vector of Entrez identifiers considered 'inactive'. For
convenience, genes contained in the \code{active} set above will automatically be removed, so this can also be
the full set of genes in the analysis. This can be a single atomic vector or a list. If list, it must either have 
the same number of components as \code{active} or the user must also supply a map from \code{active} to \code{inactive} in
\code{active2inactive}.
}
  \item{activeNames}{Names for the active sets to be used as class labels in the output tables. If \code{NULL},
the classes will all be labeled \code{"active"}.}

  \item{active2inactive}{Specifies a map from \code{active} to \code{inactive}; the inactive vector corresponding to entry
\code{i} of \code{active} is \code{active2inactive[i]}.
}
  \item{refCollection}{
     The reference collection containing gene sets to be evaluated for enrichment. All data sets within the
collection must correspond to the same organism.
}
  \item{classID2setID}{Optional mapping from identifiers used in classes (either \code{active} and \code{inactive} or in
\code{identifiers}) to identifiers used in \code{refCollection}. If given, must be a 2-column matrix or data frame in
which the first column are identifiers used for classes and the second column identifiers used in the reference
collection. Each class identifier may be present only once but multiple class identifiers may map to a single collection
identifier.}
     
  \item{useGroups}{
     Optional specification of groups within \code{refCollection} to be used.
}
  \item{useBackground}{
     Specification of what to use as the background. The choices are \code{"intersection"} (the
intersection of the given and reference identifiers; \code{"given"} (identifiers given in \code{identifiers} or union
of \code{active} and \code{inactive} above); \code{"reference"} (all identifiers present in the specified
groups in \code{refCollection}, possibly mapped through \code{classID2setID}), and \code{"custom"} (identifiers given in
\code{customBackground}. These are not mapped through \code{classID2setID}.)). 
}
  \item{customBackground}{When \code{useBackground} is \code{"custom"}, this must be a vector of identifiers to be used as
the background set in enrichment calculations. These identifiers should correspond to the ones used for classes
(\code{identifiers} or \code{active} and \code{inactive} as they are not mapped through \code{classID2setID}.
}
  \item{removeMissing}{
     Logical: should genes with missing \code{labels} be removed from the background?
}
  \item{useEvidence}{
     Character vector specifying evidence to be used. To see all available evidence codes, use
knownEvidenceCodes().
}
  \item{removeDuplicatesInDifferentClasses}{Logical: should identifiers that have duplicate entries with
different class labels (other than \code{ignoreLabels} below) be removed?}

  \item{ignoreLabels}{
     Optional specifications of label values that are considered to not be part of any proper group and 
for which enrichment should not be calculated. Note that this specification can also affect enrichment of
other labels if some input identifiers are duplicated. See Details for how duplicated identifiers are
handled and how \code{ignoreLabels} affects the removal or not of duplicated identifiers.
}

  \item{alternative}{
     Character string giving the alternative hypothesis, one of \code{"greater"}, \code{"less"} and \code{"two-sided"}.
Enrichment calculations normally only use the "greater" alternative (i.e., testing for over-enrichment) but sometimes one
may also be interested in deviations from the expected overlap size.}

  \item{nBestDataSets}{
     Number of best gene sets to include in the returned enrichment table and returned gene set details.
}
  \item{threshold}{
     Threshold for the enrichment p-value (possibly corrected for multiple testing) 
     that a gene set must pass to be included in the enrichment table and details.
}
  \item{thresholdType}{
     Multiple testing correction type for the \code{threshold} above. 
}
  \item{getBonferroniCorrection}{
     Logical: should Bonferroni-corrected p-values be returned (along with the nominal p-values)?
}
  \item{getFDR}{
     Logical: should FDR estimates be returned?
}

  \item{getOverlapIDs}{Logical: should identifiers of genes that overlap between each class and gene set be 
included in the output \code{enrichmentTable}?}

  \item{getOverlapSymbols}{Logical: should gene symbols of genes that overlap between each class and gene set be 
included in the output \code{enrichmentTable}? The symbols will be determined from the Entrez identifiers
that correspond to the organism that the collection's data sets correspond to.}

  \item{ID2symbol}{Optional map from identifiers to human-readable symbols or names, e.g., gene symbols. Needed when
\code{getOverlapSymbols} above is \code{TRUE}.}

  \item{maxReportedOverlapGenes}{If either \code{getOverlapIDs} or \code{getOverlapSymbols} are \code{TRUE},
this argument specifies the maximum number overlap genes that will be reported within the enrichment table.
If more than \code{maxReportedOverlapGenes} genes overlap, the entry will read "(More than xx genes overlap)" to
preven extremely long lines in the output table. Users who would like to see the overlap genes in the
enrichment table irrespective of their number should set this to a high-enough number. Note that this
maximum only affects reporting of genes in the component \code{enrichmentTable} of the output; the overlap
genes are reported in completeness in the component \code{dataSetDetails} of the output.}

  \item{entrySeparator}{Character string that will be used to separate entries in multi-entry output, e.g.,
in a list of groups that a set belongs to.}

  \item{groupSeparator}{Character string that will be used to separate group names in lists of groups that
a gene set belongs to. }

  \item{geneSeparator}{Character string that will be used to separate identifiers or symbols in lists of overlap
genes. }

  \item{classColName}{Character string giving the name of the class column in the output enrichment table.}

  \item{getDataSetDetails}{
     Logical: should gene set details be returned?
}
  \item{verbose}{
     Controls level of verbosity. Higher numbers mean more diagnostic output.
}
  \item{indent}{
     Controls indentation. Each unit adds two spaces.
} 
  \item{...}{Other arguments that are currently ignored.}
}
\details{
  Input labels can be specified in two ways. One can supply the \code{identifiers} and class
\code{labels}, leading to a calculation of enrichment of each of the classes in \code{labels} except for
ignored classes specified in \code{ignoreLabels}. If \code{labels} is a matrix, each column is taken to
represent a separate partition into classes. Alternatively, one can specify \code{active} and
\code{inactive} identifiers. Any genes present in both \code{active}
and \code{inactive} are removed from \code{inactive}; thus, one may specify one \code{inactive} for several \code{active}
sets if their background is the same. When multiple \code{inactive} sets are supplied, the "given" background for each
\code{active} set is the union of the set and its corresponding \code{inactive}.

 If input \code{removeDuplicatesInDifferentClasses} is \code{TRUE}, 
input \code{identifiers} are checked for duplicates that are present in different classes in the same column 
supplied in \code{labels}.
(Separate \code{active} sets are never checked for the same identifier occuring in two separate \code{active} sets, and
the checks are also not carried out between different columns of \code{labels}.)
The duplicates are 
handled as follows: if the class labels of all occurrences of a
particular identifier are the same, duplicates are removed and the class label is retained. If the class
labels are the same except for the presence of missing and ignored class labels, the (unique) non-ignored
class label is retained. Otherwise (i.e., if the class labels of a duplicated identifier are not all the
same and not in \code{ignoredLabels}), the class labels for all occurrences of the identifier are set to
missing.

The reference gene sets can be obtained in several ways. First, multiple collections are provided by the
companion package \code{anRichment}. 
Second, users can
create their own reference collection(s) using functions \code{\link{newCollection}}, \code{\link{newGeneSet}},
\code{\link{newGeneProperty}}, \code{\link{newGroup}}, and \code{\link{addToCollection}}. 

One should pay careful attention to the choice of the background (controlled by \code{useBackground}) 
for enrichment. When both the query
gene sets (i.e., \code{identifiers} or \code{active} and \code{inactive}) 
and the reference collection represent a fairly complete gene set for
the organism
(say at least 10,000 genes for human), we recommend \code{"intersection"}. This will lead to least-inflated
enrichment p-values. When either the query or the reference gene set is small but the other one is large,
we recommend using the larger of the two as background (\code{"given"} if the query set is large and
\code{"reference"} if the reference collection is large). If both query and reference sets are small, the last
resort is to use \code{"custom"} which can be, for example, all organismal genes. 
This can lead to highly inflated p-values and so should be used with
caution. Since the evidence for the presence of a particular gene in a given gene set can be based on various levels
of evidence, users can choose to restrict the genes to those with particular evidence types.

The enrichment p-values are calculated using the hypergeometric test and are equivalent to 
p-values obtained by the function \code{\link{fisher.test}}. Bonferroni correction is performed with
respect to the number of comparisons (number of classes times the number of gene sets). The FDR estimates
are obtained using Benjamini-Hochberg procedure as implemented in the \code{\link{p.adjust}} function. Note
that neither the Bonferroni nor FDR corrections take into account the (in practice usually strong) overlaps between the
various gene sets in typical collections and the resulting corrected p-values should be interpreted with caution.
}

\value{
  A list with the following components:
  \item{enrichmentIsValid}{Logical vector with one component per background, 
 indicating whether enrichment statistics are valid. Invalid statistics
result from situations with zero-size background.
}
  \item{enrichmentTable}{A data frame giving the enrichment and other statistics for highest-enriched gene
sets for each input class present in \code{classLabels}. Specifically, for each class, the data frame contains
information on the \code{nBestDataSets} gene sets, plus any additional sets that pass the appropriate
p-value threshold. The individual columns are as follows: \code{class} lists the input class, \code{rank}
the rank of the gene set within each class; \code{dataSetID} and \code{dataSetName} report the gene set ID
and name; \code{inGroups} lists the groups, separated by commas, the gene set belongs to (including groups implied by
parent relationships); columns
\code{pValue}, \code{Bonferroni} and \code{FDR} give the nominal enrichment p-value, Bonferroni
correction, and FDR estimate (\code{Bonferroni} and \code{FDR} are only present if input
\code{getBonferroni} and \code{getFDR} are \code{TRUE}); 
column \code{nCommonGenes} gives the number of genes common between the class
and the gene set, \code{fracOfEffectiveModuleSize} reports the number of common genes as a fraction of the
effective class size (i.e., number of genes in class and in the background); \code{classSize.asGiven} lists the
class sizes as they were input to the function; \code{validClassSize} gives the class size after removal of missing
identifiers and potentially identifiers that map to more than one class; \code{effectiveClassSize} lists the effective class
size; \code{fracOfEffectiveSetSize} gives the number of common genes as a fraction of the effective gene
set size; 
\code{expectedFracOfEffectiveClassSize} gives the expected fraction of the effective set size assuming the
null hypothesis of independence holds,
\code{enrichmentRatio} gives the ration of observed and expected fractions,
\code{effectiveSetSize} reports the effective gene set size (number of genes in the gene set
that are also present in the background). The effective class and set sizes are the sizes that enter the
enrichment calculation and are often smaller than the given sizes because of restriction to a common
background. Finally,  \code{overlapIDs} 
gives the overlap IDs or symbols or both for each comparison as a character string separated by the separator
\code{geneSeparator} (this column is present only if input \code{getOverlapIDs} or \code{getOverlapSymbols} is
\code{TRUE}). At most \code{maxReportedOverlapGenes} is reported here; if the overlap is larger, the field will contain
the character string "(More than xxx overlapping genes)". }

\item{pValues}{ a matrix in which every row corresponds to a gene set, and every column to an input class,
giving the enrichment p-values. }

\item{Bonferroni}{ Only present if input \code{getBonferroni} is \code{TRUE}. A matrix in which every row
corresponds to a gene set, and every column to an input class, giving the Bonferroni-corrected 
enrichment p-values. }

\item{FDR}{ Only present if input \code{getFDR} is \code{TRUE}. A matrix in which every row
corresponds to a gene set, and every column to an input class, giving the FDR estimates for the enrichments. }

\item{countsInDataSet}{ A matrix in which every row
corresponds to a gene set, and every column to an input class, giving the gene counts in the overlap of the
class and the gene set. }

\item{enrichmentRatio}{ A matrix in which every row
corresponds to a gene set, and every column to an input class, giving the ratio of observed vs. expected overlaps.}

\item{effectiveDataSetSizes}{ A (number of reference sets) times (number of backgrounds) matrix 
giving reference set sizes when restricted to the corresponding background.}

\item{effectiveClassSizes}{ A vector with one component per class giving its
size when restricted to the background used in the calculation.}

\item{effectiveClass2bg}{A vector with one component per class giving the index of the background used for enrichment
calculations for that class.}

\item{effectiveBackgroundSize}{ Number of genes in each background. }

\item{dataSetDetails}{A list in which every component corresponds to an input class. Each such component is
in turn a list in which each component corresponds to a gene set. The gene sets listed for each class are
the same as the gene sets reported in \code{enrichmentTable}. These components list more detailed
information about the highly-enriched gene sets, some of which duplicates the information given in
\code{enrichmentTable}. Specifically, the components are \code{dataSetID} (gene set ID); \code{dataSetName}
(gene set name), \code{dataSetDescription} (gene set description), \code{dataSetGroups} (groups
the gene set belongs to, including implied groups), 
\code{enrichmentP} (enrichment p value), \code{commonGeneEntrez} (a vector giving Entrez
identifiers of genes the class and gene set have in common), and \code{commonGenePositions} (positions of
the common genes in the input \code{identifiers} vector).}

\item{identifierInCollection}{Logical vector with one element per input \code{identifier} indicating whether
the identifier is part of at least one gene set in \code{refCollection}}.

\item{effectiveLabels}{A vector with one element per input \code{classLabels}, giving the "effective"
label. This equals the input class label if the corresponding identifier is present only once, or if the
corresponding identifier is present several times but the class for all occurrences is consistent (all
non-ignored class labels are the same). Otherwise, that is, when the corresponding identifier is present
more than once and the corresponding class labels are not consistent, the effective label is set to the
\code{safeToUseMissingLabel} value.}
}
\references{
For FDR claculations, Benjamini, Yoav; Hochberg, Yosef (1995). "Controlling the false discovery rate: a
practical and powerful approach to multiple testing". Journal of the Royal Statistical Society, Series B 57
(1): 289-300.
}
\author{
  Peter Langfelder
}

\seealso{
 \code{\link[anRichment]{buildGOcollection}} for Gene Ontology reference gene sets;
 \code{\link[anRichment]{internalCollection}} for a collection of gene sets from the literature.
}
\keyword{ misc }% __ONLY ONE__ keyword per line
