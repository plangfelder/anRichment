\name{newGeneProperty}
\alias{newGeneProperty}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a new gene property
}
\description{
This function creates a new numeric or ordinal gene property (characteristic) and adds it to a given
collection, or creates a new collection containing this property. 
}
\usage{
newGeneProperty(
   identifiers, properties, 
   weights = NULL, 
   discreteIndex = NULL, maxDiscreteLevels = 5, 
   propertyNames = colnames(properties), IDs, names, shortNames, descriptions, 
   source, organism, 
   internalClassificationBase, internalClassificationExtras, 
   groupBase, groupExtras, 
   lastModified = Sys.Date(), format = "\%Y-\%m-\%d",
   alternateNames = list(character(0)),
   externalDB = "",
   externalAccession = "",
   webLink = "",
   newGroups = NULL, 
   collection = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{identifiers}{
Gene identifiers (usually Entrez IDs) corresponding to the rows in \code{properties}.
}
  \item{properties}{
A vector (or data frame) in which each element (or row) corresponds to a gene identified in
\code{identifiers} and gives a numeric or discrete value of a variable for the gene. If a data frame,
each column is taken as an independent variable. A vector is treated as a 1-column data frame.
}
  \item{weights}{
If non-NULL, a single vector giving weights for each gene. The default NULL means all genes have equal
weights.
}
  \item{discreteIndex}{ Logical or numeric vector giving the columns 
in \code{properties} that are to be treated as discrete variables. If
\code{NULL}, discrete variables will be determined automatically (see details). To treat all variables as
continuous, use either a vector with all elements equal \code{FALSE} or \code{discreteIndex=numeric(0)}.
}
  \item{maxDiscreteLevels}{
Only used if \code{discreteIndex} is \code{NULL}. 
Integer giving the maximum number of unique values in a column of \code{properties} for that column to be
considered discrete.
}
  \item{propertyNames}{
Names for the properties. If not given, the names will be taken from the columns of \code{properties}. See
\code{IDs} and \code{names} below.
}
  \item{IDs}{ Specifies the identifiers for the newly created properties.
Either a single character string or a vector in which each component corresponds to a column of
\code{properties}. If a single character string, unique identifiers 
are created by concatenating \code{IDs} with \code{propertyNames} (with '.' (dot) as the separator). If
\code{IDs} a character vector already,
its elements will be used as property identifiers.
}
  \item{names}{
Either a single character string or a vector in which each component corresponds to a column of
\code{properties}. If a single character string, unique names are created by concatenating \code{names}
with \code{propertyNames} (with '.' (dot) as the separator). If \code{names} is a character vector already,
its elements will be used as property names.
}
  \item{shortNames}{Specification of short names for each property: 
Either a single character string or a vector in which each component corresponds to a column of
\code{properties}. If a single character string, unique short names are created by concatenating \code{shortNames}
with \code{propertyNames} (with '.' (dot) as the separator). If \code{shortNames} is a character vector already,
its elements will be used as property short names.
}

  \item{descriptions}{
Either a single character string or a vector in which each component corresponds to a column of
\code{properties}. If a single character string, descriptions are created by concatenating
\code{propertyNames} in parenthesis '()'
with \code{description}. If \code{description} is a character vector already,
its elements will be used as property names.

}
  \item{source}{
A character string describing the source (e.g., an article) of the data. The same source will be used for
all given properties.
}
  \item{organism}{
A character string describing the organism. See \code{organismLabels} for recognized organisms.
}
  \item{internalClassificationBase}{
A character vector giving the classification entries that are common to all columns in \code{properties}.
}
  \item{internalClassificationExtras}{
A vector or list with one component per column of \code{properties} giving additional classification
entries for each column.
}
  \item{groupBase}{
A character vector giving the groups that are common to all columns in \code{properties}.
}
  \item{groupExtras}{
A vector or list with one component per column of \code{properties} giving additional group entries for each
column.
}
  \item{lastModified}{Date of last modification of the gene set, given
either as a \code{\link{Date}} object, or as a character string specifying date. If a character string, its
format can be specified using the argument \code{format} below.
}
  \item{format}{Specification of format of character dates. See \code{\link{strptime}} for the
description of the format string; the default means <year with century>-<month>-<day>. }

  \item{alternateNames}{Optional specification of alternate names for each property. A list of character
vectors in which each component corresponds to one property and gives the alternate names for it.}

  \item{ externalDB}{Character vector specifying external database(s) from which the properies were sourced. }
  \item{externalAccession}{Character vector specifying accession numbers in the external databases above that
corresponds to the properties.}
  \item{webLink}{Character vector specifying web links to the external database (or directly to the entries
corresponding to each property.}

  \item{newGroups}{
An optional list in which each element is a group structure (such as the one returned by
\code{\link{newGroup}}). The groups will be added to the returned collection.
}
  \item{collection}{
Optional collection to which the \code{properties} with all supporting information should be added. If not
given, a new collection will be created from the supplied information.
}
}
\details{
A gene property is a numeric (i.e., continuous) or discrete (ordinal) variable with a unique value for each
gene. Additionally, for each gene a weight (e.g. giving confidence in the recorded value) can be given (the
default is equal weights). Each property is identified by a name and a unique identifier, and carries
additional information giving the description, source (for example, an article citation), internal
classification and names of the groups the property belongs to.

To make data storage more efficient, gene properties must be encapsulated in a collection.
}
\value{
A collection containing the given properties (and the content of the given input \code{collection}).
}
\author{
Peter Langfelder
}

\seealso{
\code{\link{newGeneSet}} for creating gene sets;
\code{\link{newGroup}} for creating new groups;
\code{\link{newCollection}} for creating new collections;
\code{\link{mergeCollections}} for merging existing collections;
\code{\link{addToCollection}} for adding gene sets, gene properties, groups, or other collections to an
existing collection.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}% __ONLY ONE__ keyword per line
