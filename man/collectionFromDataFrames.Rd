\name{collectionFromDataFrames}
\alias{collectionFromDataFrames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a collection from a series of data frames
}
\description{
This function creates a collection (an R-object) from a series of appropriately formatted data frames.
}
\usage{
collectionFromDataFrames(

     # Data frame giving the meta-information for gene sets
     geneSetInfoDF = NULL,
     geneSetInfoDF.IDcol = "ID",
     geneSetInfoDF.nameCol = "Name",
     geneSetInfoDF.shortNameCol = "ShortName",
     geneSetInfoDF.descriptionCol = "Description",
     geneSetInfoDF.organismCol = "Organism",
     geneSetInfoDF.sourceCol = "Source",
     geneSetInfoDF.groupsCol = "Groups",
     geneSetInfoDF.internalClassificationCol = "InternalClassification",
     geneSetInfoDF.lastModifiedCol = "LastModified",
     lastModifiedFormat = "\%Y-\%m-\%d",
     geneSetInfoDF.alternateNamesCol = "AlternateNames",
     geneSetInfoDF.externalDBCol = "ExternalDB",
     geneSetInfoDF.externalAccessionCol = "ExternalAccession",
     geneSetInfoDF.webLinkCol = "webLink",

     # Data frame giving the content of gene sets
     geneSetContentDF = NULL,
     geneSetContentDF.nameCol = "GeneSetName",
     geneSetContentDF.entrezCol = "Entrez",
     geneSetContentDF.evidenceCol = "Evidence",
     geneSetContentDF.sourceCol = "Source",

     # Data frame giving the meta-information for gene properties
     propertyInfoDF = NULL,
     propertyInfoDF.IDcol = geneSetInfoDF.IDcol,
     propertyInfoDF.nameCol = geneSetInfoDF.nameCol,
     propertyInfoDF.shortNameCol = geneSetInfoDF.shortNameCol,
     propertyInfoDF.descriptionCol = geneSetInfoDF.descriptionCol,
     propertyInfoDF.organismCol = geneSetInfoDF.organismCol,
     propertyInfoDF.sourceCol = geneSetInfoDF.sourceCol,
     propertyInfoDF.groupsCol = geneSetInfoDF.groupsCol,
     propertyInfoDF.internalClassificationCol = geneSetInfoDF.internalClassificationCol,
     propertyInfoDF.lastModifiedCol = geneSetInfoDF.lastModifiedCol,
     propertyInfoDF.alternateNamesCol = geneSetInfoDF.alternateNamesCol,
     propertyInfoDF.externalDBCol = geneSetInfoDF.externalDBCol,
     propertyInfoDF.externalAccessionCol = geneSetInfoDF.externalAccessionCol,
     propertyInfoDF.webLinkCol = geneSetInfoDF.webLinkCol,
     propertyInfoDF.weightIndexCol = "weightIndex",

     # Data frame containing the actual gene properties
     propertyContentDF = NULL,
     propertyContentDF.entrezCol = NULL,
     propertyContentDF.colnamesMatch = "ID",

     # Data frame containing gene weights corresponding to gene properties
     weightDF = NULL,

     # Data frame describing the groups
     groupDF = NULL,
     groupDF.nameCol = "Name",
     groupDF.descriptionCol = "Description",
     groupDF.sourceCol = "Source",
     groupDF.alternateNamesCol = "AlternateNames",
     groupDF.parentsCol = "Parents",

     # General options
     ignoreCase = TRUE,
     multiEntrySep = "\\\\|",
     ignoreSpacesAroundSep = TRUE,
     pruneEmptyEntries = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneSetInfoDF}{
  Data frame containing meta-information for gene sets. Each row corresponds to one gene set.
}
  \item{geneSetInfoDF.IDcol}{
    Name of the column in \code{geneSetInfoDF} that contains the gene set ID (an arbitrary character string
that should be unique to the gene set).
}
  \item{geneSetInfoDF.nameCol}{
   Name of the column in \code{geneSetInfoDF} that contains the gene set name.
}
  \item{geneSetInfoDF.shortNameCol}{
   Name of the column in \code{geneSetInfoDF} that contains a short version of the gene set name. May be
the same column as \code{geneSetInfoDF.nameCol}.
}
  \item{geneSetInfoDF.descriptionCol}{
   Name of the column in \code{geneSetInfoDF} that contains the description of the gene set.
}
  \item{geneSetInfoDF.organismCol}{
    Name of the column in \code{geneSetInfoDF} that contains the organism for which the gene set is
defined.
}
  \item{geneSetInfoDF.sourceCol}{
     Name of the column in \code{geneSetInfoDF} that contains the source of the gene set (for example, an
article reference or web site).
}
  \item{geneSetInfoDF.groupsCol}{
      Name of the column in \code{geneSetInfoDF} that contains the groups to which the gene set belongs.
Groups should be separated using the separator given in \code{multiEntrySep}.
}
  \item{geneSetInfoDF.internalClassificationCol}{
   Name of the column in \code{geneSetInfoDF} that contains the internal classification keywords for the
gene set. Keywords should be separated using the separator given in \code{multiEntrySep}.
}

  \item{geneSetInfoDF.lastModifiedCol}{Name of the column in \code{geneSetInfoDF} that contains the date of
last modification. This should be a character string in a format specified by \code{lastModifiedFormat}.}

  \item{lastModifiedFormat}{Specification of format of dates in the \code{geneSetInfoDF.lastModifiedCol}
column. See \code{\link{strptime}} for the
description of the format string; the default means <year with century>-<month>-<day>. }

  \item{geneSetInfoDF.alternateNamesCol}{Name of the column in \code{geneSetInfoDF} that contains the
alternate names for each gene set.}

  \item{geneSetInfoDF.externalDBCol}{Name of the column in \code{geneSetInfoDF} that contains an external
database name, if available, for each gene set.}
  \item{geneSetInfoDF.externalAccessionCol}{Name of the column in \code{geneSetInfoDF} that contains the
accession identifier in the external database.}
  \item{geneSetInfoDF.webLinkCol}{Name of the column in \code{geneSetInfoDF} that contains the web link (URL)
of the external database or of the entry corresponding to each gene set.}

  \item{geneSetContentDF}{
Data frame containing the information on which genes make up each gene set. Each row corresponds to a
single gene; a gene can be present on multiple rows, for example because a gene belongs to several gene
sets. A gene can be present multiple times in a single gene set as well, with different evidence codes or
with different sources.
}
  \item{geneSetContentDF.nameCol}{
  Name of the column in \code{geneSetContentDF} that contains the gene set name.
}
  \item{geneSetContentDF.entrezCol}{
  Name of the column in \code{geneSetContentDF} that contains the gene Entrez ID.
}
  \item{geneSetContentDF.evidenceCol}{
  Name of the column in \code{geneSetContentDF} that contains the evidence codes.
}
  \item{geneSetContentDF.sourceCol}{
  Name of the column in \code{geneSetContentDF} that contains the source (for example, article reference or web
site)  that claims the gene should be part of the gene set.
}
  \item{propertyInfoDF}{
  Data frame containing meta-information for gene properties. Each row of this data frame corresponds to a
gene property.
}
  \item{propertyInfoDF.IDcol}{
  Name of the column in \code{propertyInfoDF} that contains the gene property ID (an arbitrary character
string that should be unique to the gene property).
}
  \item{propertyInfoDF.nameCol}{
  Name of the column in \code{propertyInfoDF} that contains the gene property name.
}
  \item{propertyInfoDF.shortNameCol}{
  Name of the column in \code{propertyInfoDF} that contains a short version of the name. Short version
of names can be handy in display but can be the same as \code{propertyInfoDF.nameCol}.
}
  \item{propertyInfoDF.descriptionCol}{
  Name of the column in \code{propertyInfoDF} that contains the description of the gene property.
}
  \item{propertyInfoDF.organismCol}{
  Name of the column in \code{propertyInfoDF} that contains the organism  for which the gene set is
defined.
}
  \item{propertyInfoDF.sourceCol}{
     Name of the column in \code{geneSetInfoDF} that contains the source of the gene property (for example, an
article reference or web site).
}
  \item{propertyInfoDF.groupsCol}{
      Name of the column in \code{geneSetInfoDF} that contains the groups to which the gene property belongs.
Groups should be separated using the separator given in \code{multiEntrySep}.
}
  \item{propertyInfoDF.internalClassificationCol}{
   Name of the column in \code{geneSetInfoDF} that contains the internal classification keywords for the
gene set. Keywords should be separated using the separator given in \code{multiEntrySep}.
}

  \item{propertyInfoDF.lastModifiedCol}{Name of the column in \code{geneSetInfoDF} that contains the date of
last modification of the property.}

  \item{propertyInfoDF.alternateNamesCol}{Name of the column in \code{geneSetInfoDF} that contains the
alternate names for each gene property.}

  \item{propertyInfoDF.externalDBCol}{Name of the column in \code{geneSetInfoDF} that contains an external
database name, if available, for each gene property.} 

  \item{propertyInfoDF.externalAccessionCol}{Name of the column in \code{geneSetInfoDF} that contains the
accession identifier in the external database that corresponds to the gene property.}  

  \item{propertyInfoDF.webLinkCol}{Name of the column in \code{geneSetInfoDF} that contains the web link (URL)
of the external database or of the entry corresponding to each gene set.}

  \item{propertyInfoDF.weightIndexCol}{
   Name of the column in \code{geneSetInfoDF} that contains the column in the weight data frame
\code{weightDF} containing weights corresponding to this property.
}
  \item{propertyContentDF}{
   Data frame giving the actual values of the gene properties for each gene. Each row corresponds to one
gene. One of the columns (usually the first) should give the Entrez ID of the gene; other columns are
assumed to correspond each to a different gene property. Gene properties can be continuous or discrete.
}
  \item{propertyContentDF.entrezCol}{
   Name of the column in \code{propertyContentDF} that contains the Entrez ID for each gene.
}
  \item{propertyContentDF.colnamesMatch}{
   Name of the column in \code{geneSetInfoDF} whose content matches column names in \code{propertyContentDF}.
This is used to connect the actual gene property values with the meta-information.
}
  \item{weightDF}{
   A data frame giving the weights of gene properties corresponding to each gene. Each row corresponds to a
gene; the number of rows and the order is the same as in \code{propertyContentDF}. Each column hold a vector
of weights; which weight column corresponds to a property is encoded in the component
\code{propertyInfoDF.weightIndexCol} of the data frame \code{propertyInfoDF}. 
}
  \item{groupDF}{
   A data frame giving the information about groups.
}
  \item{groupDF.nameCol}{
   Name of the column in \code{groupDF} that contains the group name.
}
  \item{groupDF.descriptionCol}{
   Name of the column in \code{groupDF} that contains the group description.
}
  \item{groupDF.sourceCol}{
   Name of the column in \code{groupDF} that contains the group source (such as article reference or web
site).
}
  \item{groupDF.alternateNamesCol}{
   Name of the column in \code{groupDF} that contains alternate names.
}
  \item{groupDF.parentsCol}{
   Name of the column in \code{groupDF} that contains group "parents" or super-groups.
}
  \item{ignoreCase}{
  Logical: should case be ignored when matching names?
}
  \item{multiEntrySep}{
  A character string giving the separator of multiple values in a single character string.
}
  \item{ignoreSpacesAroundSep}{
  Logical: should spaces immediately preceding and following the \code{multiEntrySep} be stripped out of
the values?
}
  \item{pruneEmptyEntries}{
  Logical: should empty values be removed?
}
}
\details{
The aim of this function is to simplify converting gene sets encoded in tables (data frames) to the
internal collection format. 
}
\value{
A collection containing the gene sets, gene properties, and groups as described by the input data frames.

}
\author{
Peter Langfelder
}
\seealso{
\code{\link{newCollection}}, \code{\link{newGeneSet}}, \code{\link{newGeneProperty}},
\code{\link{newGroup}} for functions that create new collections, gene sets, gene properties, and groups,
respectively.

}
\keyword{ misc }% __ONLY ONE__ keyword per line
