\name{groupsFromDataFrame}
\alias{groupsFromDataFrame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create groups from a data frame with group information.
}
\description{
This function creates a list of groups from a data frame with group information.
}
\usage{
groupsFromDataFrame(
   groupDF = NULL,
   groupDF.nameCol = "Name",
   groupDF.descriptionCol = "Description",
   groupDF.sourceCol = "Source",
   groupDF.alternateNamesCol = "AlternateNames",
   groupDF.parentsCol = "Parents",

   # general options
   ignoreCase = TRUE,
   multiEntrySep = "\\\\|",
   ignoreSpacesAroundSep = TRUE,
   pruneEmptyEntries = TRUE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
   Name of the column in \code{groupDF} that contains the group source (such as article reference or 
website).
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
\value{
A list in which every component is a group. See \code{\link{newGroup}}.
}
\author{
Peter Langfelder
}
\seealso{
\code{\link{newGroup}} for creating a single new group;

\code{\link{collectionFromDataFrames}} for creating a collection from data frames.
}
\keyword{misc }
