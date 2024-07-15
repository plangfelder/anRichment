\name{checkGroups}
\alias{checkGroups}
\title{
Check groups in a collection for duplicates, groups without gene sets etc.
}
\description{
This function checks groups in a collection for duplicated groups (including groups whose names are very
similar; often happens when different collections are merged) as well as for groups that 
no data sets in the collection belong to (this often happens when a data sets in a collection are subset).
}
\usage{
checkGroups(collection, similarityCut = 0.5, ignore.case = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{collection}{
A collection of data sets.
}
  \item{similarityCut}{
Similarity cutoff for `phrase similarity' 
}
  \item{ignore.case}{
Logical: should case (uppercase vs. lowercase) be ignored when matching patterns? 
}
}
\details{
The function first checks for groups that contain no data sets (taking into account parent-child
relationships between groups) and group membership for data sets that
matches no groups.
If any are found, it tries to find approximate matches; the approximate matches should be
treated as mere suggestions and checked whether they make sense. 
}
\value{
\item{allOK}{Logical; \code{TRUE} if all groups contain at least one data set and all group names metioned in
data set information are also present in groups.}
\item{dataSetGroupsWithoutRecords}{A data frame with two columns, \code{dataSetGroup} giving the group names
from data sets that have no record among the groups in the collection, and \code{closeGroupMatches} giving
the close matches among the names of the groups in the collection.}
\item{groupsWithoutMemberDataSets}{A data frame with two columns, \code{group} giving the name of the groups
without any data set matches, and \code{closeDataSetMatches} giving the group names mentioned in data set
information that are close matches.}
}
\author{
Peter Langfelder
}
\note{
As mentioned above, treat the suggested close matches with suspicion and check them.
}
\keyword{misc}% __ONLY ONE__ keyword per line
