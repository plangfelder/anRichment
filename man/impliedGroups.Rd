\name{impliedGroups}
\alias{impliedGroups}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Groups implied by parent-children relationships.
}
\description{
The function returns, for each input query group, either the groups implied by it (parents or super-groups), or
groups that imply the query group (children or sub-groups).
}
\usage{
impliedGroups(groupList, queryGroups = NULL, includeSelf = TRUE, get = c("parents", "children"))
}
\arguments{
  \item{groupList}{
A list in which every component is a group. 
}
  \item{queryGroups}{
Names of the groups for which the implied groups are to be calculated. If \code{NULL}, implied groups will be
calculated for all groups in \code{groupList}.
}
  \item{includeSelf}{
Logical: should a group be included in th elist of groups implied by it?
}
  \item{get}{
One of \code{c("parents", "children")}. 
}
}
\value{
A list with one component per input \code{queryGroups}; each component is a vector giving the names of the
implied (parent or children) groups.
}
\author{
Peter Langfelder
}
\seealso{
\code{\link{newGroup}} for creating new groups
}
\examples{
# Create three groups A, B, C, where C is a child (subgroup) of B and B is a subgroup of A.
groupLst = list(newGroup("A"),
                newGroup("B", parents = "A"),
                newGroup("C", parents = "B"));

# Get implied groups, including self
impliedGroups(groupLst, get = "parents");
impliedGroups(groupLst, get = "children");
}
\keyword{misc}
