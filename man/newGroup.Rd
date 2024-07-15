\name{newGroup}
\alias{newGroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create a new group.
}
\description{
  Creates a new group that can be used to group gene sets together.
}
\usage{
newGroup(
  name, 
  description="", 
  source = "", 
  alternateNames = character(0),
  parents = character(0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
     Character string specifying the name of the group.
}
  \item{description}{
    Character string specifying the description of the group.
}
  \item{source}{
     Character vector specifying the source(s) of the group. Can have multiple entries if multiple sources
need to be credited.
}
  \item{alternateNames}{Character vector of alternate names for the group.}

  \item{parents}{Character vector of "parents" or super-groups for the group. Membership in this group
automatically implies membership in all parents.}
}
\details{
  This function creates a simple container with basic description of a group. It does not add any gene
sets to the group. A group itself does not contain references to which gene sets it contains;
rather, when gene sets are created, they can be made part of (possibly mutiple) groups.
}
\value{
  A new group.
}
\author{
  Peter Langfelder
}
\seealso{
   \code{\link{newGeneSet}} for creating new gene sets;

  \code{\link{newGeneProperty}} for creating new gene properties (quantitative or ordinal characteristics 
for genes);

  \code{\link{addDataSetsToGroup}} for tagging data sets with groups.

}
\examples{
  newGroup("Example", "Example group created for the purposes of this illustration");
                
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc}
