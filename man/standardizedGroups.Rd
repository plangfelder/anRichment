\name{standardizedGroups}
\alias{standardizedGroups}
\alias{standardizedGroups.df}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Retrieve standard groups
}
\description{
These accessor functions retrieve information about standard groups that are often used in anRichment/anRichmentMethods.
}
\usage{
standardizedGroups.df()
standardizedGroups(groups.df = standardizedGroups.df())
}
\arguments{
  \item{groups.df}{
A data frame containing the group information. See \code{\link{groupsFromDataFrame}} for details.
}
}
\value{
A list in which every component is a \code{group}, See \code{\link{newGroup}} for what information is contained in groups. 
}
\author{
Peter Langfelder
}
\seealso{
\code{\link{groupsFromDataFrame}} and \code{\link{collection2dataFrames}} for converting groups to data frames and vice
versa;
\code{\link{newGroup}} for creating new groups
}
\keyword{misc}
\concept{anRichment group}
\concept{anRichment}
% Use only one concept per line.
