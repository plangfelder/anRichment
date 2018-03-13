\name{genomicPositionCollection}
\alias{genomicPositionCollection}

\title{
Build genomic position collection for anRichment
}
\description{
Using annotation data from Bioconductor, this function builds 
}
\usage{
genomicPositionCollection(organism, spacings, namePrefix = "Genomic position", dataSource = "", overlapFactor = 2, membershipBy = c("start", "end", "overlap"), useUnits = "Mb", unit = 1e+06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{organism}{
%%     ~~Describe \code{organism} here~~
}
  \item{spacings}{
%%     ~~Describe \code{spacings} here~~
}
  \item{namePrefix}{
%%     ~~Describe \code{namePrefix} here~~
}
  \item{dataSource}{
%%     ~~Describe \code{dataSource} here~~
}
  \item{overlapFactor}{
%%     ~~Describe \code{overlapFactor} here~~
}
  \item{membershipBy}{
%%     ~~Describe \code{membershipBy} here~~
}
  \item{useUnits}{
%%     ~~Describe \code{useUnits} here~~
}
  \item{unit}{
%%     ~~Describe \code{unit} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (organism, spacings, namePrefix = "Genomic position", 
    dataSource = "", overlapFactor = 2, membershipBy = c("start", 
        "end", "overlap"), useUnits = "Mb", unit = 1e+06) 
{
    organism = organismLabels(organism)[1]
    if (!organism \%in\% c("human", "mouse")) 
        stop("The position collection is only available for human and mouse.")
    txNames = switch(organism, human = "TxDb.Hsapiens.UCSC.hg19.knownGene", 
        mouse = "TxDb.Mmusculus.UCSC.mm10.knownGene")
    txdb.lst = as.list(get(txNames))
    geneTxIds = txdb.lst$genes$tx_id
    geneData = txdb.lst$transcripts[match(geneTxIds, txdb.lst$transcripts$tx_id), 
        ]
    chrLevels = sort(unique(as.character(geneData$tx_chrom)))
    membershipBy = match.arg(membershipBy)
    dataSets = list()
    mainGroupName = namePrefix
    groups = list(newGroup(mainGroupName, spaste(namePrefix, 
        ": Genes by genomic location"), dataSource))
    for (chr in chrLevels) {
        printFlush(spaste("Working on chromosome ", chr))
        chrGrpName = spaste(namePrefix, ", chromosome ", chr)
        chrGrp = newGroup(name = chrGrpName, description = spaste(namePrefix, 
            ": Genes on ", chr), source = dataSource)
        keep = which(geneData$tx_chrom == chr)
        geneData1 = geneData[keep, ]
        maxPos = max(geneData1$tx_end)
        for (sp in spacings) {
            nSets = ceiling(overlapFactor * maxPos/sp)
            setStarts = c(0:(nSets - 1)) * sp/overlapFactor
            dataSets1 = lapply(1:nSets, function(i) {
                start1 = (i - 1) * sp/overlapFactor + 1
                end1 = start1 + sp
                if (membershipBy == "start") {
                  keep.int = which(geneData1$tx_start >= start1 & 
                    geneData1$tx_start < end1)
                }
                else if (membershipBy == "end") {
                  keep.int = which(geneData1$tx_end >= start1 & 
                    geneData1$tx_end < end1)
                }
                else keep.int = which(geneData1$tx_end >= start1 & 
                  geneData1$tx_start < end1)
                if (length(useUnits) > 0 && is.character(useUnits)) {
                  prettyStart = spaste(round(start1/unit, 1), 
                    " ", useUnits)
                  prettyEnd = spaste(round(end1/unit, 1), " ", 
                    useUnits)
                  prettySp = spaste(round(sp/unit, 1), " ", useUnits)
                  setName = spaste(namePrefix, ", ", chr, ": ", 
                    prettySp, " at ", prettyStart)
                  shortSetName = spaste(chr, ": ", prettySp, 
                    " starting at ", prettyStart)
                }
                else {
                  prettyStart = formatC(start1, format = "d", 
                    big.mark = ",")
                  prettyEnd = formatC(end1 - 1, format = "d", 
                    big.mark = ",")
                  setName = spaste(namePrefix, ", ", chr, ": ", 
                    prettyStart, "--", prettyEnd)
                  shortSetName = spaste(chr, ": ", prettyStart, 
                    "--", prettyEnd)
                }
                if (length(keep.int) > 0) {
                  newGeneSet(geneEntrez = txdb.lst$genes$gene_id[keep[keep.int]], 
                    geneEvidence = "other", geneSource = dataSource, 
                    ID = spaste(namePrefix, ".", chr, ".", i), 
                    name = setName, shortName = shortSetName, 
                    description = spaste(namePrefix, ": Genes on chromosome ", 
                      shortSetName), source = dataSource, organism = organism, 
                    internalClassification = c(mainGroupName, 
                      chrGrpName, setName), groups = c(mainGroupName, 
                      chrGrpName), externalDB = dataSource)
                }
                else NULL
            })
            dataSets = c(dataSets, dataSets1)
        }
        groups = c(groups, list(chrGrp))
    }
    dataSets = dataSets[sapply(dataSets, length) > 0]
    newCollection(dataSets = dataSets, groups = groups)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
