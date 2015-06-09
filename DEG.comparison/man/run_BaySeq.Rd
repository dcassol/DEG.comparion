\name{run_BaySeq}
\alias{run_BaySeq}
\title{
run_BaySeq - empirical Bayesian methods.
}
\description{
BaySeq package identifies differential expression in high-throughput 'count' data, such as that derived from next-generation sequencing machines, calculating estimated posterior likelihoods of differential expression (or more complex hypotheses) via empirical Bayesian methods.
}
\usage{
run_BaySeq(counts, mycomp3, number)
}
\arguments{
  \item{counts}{
      \code{date.frame} containing raw read counts.}
  \item{mycomp3}{
      \code{list} where comparisons are defined in a list.}
  \item{number}{
      \code{number} number the rows in the counts. ex: dim(counts).}
}
\value{
\code{data.frame} containing \code{baySeq} results from all comparisons. Comparison labels are appended to column titles for tracking.
}
\references{
Hardcastle, T.J. & Kelly, K.A., 2010. baySeq: empirical Bayesian methods for identifying differential expression in sequence count data. BMC bioinformatics, 11, p.422.
}
\author{
Daniela Cassol
}
\seealso{
\code{filterDEG_FDR}
}
\examples{
targetspath <- system.file("extdata", "targets.txt", package="DEG.comparison")
targets <- read.delim(targetspath, comment="#")
cmp <- readComp(file=targetspath, format="matrix", delim="-")
countDFeByg <- system.file("extdata", "countDFeByg.xls", package="DEG.comparison")
countDFeByg <- read.delim(countDFeByg, row.names=1)
Comp3 <- list(AP1.4_AP1.67=c("AP1.4A","AP1.4B", "AP1.67A", "AP1.67B"), AP3.4_AP3.67=c("AP3.4A","AP3.4B", "AP3.67A" ,"AP3.67B"), AG.4_AG.67=c("AG.4A", "AG.4B", "AG.67A","AG.67B"))
bayseqDF <- run_BaySeq(countDFeByg, Comp3, number=27416)
DEG_list_bayseqDF <- filterDEG_FDR(degDF=bayseqDF, filter=c(FDR=1), method="BaySeq")
DEG_list_bayseqDF$Summary[1:4,]
}
\keyword{run_DEGs}
