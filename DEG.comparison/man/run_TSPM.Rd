\name{run_TSPM}
\alias{run_TSPM}
\title{
run_TSPM - "A Two-Stage Poisson Model for Testing RNA-Seq Data"
}
\description{
Simple and powerful statistical approach, based on a two-stage Poisson model, for modeling RNA sequencing data and testing for biologically important changes in gene expression. Users are strongly encouraged to consult the Auer and Doerge (2011) for more detailed information on this topic and how to properly run \code{TSPM} on data sets with more complex experimental designs.
}
\usage{
run_TSPM(counts, mycomp3)
}
\arguments{
  \item{counts}{
\code{date.frame} containing raw read counts.}
  \item{mycomp3}{
  \code{list} where comparisons are defined in a list.}
}
\value{
\code{data.frame} containing \code{TSPM} results from all comparisons. Comparison labels are appended to column titles for tracking.
}
\references{
Paul L. Auer, Rebecca W Doerge: A Two-Stage Poisson Model for Testing RNA-Seq Data. Statistical Applications in Genetics and Molecular Biology 2011, 10(1):26.
}
\author{
Daniela Cassol
}
\seealso{
\code{TSPM}
}
\examples{
targetspath <- system.file("extdata", "targets.txt", package="DEG.comparison")
targets <- read.delim(targetspath, comment="#")
cmp <- readComp(file=targetspath, format="matrix", delim="-")
countDFeByg <- system.file("extdata", "countDFeByg.xls", package="DEG.comparison")
countDFeByg <- read.delim(countDFeByg, row.names=1)
Comp3 <- list(AP1.4_AP1.67=c("AP1.4A","AP1.4B", "AP1.67A", "AP1.67B"), AP3.4_AP3.67=c("AP3.4A","AP3.4B", "AP3.67A" ,"AP3.67B"), AG.4_AG.67=c("AG.4A", "AG.4B", "AG.67A","AG.67B"))
TSPMDF <- run_TSPM(countDFeByg, Comp3)
DEG_list_TSPM <- filterDEGnew(degDF=TSPMDF, filter=c(Fold=2, FDR=1), method="TSPM")
DEG_list_TSPM$Summary[1:4,]
}
\keyword{runDEGs}