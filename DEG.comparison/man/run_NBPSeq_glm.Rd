\name{run_NBPSeq_glm}
\alias{run_NBPSeq_glm}
\title{
run_NBPSeq_glm - Negative Binomial (NB) models for two-group comparisons and regression inferences from RNA-Sequencing Data.
}
\description{
For each row of the input data matrix, nb.glm.test fits an NB log-linear regression model and performs large-sample tests for a one-dimensional regression coefficient.
}
\usage{
run_NBPSeq_glm(counts, mycomp3)
}
\arguments{
  \item{counts}{
\code{date.frame} containing raw read counts.}
  \item{mycomp3}{
  \code{list} where comparisons are defined in a list.}
}
\value{
\code{data.frame} containing \code{NBPSeq_glm} results from all comparisons. Comparison labels are appended to column titles for tracking.
}
\references{
Di, Y. et al., The NBP Negative Binomial Model for Assessing Differential Gene Expression from RNA-Seq. Statistical applications in genetics and molecular biology, 10(1), pp.1???28.
}
\author{
Daniela Cassol
}
\seealso{
\code{run_NBPSeq_nbp} and \code{NBPSeq} vignette
}
\examples{
targetspath <- system.file("extdata", "targets.txt", package="DEG.comparison")
targets <- read.delim(targetspath, comment="#")
cmp <- readComp(file=targetspath, format="matrix", delim="-")
countDFeByg <- system.file("extdata", "countDFeByg.xls", package="DEG.comparison")
countDFeByg <- read.delim(countDFeByg, row.names=1)
Comp3 <- list(AP1.4_AP1.67=c("AP1.4A","AP1.4B", "AP1.67A", "AP1.67B"), AP3.4_AP3.67=c("AP3.4A","AP3.4B", "AP3.67A" ,"AP3.67B"), AG.4_AG.67=c("AG.4A", "AG.4B", "AG.67A","AG.67B"))
NBPSeq.glmDF <- run_NBPSeq_glm (countDFeByg, Comp3)
DEG_list_NBPSeq.glmDF <- filterDEGnew(degDF=NBPSeq.glmDF, filter=c(Fold=2, FDR=1), method="NBPSeq.glm")
DEG_list_NBPSeq.glmDF$Summary[1:4,]
}
\keyword{run_DEGs}