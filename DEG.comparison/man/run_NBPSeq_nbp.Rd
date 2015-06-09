\name{run_NBPSeq_nbp}
\alias{run_NBPSeq_nbp}
\title{
run_NBPSeq_nbp - Negative Binomial (NB) models for two-group comparisons and regression inferences from RNA-Sequencing Data.
}
\description{
nbp.test fits an NBP model to the RNA-Seq counts and performs Robinson and Smyth's exact NB test on each gene to assess differential gene expression between two groups.
}
\usage{
run_NBPSeq_nbp(counts, mycomp3)
}
\arguments{
  \item{counts}{
\code{date.frame} containing raw read counts.}
  \item{mycomp3}{
  \code{list} where comparisons are defined in a list.}
}
\value{
\code{data.frame} containing \code{NBPSeq_nbp} results from all comparisons. Comparison labels are appended to column titles for tracking.
}
\references{
Di, Y. et al., The NBP Negative Binomial Model for Assessing Differential Gene Expression from RNA-Seq. Statistical applications in genetics and molecular biology, 10(1), pp.1???28.
}
\author{
Daniela Cassol
}
\seealso{
\code{run_NBPSeq_glm} and \code{NBPSeq} vignette
}
\examples{
targetspath <- system.file("extdata", "targets.txt", package="DEG.comparison")
targets <- read.delim(targetspath, comment="#")
cmp <- readComp(file=targetspath, format="matrix", delim="-")
countDFeByg <- system.file("extdata", "countDFeByg.xls", package="DEG.comparison")
countDFeByg <- read.delim(countDFeByg, row.names=1)
Comp3 <- list(AP1.4_AP1.67=c("AP1.4A","AP1.4B", "AP1.67A", "AP1.67B"), AP3.4_AP3.67=c("AP3.4A","AP3.4B", "AP3.67A" ,"AP3.67B"), AG.4_AG.67=c("AG.4A", "AG.4B", "AG.67A","AG.67B"))
NBPSeq.nbpDF <- run_NBPSeq_nbp (countDFeByg, Comp3)
DEG_list_NBPSeq.nbpDF <- filterDEGnew(degDF=NBPSeq.nbpDF, filter=c(Fold=2, FDR=1), method="NBPSeq.nbp")
DEG_list_NBPSeq.nbpDF$Summary[1:4,]
}
\keyword{run_DEGs}
