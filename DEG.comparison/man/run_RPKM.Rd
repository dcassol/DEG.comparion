\name{run_RPKM}
\alias{run_RPKM}
\title{
run_RPKM
}
\description{
Simple Fold Change Method - RPKM  
}
\usage{
run_RPKM(counts, mycomp1, mycomp2)
}
\arguments{
  \item{counts}{
\code{date.frame} containing raw read counts.
}
  \item{mycomp1}{
code{list} where Factor, Names and groups are defined in a list.}
  \item{mycomp2}{
\code{list} where comparisons are defined in a list.}
}
\value{
\code{data.frame} containing \code{RPKM} results from all comparisons. Comparison labels are appended to column titles for tracking.}
\references{
Mortazavi, A., Williams, B. A., McCue, K., Schaeffer, L., and Wold, B. (2008). Mapping and quantifying mammalian transcriptomes by rna-seq. Nat Methods, 5(7):621-628.}
\author{
Daniela Cassol
}
\seealso{
\code{filterDEG_logFC}
}
\examples{
targetspath <- system.file("extdata", "targets.txt", package="DEG.comparison")
targets <- read.delim(targetspath, comment="#")
cmp <- readComp(file=targetspath, format="matrix", delim="-")
rpkmDFeByg <- system.file("extdata", "rpkmDFeByg.xls", package="DEG.comparison")
rpkmDFeByg <- read.delim(rpkmDFeByg, row.names=1)
#Settings
Comp1 <- list(Factor=(Reduce(union, targets$Factor)), Sample=c(colnames(rpkmDFeByg)), 
              group=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8))

Comp2 <- list(AP1.4_AP1.67=c("AP1.4", "AP1.67"), AP3.4_AP3.67=c("AP3.4", "AP3.67"), 
              AG.4_AG.67=c("AG.4", "AG.67"), AP1.4_AP3.4=c("AP1.4", "AP3.4"), 
              AP1.4_AG.4=c("AP1.4", "AG.4"), AP3.4_AG.4=c("AP3.4", "AG.4"), 
              AP1.67_AP3.67=c("AP1.67", "AP3.67"), AP1.67_AG.67=c("AP1.67", "AG.67"), 
              AP3.67_AG.67=c("AP3.67", "AG.67"))
##Compute mean values for replicates and logFC for comparisons
RPKM_FC <- run_RPKM (rpkmDFeByg, Comp1, Comp2) }
\keyword{runDEGs}