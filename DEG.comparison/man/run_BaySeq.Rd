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
\keyword{run_DEGs}
