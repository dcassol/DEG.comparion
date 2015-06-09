\name{TSPM}
\alias{TSPM}
\title{
TSPM - "A Two-Stage Poisson Model for Testing RNA-Seq Data"
}
\description{
Simple and powerful statistical approach, based on a two-stage Poisson model, for modeling RNA sequencing data and testing for biologically important changes in gene expression. Users are strongly encouraged to consult the Auer and Doerge (2011) for more detailed information on this topic and how to properly run \code{TSPM} on data sets with more complex experimental designs.
}
\usage{
TSPM(counts, x1, x0, lib.size, alpha.wh = 0.05)
}
\arguments{
  \item{counts}{
\code{date.frame} containing raw read counts
}
  \item{x1}{
\code{x1}a vector of treatment group factors (under the alternative hypothesis)
}
  \item{x0}{
\code{x0} a vector of treatment group factors (under the null hypothesis)
}
  \item{lib.size}{
\code{lib.size} a vector of RNA-Seq library sizes. This could simply be obtained by specifying lib.size <- apply(counts,2,sum). It may also be any other appropriate scaling factor.
}
  \item{alpha.wh}{
\code{alpha.wh} the significance threshold to use for deciding whether a gene is overdispersed. Defaults to 0.05.
}
}
\value{
Returns \code{list} with five components
  \item{log.fold.change}{List of a vector containing the estimated log fold changes for each gene.}
  \item{pvalues}{A vector containing the raw p-values testing differential expression for each gene.}
  \item{index.over.disp}{a vector of integer values containing the indices of the over-dispersed genes.}
  \item{index.not.over.dis}{A vector of integer values containing the indices of the non-over-dispersed genes.}
  \item{padj}{A vector containing the p-values after adjusting for multiple testing using the method of Benjamini-Hochberg.}
}
\references{
Paul L. Auer, Rebecca W Doerge: A Two-Stage Poisson Model for Testing RNA-Seq Data. Statistical Applications in Genetics and Molecular Biology 2011, 10(1):26.
}
\seealso{
\code{run_TSPM}
}
\author{
Paul Auer (plivermo@fhcrc.org) and R.W. Doerge (doerge@purdue.edu)}
\keyword{TSPM}