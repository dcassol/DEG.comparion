\name{olBarplot}
\alias{olBarplot}
\title{
olBarplot - Define Bar Plot Function
}
\description{
Detailed instructions for using the functions of this script are available on this page:
http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/R_BioCondManual.html#R_graphics_venn
}
\usage{
olBarplot(OLlist = OLlist, mycol = "default", margins = c(6, 10, 3, 2), mincount = 0, mysub = "default", ...)
}
\details{
Utilities: 
(1) Venn Intersects
    Computation of Venn intersects among 2-20 or more sample sets using the typical
    'only in' intersect logic of Venn comparisons, such as: objects present only in 
    set A, objects present only in the intersect of A & B, etc. Due to this restrictive 
    intersect logic, the combined Venn sets contain no duplicates.  
(2) Regular Intersects
    Computation of regular intersects among 2-20 or more sample sets using the
    following intersect logic: objects present in the intersect of A & B, objects present 
    in the intersect of A & B & C, etc. The approach results usually in many duplications 
    of objects among the intersect sets.
(3) Graphical Utilities
    - Venn diagrams of 2-5 sample sets. 
    - Bar plots for the results of Venn intersect and all intersect approaches derived 
      from many samples sets. 

Detailed instructions for using the functions of this script are available on this page:
    http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/R_BioCondManual.html#R_graphics_venn

Revision history:
    March 24, 2012: fixed substring problem in plotVenn function
}
\references{
Thomas Girke (2015). systemPipeR: systemPipeR: NGS workflow and report generation environment. R package version 1.0.12.
https://github.com/tgirke/systemPipeR
}
\author{
Thomas Girke
}
\seealso{
\code{vennPlot, overLapper}}