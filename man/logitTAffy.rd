\name{logitTAffy}
\alias{logitTAffy}
\title{ Testing for differential gene expression using the Logit-t algorithm}
\description{
This function takes an instance of AffyBatch and calculates t-statistics 
for tests of differential gene expression for oligonucleotide arrays using the Logit-t
algorithm.
}
\usage{
logitTAffy(object, group)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ an instance of \code{\link[affy:AffyBatch-class]{AffyBatch}} }
  \item{group}{ a vector specifying the group label for each array }
}
\details{ 
For more details see the package vignette.  
}
\value{  
A named vector containing the t-statistics for each probe set for each array.
}
\references{ William J Lemon, Sandya Liyanarachchi and Ming You (2003).  A high performance test of differential gene expression for
oligonucleotide arrays.  Genome Biology 2003, 4:R67. http://genomebiology.com/2003/4/10/R67.}
\author{ Tobias Guennel \email{tguennel@vcu.edu} }
\seealso{ \code{\link[affy:AffyBatch-class]{AffyBatch}}}
\examples{
if(require(SpikeInSubset)){
library(SpikeInSubset)
data(spikein95)
logitTex<-logitTAffy(spikein95, group=c("A","A","A","B","B","B"))
logitTex[1:10]                                                              # extract t-statistics for first ten probe sets
logitTex[grep("AFFX-BioB-5_at",names(logitTex))]                         # extract t-statistics for specific probe set
pvals<-(1-pt(abs(logitTex),df=4))*2                                         # calculate two-sided p-values
signifgenes<-names(logitTex)[pvals<0.01]                                    # find significant probe sets at 0.01 significance level
}else{
stop("Please install the SpikeInSubset package to run the example.")
}
}