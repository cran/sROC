% This is a template R document file
% Distributed under GPL 3 or later

\name{bw.CDF}
\alias{bw.CDF}
\title{Bandwidth Selectors for Kernel CDF Estimation}
\description{
  Rule-of-thumb bandwidth selectors for kernel CDF estimation using the normal CDF or PDF reference approach.}
\usage{
 bw.CDF(x, method="npdf")
}
\arguments{
\item{x}{numeric vector.}
\item{method}{either ``npdf'' (the normal PDF reference approach) or ``ncdf'' (the normal CDF reference approach).}
}
\details{
\code{bw.CDF} implements a rule-of-thumb for choosing the bandwidth of a Gaussian kernel CDF estimator. 
}
\value{
A bandwidth on a scale suitable for the \code{bw} argument of \code{\link{kCDF}}.
}
\author{
  X.F. Wang \email{wangx6@ccf.org}
}
\references{
Altman, N., and Leger, C. (1995). Bandwidth selection for kernel distribution function estimation. \emph{Journal of Statistical Planning and Inference}, 46, 195-214.
}
\seealso{
 \code{\link{kCDF}}, \code{\link{bw.CDF.pi}}.
}
\examples{

set.seed(100)
n <- 200
x <- c(rnorm(n/2, mean=-2, sd=1), rnorm(n/2, mean=3, sd=0.8))
bw.CDF(x, method="npdf")
bw.CDF(x, method="ncdf")

}
\keyword{smoothing}

