% This is a template R document file
% Distributed under GPL 3 or later

\name{bw.CDF.pi}
\alias{bw.CDF.pi}
\title{Plug-in Bandwidth Selectors for Kernel CDF Estimation}
\description{
  Plug-in bandwidth selectors for kernel CDF estimation using Altman and Leger's approach.}
\usage{
 bw.CDF.pi(x, pilot="UCV")
}
\arguments{
\item{x}{numeric vector.}
\item{pilot}{ a character string giving a rule to choose the pilot bandwidth to estimate \eqn{E(f''(x))}. 
There are 6 choices: The default, ``ucv'', implement unbiased cross-validation, ``nrd0'' is the Silverman's rule of thumb,
``nrd'' is the Scott's method (1992), ``bcv'' is the biased cross-validation, ``sj'' is  the method of Sheather and Jones (1991), 
``onestage'' is the one-stage method by Wand and Jones.}
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

Scott, D. W. (1992) \emph{Multivariate Density Estimation: Theory, Practice, and Visualization.} Wiley.

Sheather, S. J. and Jones,  M. C. (1991). A reliable data-based bandwidth selection method for kernel density estimation. \emph{Journal of Royal Statistical Society Series B}, 53, 683-690.

Silverman, B. W. (1986) \emph{Density Estimation}. London: Chapman and Hall.

Wand, M.P. and Jones, M.C. (1995) \emph{Kernel Smoothing}. New York: Chapmon and Hall.

}
\seealso{
 \code{\link{kCDF}}, \code{\link{bw.CDF}}.
}
\examples{

set.seed(100)
n <- 200
x <- c(rnorm(n/2, mean=-2, sd=1), rnorm(n/2, mean=3, sd=0.8))
bw.CDF.pi(x)
bw.CDF.pi(x, pilot="nrd0")
bw.CDF.pi(x, pilot="nrd")
bw.CDF.pi(x, pilot="bcv")
bw.CDF.pi(x, pilot="sj")
bw.CDF.pi(x, pilot="onestage")

}
\keyword{smoothing}

