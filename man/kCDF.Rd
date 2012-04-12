% This is a template R document file
% Distributed under GPL 3 or later

\name{kCDF}
\alias{kCDF}
\title{Kernel Estimation for Cumulative Distribution Function}
\description{
  To compute the nonparametric kernel estimate for cumulative distribution function (CDF).
}
\usage{
kCDF(x, bw="pi_ucv", adjust=1, kernel=c("normal", "epanechnikov"), xgrid, 
	ngrid=256, from, to, cut=3, na.rm = FALSE, ...)
}
\arguments{
\item{x}{the data from which the estimate is to be computed.}
\item{bw}{the smoothing bandwidth to be used. \code{bw} can also be a character string giving a rule to choose the bandwidth. See \code{\link{bw.CDF}} and \code{\link{bw.CDF.pi}}. The default used the Altman and Leger's plug-in approach with an unbiased cross-validation pilot bandwidth.} 
\item{adjust}{the parameter for adjusting the bandwidth. The bandwidth used for the estimate is actually \code{adjust*bw}. By default, \eqn{adjust=1}.}	
\item{kernel}{a character string giving the smoothing kernel to be used. This must be either ``normal'' or ``epanechnikov''. By default, the normal kernel is used.}
\item{xgrid}{the user-defined data points at which the CDF is to be evaluated. If missing, the CDF will be evaluated at the equally spaced points defined within the function.}
\item{ngrid}{the number of equally spaced points at which the density is to be estimated.}
\item{from}{the left-most points of the grid at which the density is to be estimated.}
\item{to}{the right-most points of the grid at which the density is to be estimated}
\item{cut}{by default, the values of from and to are cut bandwidths beyond the extremes of the data.}
\item{na.rm}{logical; if \code{TRUE}, missing values are removed from x. If \code{FALSE} any missing values cause an error.}
\item{...}{further arguments for methods.}
}
\details{estimate the nonparametric kernel cumulative distribution function.
 }
\value{
An object of class ``CDF''.
\item{x}{the points where the CDF is estimated.}
\item{Fhat}{the estimated CDF values. These will be numerical numbers between zero and one.}
\item{bw}{the bandwidth used.}
\item{n}{the sample size after elimination of missing values.}
\item{call}{the call which produced the result.}
\item{data}{the original data after elimination of missing values.}
\item{data.name}{the deparsed name of the \code{x} argument.}
\item{has.na}{logical; if \code{TRUE}, there are missing values in the original data.}
The \code{print} method reports \code{\link{summary}} values on the \code{x} and \code{Fhat} components.
}
\author{
  X.F. Wang \email{wangx6@ccf.org}
}
\references{
Nadaraya, E.A. (1964). Some new estimators for distribution functions. \emph{Theory of Probability and its Applications}, 9, 497-500.

Altman, N., and Leger, C. (1995). Bandwidth selection for kernel distribution function estimation. \emph{Journal of Statistical Planning and Inference}, 46, 195-214.
}
\seealso{
 \code{\link{bw.CDF}}, \code{\link{bw.CDF.pi}.}
}
\examples{

## --------------------
set.seed(100)
n <- 200
x <- c(rnorm(n/2, mean=-2, sd=1), rnorm(n/2, mean=3, sd=0.8))
x.CDF <- kCDF(x)
x.CDF
plot(x.CDF, alpha=0.05, main="Kernel estimate of distribution function", CI=FALSE)
curve(pnorm(x, mean=-2, sd=1)/2 + pnorm(x, mean=3, sd=0.8)/2, from =-6, to=6, add=TRUE, lty=2, col="blue")

}
\keyword{nonparametric}
\keyword{smooth}
\keyword{distribution function estimation}


