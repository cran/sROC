% This is a template R document file
% Distributed under GPL 3 or later

\name{kROC}
\alias{kROC}
\title{Kernel Estimation for ROC Curves}
\description{
  To compute the nonparametric kernel estimate of receiver operating characteristic (ROC) Curves for continuous data.
}
\usage{
kROC(x, y, bw.x="pi_ucv", bw.y="pi_ucv", adjust=1, kernel=c("normal", "epanechnikov"), xgrid, 
	ngrid=256, from, to, cut=3, na.rm = FALSE, ...)
}
\arguments{
\item{x}{numeric vector.}
\item{y}{numeric vector.}
\item{bw.x}{the smoothing bandwidth of \code{x} to be used. \code{bw} can also be a character string giving a rule to choose the bandwidth. See \code{\link{bw.CDF}} and \code{\link{bw.CDF.pi}}. The default used the Altman and Leger's plug-in approach with an unbiased cross-validation pilot bandwidth.} 
\item{bw.y}{the smoothing bandwidth of \code{y} to be used.}
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
\details{ estimate the nonparametric kernel estimate of receiver operating characteristic (ROC) Curves for continuous data
 }
\value{
An object of class ``ROC''.
\item{FPR}{the false positive rate.}
\item{TPR}{the true positive rate.}
\item{bw.x, bw.y}{the bandwidths used.}
\item{nx, ny}{the sample sizes after elimination of missing values.}
\item{call}{the call which produced the result.}
\item{x.data.name, y.data.name}{the deparsed names of the \code{x} argument.}
\item{x.has.na, y.has.na}{logical; if \code{TRUE}, there are missing values in the original data.}
The \code{print} method reports \code{\link{summary}} values on the \code{x} and \code{Fhat} components.
}
\author{
  X.F. Wang \email{wangx6@ccf.org}
}
\references{
Lloyd, C.J. (1998). Using smoothed receiver operating characteristic curves to summarize and compare diagnostic systems. \emph{Journal of the American Statistical Association}, 93(444): 1356-1364.

Zhou, X.H. and Harezlak, J. (2002). Comparison of bandwidth selection methods for kernel smoothing of ROC curves. \emph{Statistics in Medicine}, 21, 2045-2055.

Zou, K.H., Hall, W.J., and Shapiro, D.E. (1997). Smooth non-parametric receiver operating characteristic (ROC) curves for continuous diagnostic tests. \emph{Statistics in medicine}, 16(19): 2143-56.	

}
\seealso{
 \code{\link{bw.CDF}}, \code{\link{bw.CDF.pi}.}
}
\examples{

## --------------------
set.seed(100)
n <- 200
x <- rgamma(n,2,1)
y <- rnorm(n)

xy.ROC <- kROC(x,y, bw.x="pi_sj",bw.y="pi_sj")
xy.ROC

plot(xy.ROC)

}
\keyword{nonparametric}
\keyword{smooth}


