\name{monodr}
\alias{monodr}
\title{Restoration of monotonicity respecting using adaptive interpolation of
  Huang (2017)}
\description{Restoration of monotonicity respecting for dynamic regression
  models such as quantile regression model, Aalen's additive hazards model,
  and dynamic survival models of Peng and Huang (2007, \emph{Biometrika},
  pp 719--733).}
\usage{monodr(origrc,zch,initau=0.5,taus=numeric(0))}
\arguments{
  \item{origrc}{(piecewise-constant) original regression coefficient.}
  \item{zch}{collection of covariate values for which the monotonicity
   is enforced. In the presence of intercept in "origrc", unity element needs
   to be included. "zch" could be the original covariate matrix. However,
   using its convext hull, if readily available, may be computationally more
   efficient.}
  \item{initau}{starting tau value from which monotonicity-respecting is
  restored.}
  \item{taus}{(increasing) tau values at which quantile coefficient is
  of interest.}
}
\value{
  \item{airc}{regression coefficient after adaptive interpolation: each
  column corresponds to an interpolating point (the intercept is followed by
  slope coefficients, and final element is tau, the probability index.)}
  \item{bt}{adaptively interpolated coefficient at taus, only available if taus
  is specified.}
}
\examples{
## run the example in cequre first
num <- 200
beta <- c(.5, .5)
cvt.1 <- as.numeric(runif(num)<0.5)
cvt.2 <- runif(num)
resid <- rexp(num)
tres <- 1-exp(-resid)
event.t <- log(resid)+beta[1]*cvt.1*ifelse(tres<.4,tres/.4,1)+beta[2]*cvt.2
censr.t <- log(runif(num, 0, 5))
x <- pmin(event.t, censr.t)
dlt <- as.numeric(event.t<=censr.t)
fit <- cequre(x,dlt,cbind(cvt.1,cvt.2),taus=.1*seq(1,7,2),res=200)

## resotration of monotonicity respecting
zch <- matrix(c(1,1,1,1,0,0,1,1,0,1,0,1),ncol=3)
mfit <- monodr(fit$curve,zch,initau=fit$tau.bnd/2,taus=.1*seq(1,7,2))
}
\references{Huang, Y. (2017) Restoration of monotonicity respecting in dynamic regression. \emph{Journal of the American Statistical Association} 112, 613--622.}
