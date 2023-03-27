\name{cequre}
\alias{cequre}
\title{Censored quantile regression of Huang (2010)}
\description{Implementation of censored quantile regression of Huang (2010),
  with incorporation of an upper bound related to the identification limit on
  probability scale as described in Huang (2013).}
\usage{cequre(x,dlt,z,epsi=0.05,taus=numeric(0),res=0,
  resam.dist=FALSE,nbps=3*length(x))}
\arguments{
  \item{x}{follow-up time.}
  \item{dlt}{censoring indicator: 1 - event, 0 - censored.}
  \item{z}{matrix of covariates (intercept not included): each column
  corresponds to a covariate.}
  \item{epsi}{parameter for the upper bound related to the identification
  limit on probability scale.}
  \item{taus}{(increasing) tau values at which quantile coefficient is
  of interest.}
  \item{res}{number of resampling iterations for variance estimation:
  res=200 is typically sufficient for variance estimation, but res needs
  to be much larger for confidence band construction.}
   \item{resam.dist}{resampling distribution to be reported or not.}
  \item{nbps}{maximum storage size for quantile coefficient: 3*length(x
  is typically sufficient.}
}
\value{
  \item{curve}{estimated (piecewise-constant) quantile coefficient: each
  column corresponds to a jump point (the intercept is followed by slope
  coefficients, and final element is tau, the probability index.)}
  \item{tau.bnd}{upper bound of tau such that determinant of the at-risk
  matrix (for uncensored observations) is at least
  {epsi}^{# regression coefficients} times the initial value,
  subject to provided storage limit (nbps).}
  \item{bt}{estimated quantile coefficient at taus, only available if taus
  is specified.}
  \item{va}{variance estimate associated with bt, only available if taus
  is specified and res>0. As is resampling based, the variance estimate can be
  slightly different over multiple runs unless seed for the random number
  generator is reset each time.}
  \item{dist}{resampling distribution with res resampled curves: dist[ , ,1]
  through dist[ , ,res], only available if res>0 and resam.dist=TRUE.}
  \item{dist.lgth}{lengths of resampled curves, only available if res>0 and
  resam.dist=TRUE.}
}
\examples{
## simulate a dataset following Scenario 1 of Table 1 in Huang (2010)
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

## run censored quantile regression
fit <- cequre(x,dlt,cbind(cvt.1,cvt.2),taus=.1*seq(1,7,2),res=200)
}

\references{Huang, Y. (2010) Quantile calculus and censored
  regression, \emph{The Annals of Statistics} 38, 1607--1637.

  Huang, Y. (2013) Fast censored linear regression. \emph{Scandinavian Journal
  of Statistics} 40, 789--806.}
