## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install, eval=FALSE, message=FALSE, warning=FALSE------------------------
#  install.packages("cequre")

## ----cequre, eval=TRUE, message=FALSE, warning=FALSE--------------------------
## Mayo PBC data
library(survival)
pbc_analy <- as.matrix(na.omit(pbc[,c("time","status","age","edema","bili","albumin","protime")]))
# log transformation for time, bili, albumin, and protime
pbc_analy[,c(1,5:7)] <- log(pbc_analy[,c(1,5:7)])
colnames(pbc_analy)[c(1,5:7)] <- paste("log",colnames(pbc_analy)[c(1,5:7)])
# convert status to censoring indicator
pbc_analy[,2] <- pbc_analy[,2]>1

## Censored quantile regression
library(cequre)
taus <- .1*(1:9)
fit <- cequre(pbc_analy[,1],pbc_analy[,2],pbc_analy[,-c(1,2)],
taus=taus,res=200)

## ----cequreplot, eval=TRUE, message=FALSE, warning=FALSE, fig.height = 5, fig.width = 7, dpi = 300, out.width = "80%"----
var.names <- c("baseline",colnames(pbc_analy)[-(1:2)])
par(mfrow=c(3,2),mai=c(.5,.5,.1,.1))
for(k in 1:6) {
  plot(fit$curve[7,],fit$curve[k,],type="s",lwd=2,xlab="",ylab="",
  ylim=c(-1,1)*max(abs(fit$curve[k,])))
  lines(c(0,1),c(0,0),lty=2)
  text(1,max(abs(fit$curve[k,])),var.names[k],adj=c(1,1))
# pointwise 95% CI
  for(i in 1:length(taus)) lines(c(1,1)*taus[i],fit$bt[k,i]+
    c(-1,1)*qnorm(.975)*sqrt(fit$va[k,k,i]))
  }
mtext("Estimated regression coefficient",side=2,line=-1,outer=T,cex=.8)
mtext("Probability",side=1,line=-1,outer=T,cex=.8)

## ----monodr, eval=TRUE, message=FALSE, warning=FALSE, fig.height = 5, fig.width = 7, dpi = 300, out.width = "80%"----
# zch needs to include constant 1 as being the covariate for the intercept
mfit <- monodr(fit$curve,
               zch = cbind(rep(1,dim(pbc_analy)[1]),pbc_analy[,-(1:2)]),
               initau=fit$tau.bnd/2)
# plot the original regression coefficient (in black) and the one after the
# monotonicity-respecting restoration (in orange)
par(mfrow=c(3,2),mai=c(.5,.5,.1,.1))
for(k in 1:6) {
  plot(fit$curve[7,],fit$curve[k,],type="s",xlab="",ylab="",
  ylim=c(-1,1)*max(abs(fit$curve[k,])))
  lines(c(0,1),c(0,0),lty=2)
  text(1,max(abs(fit$curve[k,])),var.names[k],adj=c(1,1))
  lines(mfit$airc[7,],mfit$airc[k,],lwd=2,col="dark orange")
  }
mtext("Estimated regression coefficient",side=2,line=-1,outer=T,cex=.8)
mtext("Probability",side=1,line=-1,outer=T,cex=.8)

