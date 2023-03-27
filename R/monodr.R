"monodr" <- function(origrc,zch,initau=0.5,taus=numeric(0))
{ npred <- dim(as.matrix(origrc))[1]-2
  npt <- dim(as.matrix(origrc))[2]
  nzch <- dim(matrix(zch,ncol=npred+1))[1]
  nq <- length(taus)

  fit <- .Fortran("mono", as.double(origrc), as.integer(npred), as.integer(npt),
         as.double(initau), as.double(zch), as.integer(nzch),
         monoindx=integer(npt),nmono=integer(1),
         as.integer(nq), as.double(taus), mu=double((npred+1)*nq),
         double(npred+1),
         PACKAGE="cequre")

  obj <- list(airc=as.matrix(origrc[,fit$monoindx[1:fit$nmono]]))

  if(nq > 0) obj[["bt"]] <- matrix(fit$mu,nrow=npred+1)

  obj
}
