"cequre" <- function(x,dlt,z,epsi=0.05,taus=numeric(0),res=0,
                     resam.dist=FALSE,nbps=3*length(x))
{ size <- length(x)
  npred <- dim(as.matrix(z))[2]
  nq <- length(taus)
  ranexp <- rexp(size*res)

  fit <- .Fortran("cqr", as.double(x), as.integer(dlt), as.double(z),
         double(size), as.integer(size), as.integer(npred),
         para=double((npred+2)*nbps), as.integer(nbps), upar=integer(1),
         as.double(epsi), tauc=double(1),
         integer(npred+1), integer(size), double(size), double(npred+1),
         double(npred+1), double((npred+1)^2), double((npred+1)^2),
         double(npred+1), double(size),
         double((npred+1)^2),double((npred+1)^2),
         as.integer(res), as.double(ranexp), dist=double((npred+2)*nbps*res),
         upars=integer(res),
         as.integer(nq), as.double(taus), mu=double((npred+1)*nq),
         va=double((npred+1)^2*nq),
         PACKAGE="cequre")

  obj <- list(curve=matrix(fit$para,nrow=npred+2)[,1:fit$upar],
              tau.bnd=fit$tauc)

  if(nq > 0)
    {obj[["bt"]] <- matrix(fit$mu,nrow=npred+1)
     if(res > 0) obj[["va"]] <- array(fit$va,dim=c(npred+1,npred+1,nq))}

  if(res > 0 & resam.dist)
    {obj[["dist"]] <-
       array(fit$dist,dim=c(npred+2,nbps,res))[,1:max(fit$upars),]
     obj[["dist.lgth"]] <- fit$upars}

  obj
}
