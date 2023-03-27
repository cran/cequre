      subroutine cqr(x,dlt,z,wt,size,npred,
     &     para,npar,upar,epsi,tauc,
     &     intp,lzr,rsdl,rfrac,zbar,proj,uppm,drct,lp,atrisk,
     &     work,
     &     res,ranexp,dist,upars,
     &     nq,taus,mu,va)
      integer size,npred,dlt(size),npar,upar
      double precision x(size),z(size,npred),wt(size),
     &     para(0:npred+1,npar)
      integer intp(npred+1),lzr(size)
      double precision rsdl(size),epsi,tauc,rfrac(npred+1),
     &     zbar(0:npred),proj(0:npred,npred+1),
     &     uppm(npred+1,npred+1),drct(0:npred),lp(size),
     &     atrisk(0:npred,0:npred),work(0:npred,0:npred)
      integer res,upars(res)
      double precision ranexp(size,res),
     &     dist(0:npred+1,npar,res)
      integer nq
      double precision taus(nq),mu(0:npred,nq),va(0:npred,0:npred,nq)
      integer i,j,k,p,q
      double precision dp1,dp2,dp3
      do i=1,size
         wt(i)=1.0d0
      enddo
      call aqm(x,dlt,z,wt,size,npred,
     &     para,npar,upar,
     &     .true.,epsi,tauc,
     &     intp,lzr,rsdl,rfrac,zbar,proj,uppm,drct,lp,atrisk,
     &     work)
      do k=1,res
         call aqm(x,dlt,z,ranexp(1,k),size,npred,
     &        dist(0,1,k),npar,upars(k),
     &        .false.,epsi,tauc,
     &        intp,lzr,rsdl,rfrac,zbar,proj,uppm,drct,lp,atrisk,
     &        work)
      enddo
      do j=1,nq
         do p=0,npred
            mu(p,j)=0.0d0
            do q=0,p
               va(p,q,j)=0.0d0
            enddo
         enddo
      enddo
      do k=1,res
         i=1
         do j=1,nq
            do while(i .le. upars(k) .and.
     &           dist(npred+1,i,k) .le. taus(j))
               i=i+1
            enddo
            i=max0(i-1,1)
            do p=0,npred
               mu(p,j)=mu(p,j)+dist(p,i,k)
               do q=0,p
                  va(p,q,j)=va(p,q,j)+dist(p,i,k)*dist(q,i,k)
               enddo
            enddo
         enddo
      enddo
      if(res .gt. 0) then
         do j=1,nq
            do p=0,npred
               mu(p,j)=mu(p,j)/dble(res)
            enddo
            do p=0,npred
               do q=0,p
                  va(p,q,j)=va(p,q,j)/dble(res)-mu(p,j)*mu(q,j)
               enddo
            enddo
            do p=0,npred
               do q=p+1,npred
                  va(p,q,j)=va(q,p,j)
               enddo
            enddo
         enddo
      endif

      do j=1,nq
         do p=0,npred
            mu(p,j)=0.0d0
         enddo
      enddo
      i=1
      do j=1,nq
         do while(i .le. upar .and.
     &        para(npred+1,i) .le. taus(j))
            i=i+1
         enddo
         i=max0(i-1,1)
         do p=0,npred
            mu(p,j)=para(p,i)
         enddo
      enddo

      end
      subroutine aqm(x,dlt,z,wt,size,npred,
     &     para,npar,upar,
     &     tfil,epsi,tauc,
     &     intp,lzr,rsdl,rfrac,zbar,proj,uppm,drct,lp,atrisk,
     &     work)
      logical tfil
      integer size,npred,dlt(size),npar,upar
      double precision x(size),z(size,npred),wt(size),
     &     para(0:npred+1,npar)
      integer intp(npred+1),lzr(size)
      double precision rsdl(size),epsi,tauc,rfrac(npred+1),
     &     zbar(0:npred),proj(0:npred,npred+1),
     &     uppm(npred+1,npred+1),drct(0:npred),lp(size),
     &     atrisk(0:npred,0:npred),work(0:npred,0:npred)

      double precision logdet

      integer nintp,eintp,bd,ninf

      logical cont,uprsk
      integer i,j,p,q
      double precision irisk,dtmp

      para(0,1)=x(1)
      do i=2,size
         para(0,1)=dmin1(para(0,1),x(i))
      enddo
      do p=1,npred
         para(p,1)=0.0d0
      enddo
      do i=1,size
         rsdl(i)=x(i)-para(0,1)
         if(rsdl(i) .gt. 1.0d-10) then
            lzr(i)=2
         else
            lzr(i)=1
         endif
      enddo
      do p=0,npred
         zbar(p)=0.0d0
      enddo
      do i=1,size
         zbar(0)=zbar(0)+wt(i)
         do p=1,npred
            zbar(p)=zbar(p)+z(i,p)*wt(i)
         enddo
      enddo
      nintp=0
      uprsk=tfil
      if(tfil) then
         do p=0,npred
            do q=0,p
               atrisk(p,q)=0.0d0
            enddo
         enddo
         do i=1,size
            if(dlt(i) .eq. 1) then
               p=0
               q=0
               atrisk(p,q)=atrisk(p,q)+wt(i)
               do p=1,npred
                  q=0
                  atrisk(p,q)=atrisk(p,q)+z(i,p)*wt(i)
                  do q=1,p
                     atrisk(p,q)=atrisk(p,q)+z(i,p)*z(i,q)*wt(i)
                  enddo
               enddo
            endif
         enddo
         irisk=logdet(atrisk,npred+1,work,ninf)
         if(ninf .eq. 1) then
            tauc=0.0d0
            uprsk=.false.
         endif
      endif
      upar=0
      cont=.true.

      do while(cont .and. upar .lt. npar)
         upar=upar+1

         if(upar .eq. 1) then
            para(npred+1,1)=0.0d0
         else
            eintp=0
            do i=1,nintp
               eintp=eintp+dlt(intp(i))
            enddo

            if(eintp .eq. 0) then
               para(npred+1,upar)=1.0d0
            else
               call incstep(dlt,z,wt,size,npred,zbar,
     &              intp,rfrac,eintp,nintp,lzr,
     &              para(npred+1,upar),proj,uppm,
     &              uprsk,atrisk,drct)
            endif

            if(upar .gt. 2) then
               para(npred+1,upar)=1.0d0
     &              -(1.0d0-para(npred+1,upar))*
     &              (1.0d0-para(npred+1,upar-1))
               if(bd .eq. 0) then
                  upar=upar-1
                  para(npred+1,upar)=para(npred+1,upar+1)
               endif
            endif

            if(para(npred+1,upar) .gt. 1.0d0-1.0d-10)
     &           cont=.false.

            if(uprsk) then
               dtmp=logdet(atrisk,npred+1,work,ninf)-irisk
               if(ninf .eq. 1 .or.
     &              dtmp .le. dlog(epsi)*dble(npred+1)) then
                  tauc=para(npred+1,upar)
                  uprsk=.false.
               endif
            endif

            do p=0,npred
               para(p,upar)=para(p,upar-1)
            enddo
         endif

         if(cont) call minstep(rsdl,dlt,z,wt,size,npred,zbar,
     &        intp,rfrac,nintp,lzr,para(0,upar),bd,
     &        proj,uppm,drct,lp)

      enddo

      if(tfil .and. upar .gt. 0 .and. tauc .lt. -5.0d-1)
     &     tauc=para(npred+1,upar)

      end
      double precision function logdet(atrisk,dim,mat,ninf)
      integer dim,ninf
      double precision atrisk(dim,dim),mat(dim,dim)
      integer i,j,k

      do i=1,dim
         do j=1,i
            mat(i,j)=atrisk(i,j)
         enddo
      enddo
      do i=1,dim
         do j=i+1,dim
            mat(i,j)=mat(j,i)
         enddo
      enddo

      ninf=0
      i=1
      do while(i .le. dim)
         if(i .ge. 2) then
            do j=1,dim
               do k=1,min0(i,j)-1
                  mat(i,j)=mat(i,j)-mat(i,k)*mat(k,j)
               enddo
               if(j .lt. i) mat(i,j)=mat(i,j)/mat(j,j)
            enddo
         endif
         if(mat(i,i) .lt. 1.0d-10) then
            ninf=1
            return
         endif
         i=i+1
      enddo

      logdet=dlog(mat(1,1))
      do i=2,dim
         logdet=logdet+dlog(mat(i,i))
      enddo

      end
      subroutine incstep(dlt,z,wt,size,npred,zbar,
     &     intp,rfrac,eintp,nintp,lzr,inc,
     &     proj,uppm,uprsk,atrisk,drct)
      logical uprsk
      integer size,npred,dlt(size),intp(npred+1),eintp,
     &     nintp,lzr(size)
      double precision z(size,npred),wt(size),rfrac(npred+1),
     &     zbar(0:npred),inc
      double precision proj(0:npred,npred+1),
     &     uppm(npred+1,npred+1),atrisk(0:npred,0:npred),
     &     drct(npred+1)
      logical first
      integer i,j,p,q,k,l
      double precision dtmp
      do i=1,nintp
         drct(i)=0.0d0
         do p=0,npred
            drct(i)=drct(i)+zbar(p)*proj(p,i)
         enddo
      enddo
      k=0
      i=nintp
      do while(k .lt. eintp)
         do j=i+1,nintp
            drct(i)=drct(i)-uppm(i,j)*drct(j)
         enddo
         drct(i)=drct(i)/uppm(i,i)
         k=k+dlt(intp(i))
         i=i-1
      enddo
      inc=1.0d0
      i=1
      do while(i .le. nintp)
         if(dlt(intp(i)) .eq. 1 .and. dabs(drct(i)) .gt. 1.0d-10) then
            if(drct(i) .gt. 0.0d0) then
               dtmp=rfrac(i)*wt(intp(i))/drct(i)
            else
               dtmp=(rfrac(i)-1.0d0)*wt(intp(i))/drct(i)
            endif
            if(dtmp .gt. 1.0d-50) then
               inc=dmin1(inc,dtmp)
            else
               do j=i+1,nintp
                  intp(j-1)=intp(j)
                  rfrac(j-1)=rfrac(j)
                  drct(j-1)=drct(j)
               enddo
               nintp=nintp-1
               i=i-1
            endif
         endif
         i=i+1
      enddo
      do p=0,npred
         zbar(p)=zbar(p)*(1.0d0-inc)
      enddo
      first=.true.
      i=1
      l=nintp+1
      do while(i .le. nintp)
         if(dlt(intp(i)) .eq. 1) then
            rfrac(i)=rfrac(i)-drct(i)/wt(intp(i))*inc
            if(uprsk) then
               p=0
               q=0
               atrisk(p,q)=atrisk(p,q)-drct(i)*inc
               do p=1,npred
                  q=0
                  atrisk(p,q)=atrisk(p,q)-z(intp(i),p)
     &                 *drct(i)*inc
                  do q=1,p
                     atrisk(p,q)=atrisk(p,q)-z(intp(i),p)
     &                    *z(intp(i),q)*drct(i)*inc
                  enddo
               enddo
            endif
            if(rfrac(i) .lt. 1.0d-10 .or.
     &           rfrac(i) .gt. 1.0d0-1.0d-10) then
               if(first) l=i
               first=.false.
               lzr(intp(i))=1
               if(rfrac(i) .lt. 1.0d-10) lzr(intp(i))=-1
               do j=i,nintp-1
                  drct(j)=drct(j+1)
                  intp(j)=intp(j+1)
                  rfrac(j)=rfrac(j+1)
               enddo
               nintp=nintp-1
               i=i-1
            endif
         endif
         i=i+1
      enddo

      i=1
      j=nintp
      first=.true.
      do while(i .lt. j)
         if(dlt(intp(i)) .eq. 0) then
            do while(j .gt. i .and. dlt(intp(j)) .eq. 0)
               j=j-1
            enddo
            if(i .lt. j) then
               if(first) then
                  first=.false.
                  if(i .lt. l) l=i
               endif
               k=intp(i)
               intp(i)=intp(j)
               intp(j)=k
               dtmp=rfrac(i)
               rfrac(i)=rfrac(j)
               rfrac(j)=dtmp
               i=i+1
               j=j-1
            endif
         else
            i=i+1
         endif
      enddo

      call orth(z,size,npred,intp,proj,uppm,l,nintp)

      end
      subroutine minstep(rsdl,dlt,z,wt,size,npred,zbar,
     &     intp,rfrac,nintp,lzr,beta,bd,proj,uppm,
     &     drct,lp)
      integer size,npred,dlt(size),intp(npred+1),nintp,
     &     lzr(size),bd
      double precision rsdl(size),z(size,npred),wt(size),
     &     zbar(0:npred),rfrac(npred+1),beta(0:npred)
      double precision proj(0:npred,npred+1),
     &     uppm(npred+1,npred+1),drct(0:npred),lp(size)
      integer ocode(2),i,p,mark,cyc,nconst

      bd=0
      mark=nintp
      do while(mark .gt. 0 .and. dlt(intp(mark)) .eq. 0)
         mark=mark-1
      enddo

      cyc=npred+1
      nconst=0
      do while(nconst .lt. cyc-mark)
         if(nintp .eq. cyc) then
            lzr(intp(mark+1))=1
            if(dlt(intp(mark+1)) .eq. 0) then
               zbar(0)=zbar(0)+wt(intp(mark+1))
     &              *(1-rfrac(mark+1))
               do p=1,npred
                  zbar(p)=zbar(p)+z(intp(mark+1),p)
     &                 *wt(intp(mark+1))*(1-rfrac(mark+1))
               enddo
            else
               if(rfrac(mark+1) .lt. 1.0d-10)
     &              lzr(intp(mark+1))=-1
            endif
            do i=mark+1,nintp-1
               intp(i)=intp(i+1)
               rfrac(i)=rfrac(i+1)
            enddo
            call orth(z,size,npred,intp,proj,uppm,mark+1,cyc-1)
            nintp=nintp-1
         endif

         call line(rsdl,dlt,z,wt,size,npred,proj,zbar,
     &        intp,rfrac,nintp,lzr,beta,ocode,drct,lp)
         if(ocode(1) .eq. 1) bd=1
         if(ocode(2) .eq. 0) then
            cyc=nintp
         else
            call orth(z,size,npred,intp,proj,uppm,nintp,nintp)
         endif

         if(ocode(1) .eq. 0) then
            nconst=nconst+1
         else
            nconst=0
            if(ocode(2) .eq. 1 .and. cyc .eq. nintp) nconst=1
         endif

         if(cyc .lt. npred+1 .and. cyc .gt. mark .and.
     &        nconst .eq. cyc-mark) then
            intp(nintp+1)=-1
            call line(rsdl,dlt,z,wt,size,npred,proj,zbar,
     &           intp,rfrac,nintp,lzr,beta,ocode,drct,lp)
            if(ocode(1) .eq. 1) then
               bd=1
               if(ocode(2) .eq. 1) then
                  cyc=npred+1
                  call orth(z,size,npred,intp,proj,uppm,nintp,nintp)
                  nconst=0
                  if(cyc .eq. nintp) nconst=1
               endif
            else
               if(ocode(2) .eq. 1) nintp=nintp-1
            endif
         endif

      enddo

      end
      subroutine orth(z,size,npred,intp,proj,uppm,m,n)
      integer size,npred,intp(npred+1),m,n
      double precision z(size,npred),proj(0:npred,npred+1),
     &     uppm(npred+1,npred+1)
      integer i,j,p

      do i=m,n
         proj(0,i)=1.0d0
         do p=1,npred
            proj(p,i)=z(intp(i),p)
         enddo
         do j=1,i-1
            uppm(j,i)=0.0d0
            do p=0,npred
               uppm(j,i)=uppm(j,i)+proj(p,j)*proj(p,i)
            enddo
            do p=0,npred
               proj(p,i)=proj(p,i)-uppm(j,i)*proj(p,j)
            enddo
         enddo
         uppm(i,i)=0.0d0
         do p=0,npred
            uppm(i,i)=uppm(i,i)+proj(p,i)**2
         enddo
         uppm(i,i)=dsqrt(uppm(i,i))
         do p=0,npred
            proj(p,i)=proj(p,i)/uppm(i,i)
         enddo
      enddo

      end
      subroutine line(rsdl,dlt,z,wt,size,npred,proj,zbar,
     &     intp,rfrac,nintp,lzr,beta,ocode,drct,lp)
      integer size,npred,dlt(size),intp(npred+1),nintp,lzr(size),
     &     ocode(2)
      double precision rsdl(size),z(size,npred),wt(size),
     &     proj(0:npred,npred+1),zbar(0:npred),rfrac(npred+1),
     &     beta(0:npred)
      double precision drct(0:npred),lp(size)
      logical first
      integer i,j,p,act
      double precision dtmp,dcmp

      ocode(1)=0
      ocode(2)=0

 10   continue
      do p=0,npred
         drct(p)=zbar(p)
      enddo
      do i=1,nintp
         dtmp=0.0d0
         do p=0,npred
            dtmp=dtmp+proj(p,i)*drct(p)
         enddo
         do p=0,npred
            drct(p)=drct(p)-dtmp*proj(p,i)
         enddo
      enddo

      dtmp=0.0d0
      do p=0,npred
         dtmp=dtmp+drct(p)**2
      enddo
      if(dtmp .lt. 1.0d-20) return

 20   continue
      do i=1,size
         if(iabs(lzr(i)) .eq. 1 .and. dlt(i) .eq. 1) then
            call caseone(dlt,z,wt,size,npred,proj,zbar,
     &           intp,rfrac,nintp,lzr,ocode,drct,lp,i,act)
            if(act .eq. 1) return
         endif
      enddo
c     check on censored cases of lzr=-1,1
      do i=1,size
         if(iabs(lzr(i)) .eq. 1 .and. dlt(i) .eq. 0) then
            call caseone(dlt,z,wt,size,npred,proj,zbar,
     &           intp,rfrac,nintp,lzr,ocode,drct,lp,i,act)
            if(act .eq. 1) return
            if(act .eq. 2) goto 10
         endif
      enddo

      first=.true.
      do i=1,size
         if(iabs(lzr(i)) .eq. 2) then
            lp(i)=drct(0)
            do p=1,npred
               lp(i)=lp(i)+z(i,p)*drct(p)
            enddo
            if(lp(i)*dble(lzr(i)) .gt. 2.0d-10) then
               dcmp=rsdl(i)/lp(i)
               if(first .or. dcmp .lt. dtmp) then
                  first=.false.
                  dtmp=dcmp
               endif
            endif
         endif
      enddo
      if(.not. first) then
         do i=1,size
            if(iabs(lzr(i)) .eq. 1 .and.
     &           (rsdl(i)-dtmp*lp(i))*dble(lzr(i))
     &           .le. -1.0d-10) return
         enddo
         do i=1,size
            if(lzr(i) .ne. 0) then
               rsdl(i)=rsdl(i)-dtmp*lp(i)
               j=1
               if(lzr(i) .lt. 0) j=-1
               if(dabs(rsdl(i)) .lt. 1.0d-10) then
                  lzr(i)=1
               else
                  lzr(i)=2
               endif
               lzr(i)=j*lzr(i)
            endif
         enddo
         do p=0,npred
            beta(p)=beta(p)+dtmp*drct(p)
         enddo
         ocode(1)=1

         goto 20
      endif

      end
      subroutine caseone(dlt,z,wt,size,npred,proj,zbar,
     &     intp,rfrac,nintp,lzr,ocode,drct,lp,i,act)
      integer size,npred,dlt(size),intp(npred+1),nintp,lzr(size),
     &     ocode(2),i,act
      double precision z(size,npred),wt(size),
     &     proj(0:npred,npred+1),zbar(0:npred),rfrac(npred+1)
      double precision drct(0:npred),lp(size)

      integer j,p
      double precision dtmp,dcmp

      act=0
      lp(i)=drct(0)
      do p=1,npred
         lp(i)=lp(i)+z(i,p)*drct(p)
      enddo
      dcmp=lp(i)*dble(lzr(i))
      if(dcmp .gt. 1.0d-10) then
         if(dlt(i) .eq. 1) then
            act=1
            nintp=nintp+1
            rfrac(nintp)=dble(lzr(i)+1)/2.0d0
            intp(nintp)=i
            lzr(i)=0
            ocode(2)=1
         else
            proj(0,npred+1)=1.0d0
            do p=1,npred
               proj(p,npred+1)=z(i,p)
            enddo
            do j=1,nintp
               dtmp=0.0d0
               do p=0,npred
                  dtmp=dtmp+proj(p,j)*proj(p,npred+1)
               enddo
               do p=0,npred
                  proj(p,npred+1)=proj(p,npred+1)
     &                 -dtmp*proj(p,j)
               enddo
            enddo
            dtmp=proj(0,npred+1)
            do p=1,npred
               dtmp=dtmp+proj(p,npred+1)*z(i,p)
            enddo
            dtmp=dtmp*wt(i)
            if(dcmp-dtmp .lt. -1.0d-10) then
               act=1
               nintp=nintp+1
               intp(nintp)=i
               dtmp=dcmp/dtmp
               zbar(0)=zbar(0)-dble(lzr(i))*wt(i)*dtmp
               do p=1,npred
                  zbar(p)=zbar(p)-dble(lzr(i))*z(i,p)
     &                 *wt(i)*dtmp
               enddo
               if(lzr(i) .eq. 1) then
                  rfrac(nintp)=1.0d0-dtmp
               else
                  rfrac(nintp)=dtmp
               endif
               lzr(i)=0
               ocode(2)=1
            else
               act=2
               zbar(0)=zbar(0)-wt(i)*dble(lzr(i))
               do p=1,npred
                  zbar(p)=zbar(p)-z(i,p)*wt(i)*dble(lzr(i))
               enddo
               lzr(i)=-lzr(i)
            endif
         endif
      endif

      end
      subroutine mono(para,npred,npt,
     &     initau,zch,nzch,
     &     monoindx,nmono,
     &     nq,taus,mu,
     &     diffhb)
      integer npred,npt,nzch,monoindx(npt),nmono,nq
      double precision para(0:npred+1,npt),
     &     initau,zch(nzch,0:npred),taus(nq),mu(0:npred,nq),
     &     diffhb(0:npred)

      logical mnt
      integer i,j,p,r,ini,mk
      double precision dtmp
      do i=1,npt
         monoindx(i)=0
      enddo
      ini=1
      do while(ini .le. npt .and.
     &     para(npred+1,ini) .lt. initau)
         ini=ini+1
      enddo
      if(ini .gt. npt) ini=ini-1
      monoindx(ini)=1

      mk=ini
      i=ini-1
      do while(i .ge. 1)
         do p=0,npred
            diffhb(p)=-para(p,i)+para(p,mk)
         enddo
         mnt=.true.
         r=1
         do while(mnt .and. r .le. nzch)
            dtmp=0.0d0
            do p=0,npred
               dtmp=dtmp+zch(r,p)*diffhb(p)
            enddo
            if(dtmp .lt. 0.0d0) mnt=.false.
            r=r+1
         enddo
         if(mnt) then
            mk=i
            monoindx(i)=1
         endif
         i=i-1
      enddo

c     forward
      mk=ini
      i=ini+1
      do while(i .le. npt)
         do p=0,npred
            diffhb(p)=para(p,i)-para(p,mk)
         enddo
         mnt=.true.
         r=1
         do while(mnt .and. r .le. nzch)
            dtmp=0.0d0
            do p=0,npred
               dtmp=dtmp+zch(r,p)*diffhb(p)
            enddo
            if(dtmp .lt. 0.0d0) mnt=.false.
            r=r+1
         enddo
         if(mnt) then
            mk=i
            monoindx(i)=1
         endif
         i=i+1
      enddo

      nmono=0
      do i=1,npt
         if(monoindx(i) .eq. 1) then
            nmono=nmono+1
            monoindx(nmono)=i
         endif
      enddo

      j=1
      do i=1,nmono
         do while(j .le. nq .and.
     &        (para(npred+1,monoindx(i)) .ge. taus(j) .or.
     &        i .eq. nmono))
            if(i .eq. 1 .or.
     &           para(npred+1,monoindx(i)) .lt. taus(j)) then
               do p=0,npred
                  mu(p,j)=para(p,monoindx(i))
               enddo
            else
               dtmp=(taus(j)-para(npred+1,monoindx(i-1)))/
     &              (para(npred+1,monoindx(i))
     &              -para(npred+1,monoindx(i-1)))
               do p=0,npred
                  mu(p,j)=(1.0d0-dtmp)*para(p,monoindx(i-1))
     &                 +dtmp*para(p,monoindx(i))
               enddo
            endif
            j=j+1
         enddo
      enddo

      end
