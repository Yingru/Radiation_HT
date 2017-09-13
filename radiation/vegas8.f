c    This is the subroutine to integrate a give function by Monto Carlo method.

      SUBROUTINE vegas8(region,ndim,fxn,init,ncall,itmx,nprn,
     *	tgral,sd,chi2a)
  
      INTEGER init,itmx,ncall,ndim,nprn,NDMX,MXDIM
      DOUBLE PRECISION tgral,chi2a,sd,region(2*ndim),fxn,ALPH,TINY
      PARAMETER (ALPH=1.5,NDMX=50,MXDIM=10,TINY=1.e-30)
      EXTERNAL fxn
CU    	USES fxn,ran28,rebin8
      INTEGER i,idum,it,j,k,mds,nd,ndo,ng,npg,ia(MXDIM),kg(MXDIM)
      DOUBLE PRECISION calls,dv2g,dxg,f,f2,f2b,fb,rc,ti,tsi,wgt,
     *  xjac,xn,xnd,xo,
     *	d(NDMX,MXDIM),di(NDMX,MXDIM),dt(MXDIM),dx(MXDIM),r(NDMX),
     *	x(MXDIM),xi(NDMX,MXDIM),xin(NDMX),ran28,trial,pp(MXDIM),qq
      DOUBLE PRECISION schi,si,swgt,temp
      COMMON /ranno/ idum
      SAVE

c	write (*,*) "1:", x(1), x(2), x(3)

      if(init.le.0)then
        mds=1
        ndo=1
        do 11 j=1,ndim
          xi(1,j)=1.
11      continue
      endif

c	write (*,*) "2:", x(1), x(2), x(3)

      if (init.le.1)then
        si=0.
        swgt=0.
        schi=0.
      endif
      if (init.le.2)then
        nd=NDMX
        ng=1
        if(mds.ne.0)then
          ng=(ncall/2.+0.25)**(1./ndim)
          mds=1
          if((2*ng-NDMX).ge.0)then
            mds=-1
            npg=ng/NDMX+1
            nd=ng/npg
            ng=npg*nd
          endif
        endif
        k=ng**ndim
        npg=max(ncall/k,2)
        calls=npg*k
        dxg=1./ng
        dv2g=(calls*dxg**ndim)**2/npg/npg/(npg-1.)
        xnd=nd
        dxg=dxg*xnd
        xjac=1./calls
        do 12 j=1,ndim
          dx(j)=region(j+ndim)-region(j)
          xjac=xjac*dx(j)
12      continue
        if(nd.ne.ndo)then
          do 13 i=1,nd
            r(i)=1.
13        continue
          do 14 j=1,ndim
		  temp=ndo/xnd
            call rebin8(temp,nd,r,xin,xi(1,j))
c            call rebin8(ndo/xnd,nd,r,xin,xi(1,j))

14        continue
          ndo=nd
        endif
     
      endif
 
c	write (*,*) "3:", x(1), x(2), x(3)

      do 28 it=1,itmx

        ti=0.
        tsi=0.
        do 16 j=1,ndim
          kg(j)=1
          do 15 i=1,nd
            d(i,j)=0.
            di(i,j)=0.
15        continue
16      continue
10      continue
        fb=0.
        f2b=0.
        do 19 k=1,npg
          wgt=xjac
          do 17 j=1,ndim
            xn=(kg(j)-ran28(idum))*dxg+1.
            ia(j)=max(min(int(xn),NDMX),1)
            if(ia(j).gt.1)then
              xo=xi(ia(j),j)-xi(ia(j)-1,j)
              rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
            else
              xo=xi(ia(j),j)
              rc=(xn-ia(j))*xo
            endif
          x(j)=region(j)+rc*dx(j)
          wgt=wgt*xo*xnd

17        continue


          f=wgt*fxn(x,wgt)
c          write(6,*) f/wgt,wgt,it
c          write (*,*) "4:", x(1), x(2), x(3), f
	if (f .gt. 0.d0 .OR. f.le. 0.d0) then
	else
	stop
	end if

          f2=f*f
          fb=fb+f
          f2b=f2b+f2
          do 18 j=1,ndim
            di(ia(j),j)=di(ia(j),j)+f
            if(mds.ge.0) d(ia(j),j)=d(ia(j),j)+f2
18        continue
19      continue
        f2b=sqrt(f2b*npg)
        f2b=(f2b-fb)*(f2b+fb)
        if (f2b.le.0.) f2b=TINY
          ti=ti+fb
          tsi=tsi+f2b
          if(mds.lt.0)then
            do 21 j=1,ndim
              d(ia(j),j)=d(ia(j),j)+f2b
21          continue
          endif
          do 22 k=ndim,1,-1
            kg(k)=mod(kg(k),ng)+1
            if(kg(k).ne.1) goto 10
22        continue
          tsi=tsi*dv2g
          wgt=1./tsi
          si=si+dble(wgt)*dble(ti)
          schi=schi+dble(wgt)*dble(ti)**2
          swgt=swgt+dble(wgt)
          tgral=si/swgt
          chi2a=max((schi-si*tgral)/(it-.99d0),0.d0)
          sd=sqrt(1./swgt)
          tsi=sqrt(tsi)

          do 25 j=1,ndim
            xo=d(1,j)
            xn=d(2,j)
            d(1,j)=(xo+xn)/2.
            dt(j)=d(1,j)
            do 24 i=2,nd-1
              rc=xo+xn
              xo=xn
              xn=d(i+1,j)
              d(i,j)=(rc+xn)/3.
              dt(j)=dt(j)+d(i,j)

24          continue
            d(nd,j)=(xo+xn)/2.
            dt(j)=dt(j)+d(nd,j)
25        continue


          do 27 j=1,ndim
            rc=0.
            do 26 i=1,nd
              if(d(i,j).lt.TINY) d(i,j)=TINY
                r(i)=((1.-d(i,j)/dt(j))/(log(dt(j))-log(d(i,j))))**ALPH
                rc=rc+r(i)
26        continue
		temp=rc/xnd
          call rebin8(temp,nd,r,xin,xi(1,j))
c	    call rebin8(rc/xnd,nd,r,xin,xi(1,j))
27      continue


28    continue
      return

      END
C  (C)	Copr. 1986-92 Numerical Recipes Software Dt+;39.


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      FUNCTION ran28(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      DOUBLE PRECISION ran28,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *	IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *	NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
      idum=max(-idum,1)
      idum2=idum
      do 11 j=NTAB+8,1,-1
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      if (j.le.NTAB) iv(j)=idum
11    continue
      iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran28=min(AM*iy,RNMX)
      return
      END
C  (C)	Copr. 1986-92 Numerical Recipes Software Dt+;39.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE rebin8(rc,nd,r,xin,xi)
      INTEGER nd
      DOUBLE PRECISION rc,r(*),xi(*),xin(*)
      INTEGER i,k
      DOUBLE PRECISION dr,xn,xo
      k=0
      xn=0.
      dr=0.
      do 11 i=1,nd-1
1       if(rc.gt.dr)then
          k=k+1
          dr=dr+r(k)
          xo=xn
          xn=xi(k)
          goto 1
        endif
        dr=dr-rc
        xin(i)=xn-(xn-xo)*dr/r(k)
11    continue
      do 12 i=1,nd-1
        xi(i)=xin(i)
12    continue
      xi(nd)=1.
      return
      END
C  (C) 	Copr. 1986-92 Numerical Recipes Software Dt+;39.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
