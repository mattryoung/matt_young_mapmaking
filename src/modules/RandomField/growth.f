      SUBROUTINE growth(z, omegam, omegal, w, D)
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     %     growth.F                                        %
c     %                                                     %
c     %     This subroutine calculates the growth factor    %
c     %     of linear perturbations at z, D(z), which is    %
c     %     normalized to D(0)=1.                           %
c     %                                                     %
c     %     References:                                     %
c     %        Carrol, Press & Turner (1992)                %
c     %        Ma, Caldwell, Bode & Wang (1999) for w/=-1   %
c     %        Silveira & Waga (1994) for exact solutions   %
c     %                                                     %
c     %                           1999. 6.28     E.Komatsu  %
c     %		(subroutine)	  2000.11.18     E.Komatsu  %
c     %         (quintessense)    2002. 1.22     E.Komatsu  %
c     %         (add selectors)   2002. 1.23     E.Komatsu  %
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      implicit none
c
c output
c
      real*8 D ! growth factor of linear perturbations at z
c
c inputs
c
      real*8 z                 ! redshift
      real*8 omegam, omegal, w ! cosmological parameters
c
      real*8 x, x2, x3, x3w
      real*8 omega, lambda     ! Omega_m(z), Omega_lambda(z)
      real*8 g, g0, t, t0
      complex hypgeo, a, b, c, y, y0 ! for Hypergeometric functions
c
      x  = 1d0+z
      x2 = x**2d0
      x3 = x**3d0
      x3w= x**(3d0*w)
c
c calculate Omega_m(z) and Omega_lambda(z)
c
      omega = omegam*x3
     &     /( omegam*x3 + (1d0-omegam-omegal)*x2 + omegal*x3*x3w )
      lambda= omegal*x3*x3w
     &     /( omegam*x3 + (1d0-omegam-omegal)*x2 + omegal*x3*x3w )
c
c Fitting formulae for a w=-1 universe [Carrol, Press & Turner (1992)]
c
      g = 2.5d0*omega
     &     /( omega**(4d0/7d0)-lambda
     &       + ( 1d0+omega/2d0 )*( 1d0+lambda/70d0 ) )
      g0= 2.5d0*omegam
     &     /( omegam**(4d0/7d0)-omegal
     &       + ( 1d0+omegam/2d0 )*( 1d0+omegal/70d0 ) )
c
c Fitting formulae for a w/=-1 universe [Ma, Caldwell, Bode & Wang (1999)] 
c                << for spatially flat universe only >>
c
      t = -( 0.255+0.305*w+0.0027/w )*( 1-omega )
     &    -( 0.366+0.266*w-0.07/w )*dlog(omega)
      t0= -( 0.255+0.305*w+0.0027/w )*( 1-omegam )
     &    -( 0.366+0.266*w-0.07/w )*dlog(omegam)
      g = g*(-w)**t
      g0= g0*(-w)**t0
c
c Calculate the growth factor
c
      D= (g/x)/g0
c
#ifdef GROWTH_EXACT
c
c Exact solution for a spatially flat universe [Silveira & Waga (1994)]
c
      a = -1d0/(3d0*w)
      b = (w-1d0)/(2d0*w)
      c = 1d0-5d0/(6d0*w)
      y = -x3w*(1d0-omegam)/omegam
      y0= -(1d0-omegam)/omegam
c
      D = ( hypgeo(a,b,c,y)/hypgeo(a,b,c,y0) )/x
#endif

      return
      END

c*************************************************************************
      FUNCTION hypgeo(a,b,c,z)
      COMPLEX hypgeo,a,b,c,z
      REAL EPS
      PARAMETER (EPS=1.e-6)
CU    USES bsstep,hypdrv,hypser,odeint
      INTEGER kmax,nbad,nok
      EXTERNAL bsstep,hypdrv
      COMPLEX z0,dz,aa,bb,cc,y(2)
      COMMON /hypg/ aa,bb,cc,z0,dz
      COMMON /path/ kmax
      kmax=0
      if (real(z)**2+aimag(z)**2.le.0.25) then
        call hypser(a,b,c,z,hypgeo,y(2))
        return
      else if (real(z).lt.0.) then
        z0=cmplx(-0.5,0.)
      else if (real(z).le.1.0) then
        z0=cmplx(0.5,0.)
      else
        z0=cmplx(0.,sign(0.5,aimag(z)))
      endif
      aa=a
      bb=b
      cc=c
      dz=z-z0
      call hypser(aa,bb,cc,z0,y(1),y(2))
      call odeint(y,4,0.,1.,EPS,.1,.0001,nok,nbad,hypdrv,bsstep)
      hypgeo=y(1)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software D041&0(9p#3.

c*************************************************************************
      SUBROUTINE hypser(a,b,c,z,series,deriv)
      INTEGER n
      COMPLEX a,b,c,z,series,deriv,aa,bb,cc,fac,temp
      deriv=cmplx(0.,0.)
      fac=cmplx(1.,0.)
      temp=fac
      aa=a
      bb=b
      cc=c
      do 11 n=1,1000
        fac=fac*aa*bb/cc
        deriv=deriv+fac
        fac=fac*z/n
        series=temp+fac
        if (series.eq.temp) return
        temp=series
        aa=aa+1.
        bb=bb+1.
        cc=cc+1.
11    continue
      pause 'convergence failure in hypser'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software D041&0(9p#3.

c*************************************************************************
      SUBROUTINE hypdrv(s,y,dyds)
      REAL s
      COMPLEX y(2),dyds(2),aa,bb,cc,z0,dz,z
      COMMON /hypg/ aa,bb,cc,z0,dz
      z=z0+s*dz
      dyds(1)=y(2)*dz
      dyds(2)=(aa*bb*y(1)-(cc-(aa+bb+1.)*z)*y(2))*dz/(z*(1.-z))
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software D041&0(9p#3.

c*************************************************************************
      SUBROUTINE bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER nv,NMAX,KMAXX,IMAX
      REAL eps,hdid,hnext,htry,x,dydx(nv),y(nv),yscal(nv),SAFE1,SAFE2,
     *REDMAX,REDMIN,TINY,SCALMX
      PARAMETER (NMAX=50,KMAXX=8,IMAX=KMAXX+1,SAFE1=.25,SAFE2=.7,
     *REDMAX=1.e-5,REDMIN=.7,TINY=1.e-30,SCALMX=.1)
CU    USES derivs,mmid,pzextr
      INTEGER i,iq,k,kk,km,kmax,kopt,nseq(IMAX)
      REAL eps1,epsold,errmax,fact,h,red,scale,work,wrkmin,xest,xnew,
     *a(IMAX),alf(KMAXX,KMAXX),err(KMAXX),yerr(NMAX),ysav(NMAX),
     *yseq(NMAX)
      LOGICAL first,reduct
      SAVE a,alf,epsold,first,kmax,kopt,nseq,xnew
      EXTERNAL derivs
      DATA first/.true./,epsold/-1./
      DATA nseq /2,4,6,8,10,12,14,16,18/
      if(eps.ne.epsold)then
        hnext=-1.e29
        xnew=-1.e29
        eps1=SAFE1*eps
        a(1)=nseq(1)+1
        do 11 k=1,KMAXX
          a(k+1)=a(k)+nseq(k+1)
11      continue
        do 13 iq=2,KMAXX
          do 12 k=1,iq-1
            alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.)*(2*k+
     *1)))
12        continue
13      continue
        epsold=eps
        do 14 kopt=2,KMAXX-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1
14      continue
1       kmax=kopt
      endif
      h=htry
      do 15 i=1,nv
        ysav(i)=y(i)
15    continue
      if(h.ne.hnext.or.x.ne.xnew)then
        first=.true.
        kopt=kmax
      endif
      reduct=.false.
2     do 17 k=1,kmax
        xnew=x+h
        if(xnew.eq.x)pause 'step size underflow in bsstep'
        call mmid(ysav,dydx,nv,x,h,nseq(k),yseq,derivs)
        xest=(h/nseq(k))**2
        call pzextr(k,xest,yseq,y,yerr,nv)
        if(k.ne.1)then
          errmax=TINY
          do 16 i=1,nv
            errmax=max(errmax,abs(yerr(i)/yscal(i)))
16        continue
          errmax=errmax/eps
          km=k-1
          err(km)=(errmax/SAFE1)**(1./(2*km+1))
        endif
        if(k.ne.1.and.(k.ge.kopt-1.or.first))then
          if(errmax.lt.1.)goto 4
          if(k.eq.kmax.or.k.eq.kopt+1)then
            red=SAFE2/err(km)
            goto 3
          else if(k.eq.kopt)then
            if(alf(kopt-1,kopt).lt.err(km))then
              red=1./err(km)
              goto 3
            endif
          else if(kopt.eq.kmax)then
            if(alf(km,kmax-1).lt.err(km))then
              red=alf(km,kmax-1)*SAFE2/err(km)
              goto 3
            endif
          else if(alf(km,kopt).lt.err(km))then
            red=alf(km,kopt-1)/err(km)
            goto 3
          endif
        endif
17    continue
3     red=min(red,REDMIN)
      red=max(red,REDMAX)
      h=h*red
      reduct=.true.
      goto 2
4     x=xnew
      hdid=h
      first=.false.
      wrkmin=1.e35
      do 18 kk=1,km
        fact=max(err(kk),SCALMX)
        work=fact*a(kk+1)
        if(work.lt.wrkmin)then
          scale=fact
          wrkmin=work
          kopt=kk+1
        endif
18    continue
      hnext=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
        fact=max(scale/alf(kopt-1,kopt),SCALMX)
        if(a(kopt+1)*fact.le.wrkmin)then
          hnext=h/fact
          kopt=kopt+1
        endif
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software D041&0(9p#3.

c*************************************************************************
      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,
     *rkqs)
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      REAL eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs,rkqs
      PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=200,TINY=1.e-30)
      INTEGER i,kmax,kount,nstp
      REAL dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),
     *yp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.*dxsav
      do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do 12 i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
12      continue
        if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          return
        endif
        if(abs(hnext).lt.hmin) pause
     *'stepsize smaller than minimum in odeint'
        h=hnext
16    continue
      pause 'too many steps in odeint'
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software D041&0(9p#3.

c*************************************************************************
      SUBROUTINE mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs)
      INTEGER nstep,nvar,NMAX
      REAL htot,xs,dydx(nvar),y(nvar),yout(nvar)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
      INTEGER i,n
      REAL h,h2,swap,x,ym(NMAX),yn(NMAX)
      h=htot/nstep
      do 11 i=1,nvar
        ym(i)=y(i)
        yn(i)=y(i)+h*dydx(i)
11    continue
      x=xs+h
      call derivs(x,yn,yout)
      h2=2.*h
      do 13 n=2,nstep
        do 12 i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
12      continue
        x=x+h
        call derivs(x,yn,yout)
13    continue
      do 14 i=1,nvar
        yout(i)=0.5*(ym(i)+yn(i)+h*yout(i))
14    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software D041&0(9p#3.

c*************************************************************************
      SUBROUTINE pzextr(iest,xest,yest,yz,dy,nv)
      INTEGER iest,nv,IMAX,NMAX
      REAL xest,dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=13,NMAX=50)
      INTEGER j,k1
      REAL delta,f1,f2,q,d(NMAX),qcol(NMAX,IMAX),x(IMAX)
      SAVE qcol,x
      x(iest)=xest
      do 11 j=1,nv
        dy(j)=yest(j)
        yz(j)=yest(j)
11    continue
      if(iest.eq.1) then
        do 12 j=1,nv
          qcol(j,1)=yest(j)
12      continue
      else
        do 13 j=1,nv
          d(j)=yest(j)
13      continue
        do 15 k1=1,iest-1
          delta=1./(x(iest-k1)-xest)
          f1=xest*delta
          f2=x(iest-k1)*delta
          do 14 j=1,nv
            q=qcol(j,k1)
            qcol(j,k1)=dy(j)
            delta=d(j)-q
            dy(j)=f1*delta
            d(j)=f2*delta
            yz(j)=yz(j)+dy(j)
14        continue
15      continue
        do 16 j=1,nv
          qcol(j,iest)=dy(j)
16      continue
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software D041&0(9p#3.
