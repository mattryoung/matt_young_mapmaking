MODULE psubs_Dlinear
END MODULE psubs_Dlinear

!c -------------- THIS IS THE DLINEAR PACKAGE ------------

subroutine read_cosmology(Omt,h,Omvac,iamcurved,dcurv)
  USE intreal_types

  Omnr=Omx+OmB
  Omt=Omnr+Omvac
  Omcur=1-Omt
  if(Omcurv.eq.0.0) then
     iamcurved=0
     dcurv=1.e10
  elseif(Omcurv.gt.0.0) then
     iamcurved=-1
     dcurv=3000.0/sqrt(Omcurv)/h
  elseif(Omcurv.lt.0.0) then
     iamcurved=1
     dcurv=3000.0/sqrt(abs(Omcurv))/h
  endif
  return

end subroutine read_cosmology

subroutine get_cosmology
  USE intreal_types
  use univer
  use univer2
  use fac
  use univercurv
  use EdeS
  use h_var

  write(*,*) 'Omx,OmB,Omhdm,fhdmclus,Omvac,Omcurv,h '
  read(999,*)  Omx,OmB,Omhdm,fhdmclus,Omvac,Omcurv,h
  Omnr=OmB+OmX+Omhdm
  Omt=(Omnr+Omvac)

  if(Omcurv.eq.0.0) then
     iamcurved=0
     dcurv=1.e10
  elseif(Omcurv.gt.0.0) then
     iamcurved=-1
     dcurv=3000.0/sqrt(Omcurv)/h
  elseif(Omcurv.lt.0.0) then
     iamcurved=1
     dcurv=3000.0/sqrt(abs(Omcurv))/h
  endif
  iEdeS=0
  if(Omnr.eq.1.and.Omt.eq.1.and.Omhdm.eq.0) iEdeS=1
  Omer=4.1e-05/h/h
  Omer=0.0
  Omnr=Omnr-Omer
  call Dlinear_setup
  return
end subroutine get_cosmology

subroutine Dlinear_cosmology(om0,h0,omvac0,iamcurved0,dvurv0)
  USE intreal_types
  use univer
  use univer2
  use fac
  use univercurv
  use h_var
  use EdeS

  !c     Note: Dlinear_setup should make tables unless
  !c     we have a a=t^(2/3) universe. This is the only universe
  !c     for which all the t(a), a(t) conversions are analytic
  !c     i.e. iEdeS eq 1 iff we are flat, zero Lambda, zero relativistic
  !c     Note: omt refers to non-curvature  omega
  
  omnr0     = om0
  omx0      = om0
  omb0      = 0
  omcurv0   = 1 - om0 - omvac0
  omhdm0    = 0
  fhdmclus0 = 0
  omtvac0   = om0 + omvac0
  
  OmB=omb0
  Omx=omx0
  Omvac=omvac0
  Omcurv=omcurv0
  fhdmclus=fhdmclus0
  h=h0
  Omhdm=Omhdm0

  Omnr=OmB+OmX+Omhdm
  Omt=(Omnr+Omvac)
  if (abs(Omt-1.0).lt.0.0001) Omt=1.0
  if (Omtvac0.ne.Omt) then 
     write(6,*) 'Omt''s different ',Omt,Omtvac0
     stop
  endif
  if(Omnr0.ne.Omnr) then
     write(*,*) ' Omnr0.ne.Omnr '
     Omnr0=Omnr
  endif
  if(Omcurv.eq.0.0) then
     iamcurved=0
     dcurv=1.e10
  elseif(Omcurv.gt.0.0) then
     iamcurved=-1
     dcurv=3000.0/sqrt(Omcurv)/h
  elseif(Omcurv.lt.0.0) then
     iamcurved=1
     dcurv=3000.0/sqrt(abs(Omcurv))/h
  endif
  if(iamcurved0.ne.iamcurved) then
     iamcurved0=iamcurved
     dcurv0=dcurv
  endif
  !c      Omer=0.0
  iEdeS=0
  if(Omnr.eq.1.and.Omt.eq.1.and.Omhdm.eq.0) iEdeS=1
  Omer=4.1e-05/h/h
  Omer=0.0
  Omnr=Omnr-Omer
  call Dlinear_setup

  return
end subroutine Dlinear_cosmology

function chifn_Sigma(Sigma,iamcurved,dcurv)
  USE intreal_types
  if(iamcurved.eq.0) then
     chifn_Sigma=Sigma
  elseif(iamcurved.eq.-1) then
     xx=Sigma/dcurv
     chifn_Sigma=dcurv*alog(xx+sqrt(xx*xx+1.0))
  elseif(iamcurved.eq.1) then
     xx=Sigma/dcurv
     chifn_Sigma=dcurv*asin(xx)
  endif
  return
end function chifn_Sigma

function Sigmafn_chi(chi,iamcurved,dcurv)
  USE intreal_types
  if(iamcurved.eq.0) then
     Sigmafn_chi=chi
  elseif(iamcurved.eq.-1) then
     xx=chi/dcurv
     exx=exp(xx)
     Sigmafn_chi=dcurv*(exx-1.0/exx)*0.5
  elseif(iamcurved.eq.1) then
     xx=chi/dcurv
     Sigmafn_chi=dcurv*sin(xx)
  endif
  return
end function Sigmafn_chi

function afn(chi)
  USE intreal_types
  use array_bounds
  use Dlin_tab
  use Dlin_params
  use univer
  use univer2
  use h_var
  use EdeS
  !C  **   THIS GIVES THE EXPANSION FACTOR AS A FUNCTION OF COMOVING DISTANCE FROM THE ORIGIN
  parameter (chin0_1=1.0/6000.0)
  dimension chixa(4),ptab(4)
  if(imaketab_Dlin.eq.0) then
     !c         if((Omt.eq.1.0.and.Omnr.eq.1.0)) then
     if(iEdeS.eq.1) then
        afn=(1.0-chi*h*chin0_1)**2
        return
     endif
     if(Omt.lt.1.0) then
        write(*,*) ' THE ANALYTIC OPEN NEEDS WORK - MAKE TABLE '
        stop
     endif
     if(Omt.gt.1.0) then
        write(*,*) ' THE ANALYTIC CLOSED NEEDS WORK - MAKE TABLE '
        stop
     endif
  endif

  if(chi.gt.chimintab.and.chi.lt.chimaxtab) then
     call HUNT(chitab,ntab_Dlin,chi,kk)
     !c         mm=2
     if(kk.lt.2.or.kk.gt.ntab_Dlin-2) then
        nlint=2
        kkl=kk-1
     else 
        nlint=4
        kkl=kk-2
     endif
     do ikk=1,nlint
        chixa(ikk)=chitab(kkl+ikk)
        ptab(ikk)=atab(kkl+ikk)
     enddo
     call lagrint_one(chixa,ptab,pint,nlint,chi)
     afn=pint
     return
  elseif(chi.ge.chimaxtab) then
     afn= ((chin0-chi)/chininfty)**2
     return
  elseif(chi.le.chimintab) then
     WRITE(*,*) 'chi OFF END OF chiTABLE ',chi,chimintab
     afn=amaxtab
     return
  endif
  return
end function afn

real function chifn(a)
  USE intreal_types
  use array_bounds
  use Dlin_tab
  use Dlin_params
  use univer
  use univer2
  use h_var
  use EdeS
  !C **   THIS GIVES COMOVING DISTANCE FROM THE ORIGIN AS A FUNCTION OF THE EXPANSION FACTOR 
  parameter (chin0_1=1.0/6000.0)
  dimension axa(4),ptab(4)
  if(imaketab_Dlin.eq.0) then
     !c         if((Omt.eq.1.0.and.Omnr.eq.1.0)) then
     if(iEdeS.eq.1) then
        chifn=6000.0/h*(1.0-sqrt(a))
        return
     endif
     if(Omt.lt.1.0) then
        write(*,*) ' THE ANALYTIC OPEN NEEDS WORK - MAKE TABLE '
        stop
     endif
     if(Omt.gt.1.0) then
        write(*,*) ' THE ANALYTIC CLOSED NEEDS WORK - MAKE TABLE '
        stop
     endif
  endif

  if(a.gt.amintab.and.a.lt.amaxtab) then
     call HUNT(atab,ntab_Dlin,a,kk)
     !c         mm=2
     if(kk.lt.2.or.kk.gt.ntab_Dlin-2) then
        nlint=2
        kkl=kk-1
     else 
        nlint=4
        kkl=kk-2
     endif
     do ikk=1,nlint
        axa(ikk)=atab(kkl+ikk)
        ptab(ikk)=chitab(kkl+ikk)
     enddo
     call lagrint_one(axa,ptab,pint,nlint,a)
     chifn=pint
     return
  elseif(a.ge.amaxtab) then
     WRITE(*,*) ' a OFF END OF chiTABLE ',a,amaxtab
     chifn=0.0
     return
  elseif(a.le.amintab) then
     chifn=chin0-chininfty*sqrt(a)
     return
  endif
  return
end function chifn

real function afnofD(D)
  USE intreal_types
  use array_bounds
  use Dlin_tab
  use Dlin_params
  use univer
  use univer2
  use h_var
  use EdeS

  !C  **   THIS GIVES THE EXPANSION FACTOR AS A FUNCTION OF LINEAR GROWTH FACTOR
  dimension Dxa(4),ptab(4)
  if(imaketab_Dlin.eq.0) then
     !c         if((Omt.eq.1.0.and.Omnr.eq.1.0)) then
     if(iEdeS.eq.1) then
        afnofD=D
        return
     endif
     if(Omt.lt.1.0) then
        write(*,*) ' THE ANALYTIC OPEN NEEDS WORK - MAKE TABLE '
        stop
     endif
     if(Omt.gt.1.0) then
        write(*,*) ' THE ANALYTIC CLOSED NEEDS WORK - MAKE TABLE '
        stop
     endif
  endif

  if(D.gt.Dmintab.and.D.lt.Dmaxtab) then
     call HUNT(Dtab,ntab_Dlin,D,kk)
     !c         mm=2
     if(kk.lt.2.or.kk.gt.ntab_Dlin-2) then
        nlint=2
        kkl=kk-1
     else 
        nlint=4
        kkl=kk-2
     endif
     do ikk=1,nlint
        Dxa(ikk)=Dtab(kkl+ikk)
        ptab(ikk)=atab(kkl+ikk)
     enddo
     call lagrint_one(Dxa,ptab,pint,nlint,D)
     afnofD=pint
     return
  elseif(D.ge.Dmaxtab) then
     !c         WRITE(*,*) 'D OFF END OF chiTABLE ',D,Dmaxtab
     afnofD=D/D_atab(ntab_Dlin)
     return
  elseif(D.le.Dmintab) then
     afnofD=D/D_atab(1)
     return
  endif
  return
end function afnofD

real function Dfnofa(a)
  USE intreal_types
  use array_bounds
  use Dlin_tab
  use Dlin_params
  use univer
  use univer2
  use h_var
  use EdeS

  !C **   THIS GIVES LINEAR GROWTH FACTOR AS A FUNCTION OF THE EXPANSION FACTOR 

  dimension axa(4),ptab(4)
  if(imaketab_Dlin.eq.0) then
     !c         if((Omt.eq.1.0.and.Omnr.eq.1.0)) then
     if(iEdeS.eq.1) then
        Dfnofa=a
        return
     endif
     if(Omt.lt.1.0) then
        write(*,*) ' THE ANALYTIC OPEN NEEDS WORK - MAKE TABLE '
        stop
     endif
     if(Omt.gt.1.0) then
        write(*,*) ' THE ANALYTIC CLOSED NEEDS WORK - MAKE TABLE '
        stop
     endif
  endif

  if(a.gt.amintab.and.a.lt.amaxtab) then
     call HUNT(atab,ntab_Dlin,a,kk)
     !c         mm=2
     if(kk.lt.2.or.kk.gt.ntab_Dlin-2) then
        nlint=2
        kkl=kk-1
     else 
        nlint=4
        kkl=kk-2
     endif
     do ikk=1,nlint
        axa(ikk)=atab(kkl+ikk)
        !c            ptab(ikk)=Dtab(kkl+ikk)
        ptab(ikk)=D_atab(kkl+ikk)
     enddo
     call lagrint_one(axa,ptab,pint,nlint,a)
     !c         Dfnofa=pint
     Dfnofa=a*pint
     return
  elseif(a.ge.amaxtab) then
     Dfnofa=a*D_atab(ntab_Dlin)
     return
  elseif(a.le.amintab) then
     Dfnofa=a*D_atab(1)
     return
  endif
  return
end function Dfnofa

real function geta(t)
  USE intreal_types
  use array_bounds
  use Dlin_tab
  use Dlin_params
  use univer
  use univer2
  use h_var
  use EdeS

  !c You must send t in 'dlinear' units
  !c To this end, call gett(a) and calculate a conversion for your
  !c t to t_Dlin units and perform it for all calls.
  !c gett(a) is similar

  dimension axa(4),txt(4)
  if(imaketab_Dlin.eq.0) then
     if(iEdeS.eq.1) then
        !c Note we assume iEdeS is 1 only if flat, omnr=omt=1, i.e. a=t**(2/3)
        !c            if (omnr.lt.1.0) then
        !c               tstar=sqrt(2*omnr/(1.0-omnr))/3.0*(1+z0)**1.5
        !c               geta=((omnr/(1.0-omnr))**0.3333333*
        !c     $       (log(t1/tstar+sqrt(1+(t1/tstar)**2)))**0.666666666667)
        !c            else
        t1=1./h
        geta=(t*1.5/t1)**(2./3.)
        return
     endif
     write(6,*) 'No Tables installed for this cosmology'
     stop
  endif
  if(t.gt.tmintab.and.t.lt.tmaxtab) then
     call HUNT(ttab,ntab_Dlin,t,kk)
     if(kk.lt.2.or.kk.gt.ntab_Dlin-2) then
        nlint=2
        kkl=kk-1
     else 
        nlint=4
        kkl=kk-2
     endif
     do  ikk=1,nlint
        axa(ikk)=atab(kkl+ikk)
        txt(ikk)=ttab(kkl+ikk)
     enddo
     call lagrint_one(txt,axa,a,nlint,t)
     geta=a
     return
  elseif(t.ge.tmaxtab) then
     WRITE(*,*) ' OFF END OF tTABLE ',t,tmintab,tmaxtab
     stop
  elseif(t.le.tmintab) then
     geta=amintab*(t/tmintab)**(2./3.)
     return
  endif
end function geta

real function gett(a)
  USE intreal_types
  use array_bounds
  use Dlin_tab
  use Dlin_params
  use univer
  use univer2
  use h_var
  use EdeS

  !c You must send t in 'dlinear' units
  !c To this end, call gett(a) and calculate a conversion for your
  !c t to t_Dlin units and perform it for all calls.
  !c gett(a) is similar

  dimension axa(4),txt(4)
  if(imaketab_Dlin.eq.0) then
     if(iEdeS.eq.1) then
        !c Note we assume iEdeS is 1 only if flat, omnr=omt=1, i.e. a=t**(2/3)
        !c Non-zero Lambda:
        !c            if (omnr.lt.1.0) then
        !c               tstar=sqrt(2*omnr/(1.0-omnr))/3.0*(1+z0)**1.5
        !c               geta=((omnr/(1.0-omnr))**0.3333333*
        !c     $       (log(t1/tstar+sqrt(1+(t1/tstar)**2)))**0.666666666667)
        !c            else
        t1=1./(1.5*h)
        gett=t1*a**1.5
        return
     endif
     write(6,*) 'No Tables installed for this cosmology'
     stop
  endif
  if(a.gt.amintab.and.a.lt.amaxtab) then
     call HUNT(atab,ntab_Dlin,a,kk)
     if(kk.lt.2.or.kk.gt.ntab_Dlin-2) then
        nlint=2
        kkl=kk-1
     else 
        nlint=4
        kkl=kk-2
     endif
     do  ikk=1,nlint
        axa(ikk)=atab(kkl+ikk)
        txt(ikk)=ttab(kkl+ikk)
     enddo
     call lagrint_one(axa,txt,t,nlint,a)
     gett=t
     return
  elseif(a.ge.amaxtab) then
     WRITE(*,*) ' OFF END OF aTABLE ',a,amintab,amaxtab
     stop
  elseif(a.le.amintab) then
     gett=tmintab*(a/amintab)**1.5
     return
  endif
end function gett

real function Dlinear(a,chi,HD_Ha,D_a)
  USE intreal_types
  use array_bounds
  use Dlin_tab
  use Dlin_params
  use univer
  use univer2
  use h_var
  use EdeS

  !C **   THIS GIVES THE LINEAR GROWTH FACTOR FOR OPEN UNIVERSES IF THE U IS CLOSED OR IF THERE IS VACUUM ENERGY etc

  dimension axa(4),ptab(4,3),pint(3)
  if(imaketab_Dlin.eq.0) then
     if(iEdeS.eq.1) then
        D_a=1.0
        Dlinear=D_a*a
        HD_Ha=1.0
        chi=6000.0/h*(1.0-sqrt(a))
        return
     endif
     if(Omt.lt.1.0) then
        write(*,*) ' THE ANALYTIC OPEN NEEDS WORK - MAKE TABLE '
        stop
        aomt=1.0/(1.0/omt-1.0)
        y0=1.0/aomt
        y=a/aomt
        uu0=sqrt(1.0/y0+1.0)*alog(sqrt(1.0+y0)-sqrt(y0))
        y0=1.0+3.0*(1.0+uu0)/y0
        dainfty=0.4/y0/aomt
        if(y.lt.0.01) then
           D_a=dainfty
           Dlinear=D_a*a
           HD_Ha=1.0
           return
        else
           uu=sqrt(1.0/y+1.0)*alog(sqrt(1.0+y)-sqrt(y))
           Dlinear=1.0+3.0*(1.0+uu)/y
           Dlinear=Dlinear/y0
           D_a=Dlinear/a
           HD_Ha=(1.5*uu/(1.0+y)-4.5/y*(1.0+uu))/y*2.5*dainfty/D_a
           return
        endif
     endif
     if(Omt.gt.1.0) then
        write(*,*) ' THE ANALYTIC CLOSED NEEDS WORK - MAKE TABLE '
        stop
        aomt=1.0/abs((1.0/omt-1.0))
        y0=1.0/aomt
        y=a/aomt
        uu0=sqrt(1.0/y0-1.0)*atan(sqrt(y0/(1.0-y0)))
        y0=-1.0+3.0*(1.0-uu0)/y0
        dainfty=0.4/y0/aomt
        if(y.lt.0.01) then
           D_a=dainfty
           Dlinear=D_a*a
           HD_Ha=1.0
           return
        else
           uu=sqrt(1.0/y-1.0)*atan(sqrt(y/(1.0-y)))
           Dlinear=-1.0+3.0*(1.0-uu)/y
           Dlinear=Dlinear/y0
           D_a=Dlinear/a
           HD_Ha=(1.5*uu/(1.0-y)-4.5/y*(1.0-uu))/y*2.5*dainfty/D_a
           Dlinear=D_a*a
           return
        endif
     endif
  endif

  if(a.gt.amintab.and.a.lt.amaxtab) then
     call HUNT(atab,ntab_Dlin,a,kk)
     if(kk.lt.2.or.kk.gt.ntab_Dlin-2) then
        nlint=2
        kkl=kk-1
     else 
        nlint=4
        kkl=kk-2
     endif
     do ikk=1,nlint
        axa(ikk)=atab(kkl+ikk)
        ptab(ikk,1)=D_atab(kkl+ikk)
        ptab(ikk,2)=dlnD_dlnatab(kkl+ikk)
        ptab(ikk,3)=chitab(kkl+ikk)
     enddo
     call lagrint_Dlin(axa,ptab,pint,nlint,a)
     D_a=pint(1)
     HD_Ha=pint(2)
     chi=pint(3)
     Dlinear=D_a*a
     return
  elseif(a.ge.amaxtab) then
     D_a=D_a0
     HD_Ha=HD_Ha0
     Dlinear=D_a*a
     chi=0.0
     return
  elseif(a.le.amintab) then
     D_a=D_ainfty
     HD_Ha=HD_Hainfty
     Dlinear=D_a*a
     chi=chin0-chininfty*sqrt(a)
     return
  endif
  return
end function Dlinear

subroutine Dlinear_setup
  USE intreal_types
  use array_bounds
  use Dlin_tab
  use Dlin_params
  use count
  use fac
  use qtime
  use univer
  use univer2
  use decay
  use tau_vars
  use h_var
  use EdeS

  parameter(H0_hinv=3.0856775806e17,H0_hinvyr=H0_hinv/3.155760e07)
  dimension y(50)

  !c     Note: Dlinear_setup should be hardwired to make tables unless
  !c     we have a a=t^(2/3) universe. This is the only universe
  !c     for which all the t(a), a(t) conversions are analytic
  !c     i.e. iEdeS eq 1 iff we are flat, zero Lambda, zero relativistic
  if(iEdeS.eq.1) then
     D_ainfty=1.0
     D_a0=1.0
     HD_Hainfty=1.0
     HD_Ha0=1.0
     chininfty=6000.0/h
     chin0=6000.0/h 
     imaketab_Dlin=0
     return
  endif
  imaketab_Dlin=1
  amintab=1001
  amaxtab=0.5
  amintab=1.0/amintab
  amaxtab=1.0/amaxtab
  fnr=1.0
  fer=0.0
  ldecay=0
  Omnrdi=0.0
  taud=1.0
  Omnrd=Omnrdi
  nmax=itabmax-1
  now=1
  ainit=1.0e-07
  alogn=-16.11809
  a = ainit
  amlog = log(amaxtab)
  hh=(amlog-alogn)/nmax
  !c   ** integration var is ln(a), y1=time in 10**10y units, 
  !c     y2=fb, y3=dlnb/dlna
  !c   **   y4=omerd*a, y5=tau, ftau*2/H0/sqrt(Omnr)*a**1.5=tau 
  y(1) = tin_Dlin(a,tauin_Dlin)

  y(2) = log(a)
  !c y(3)=0 is the initial value of dlnD/dlna for radiation dominated (ainit<<aeq)
  !c We set it to 1 (matter dominated) because we assume a no-radiation past
  !c history for calculating t
  y(3) = 1.0
  !C THIS IS MORE APPROPRIATE BACK AT THIS EARLY TIME         
  !C AND RADIATION CAUSES EVEN MORE OF A DRAG TO BE INCLUDED
  fbclus=0.0
  fb=0.25*(sqrt(24.0*(omx+fbclus*omb+fhdmclus*omhdm)/omnr+1.0)-1.0)
  y(3) = fb

  y(4) = 0.0
  if(ldecay.eq.1) then
     y(4)=Omnrdi*0.6666667*y(1)/taud*a
  endif
  y(5)=tauin_Dlin
  ftau=0.5*y(5)/t1
  i=0
  lnorm=1
  anorm=0.001
  nn=5
  ilm = 5
  ip=0
  do i = 1, nmax
     call rkstep_Dlin(y,hh)
     Omnrd=fomnrd_Dlin(y(1))
     Omerd=y(4)/a
     if(lnorm.eq.1.and.a.ge.anorm) then
        lnorm=0
        b=exp(y(2))
        bnorm=b/a
     endif

     if(a.ge.amintab.and.a.le.amaxtab) then
        ip=ip+1
        zz=1.0/a
        hu=hub_Dlin(a)
        t=y(1)
        tau=y(5)
        ftau=0.5*tau/t1/sqrt(a)
        !c     UNITS OF T and TAU ARE 10**10 YEARS
        b=exp(y(2))
        fb=y(3)
        atab(ip)=1.0/zz
        !C THIS IS TEMPORARY STORAGE FOR D_atab
        ttab(ip)=t
        D_atab(ip)=b/a
        Dtab(ip)=b
        dlnD_dlnatab(ip)=fb
        chitab(ip)=tau
        if(a.ge.1.0.and.now.eq.1) then 
           now=0
           ip0=ip
           anow=a
           tnow=t
           taunow=tau
           ftaunow=ftau
        endif
     endif
  enddo
  ainterp=(atab(ip0)-1.0)/(atab(ip0)-atab(ip0-1))
  anow=atab(ip0)*(1.0-ainterp)+atab(ip0-1)*ainterp
  taunow=chitab(ip0)*(1.0-ainterp)+chitab(ip0-1)*ainterp
  D_agrowth=(D_atab(ip0)*(1.0-ainterp)+D_atab(ip0-1)*ainterp)/bnorm
  D_agrowth_1fac=1.0/(D_agrowth*bnorm)
  HD_Hanow=dlnD_dlnatab(ip0)*(1.0-ainterp)+dlnD_dlnatab(ip0-1)*ainterp
  ntab_DLin = ip     
  do ip=1,ntab_DLin
     D_atab(ip)=D_atab(ip)*D_agrowth_1fac
     Dtab(ip)=Dtab(ip)*D_agrowth_1fac
     chitab(ip)=3000.0*(taunow-chitab(ip))
  enddo
22 format(1x,8(1pe10.3))
12 continue
  amintab=atab(1)
  amaxtab=atab(ntab_DLin)
  D_ainfty=D_atab(1)
  D_a0=D_atab(ntab_DLin)
  HD_Hainfty=dlnD_dlnatab(1)
  HD_Ha0=dlnD_dlnatab(ntab_DLin)
  chin0=3000.0*taunow
  chininfty=6000.0/h/sqrt(Omnr)
  chimaxtab=chitab(1)
  chimintab=chitab(ntab_DLin)
  Dmintab=Dtab(1)
  Dmaxtab=Dtab(ntab_DLin)
  tmintab=ttab(1)
  tmaxtab=ttab(ntab_DLin)
!  write(*,*) 'D_ainfty,D_a0,HD_Hainfty,HD_Ha0 '
!  write(*,*) D_ainfty,D_a0,HD_Hainfty,HD_Ha0
!  write(*,*) 'chin0,chininfty,chimaxtab,chimintab '
!  write(*,*) chin0,chininfty,chimaxtab,chimintab 
!  write(*,*) 'amaxtab,amintab,Dmaxtab,Dmintab,tmaxtab,tmintab'
!  write(*,*) amaxtab,amintab,Dmaxtab,Dmintab,tmaxtab,tmintab

  return
end subroutine Dlinear_setup

subroutine rkstep_Dlin(yy,hh)
  USE intreal_types
  use count
  use qtime
  use decay
  dimension y(50),yy(50)
  dimension w1(50),w2(50),w3(50),w4(50)
  !c   ** integration var is ln(a), y1=time in 10**10y units, 
  !c     y2=fb=ln(D), y3=dlnb/dlna
  !c   **   y4=omerd*a, y5=tau, ftau*2/H0/sqrt(Omnr)*a**1.5=tau 

  do  j = 1, nn
     y(j) = yy(j)
  enddo
  call deriv_Dlin(y,w1)
  alogn = alogn + hh*0.5
  a = exp(alogn)

  do  j = 1, nn
     y(j) = yy(j) + hh*w1(j)*0.5
  enddo
  call deriv_Dlin(y,w2)

  do j = 1, nn
     y(j) = yy(j) + hh*w2(j)*0.5
  enddo
  call deriv_Dlin(y,w3)
  alogn = alogn + hh*0.5
  a = exp(alogn)

  do  j = 1, nn
     y(j) = yy(j) + hh*w3(j)
  enddo

  call deriv_Dlin(y,w4)
  do  j = 1, nn
     yy(j) = yy(j) + hh/6.*(w1(j) + w4(j) + 2.*(w2(j) + w3(j)))
  enddo

  return
end subroutine rkstep_Dlin

subroutine deriv_Dlin(y,d)
  USE intreal_types
  use qtime
  use decay
  dimension y(50),d(50)

  Omnrd=fomnrd_Dlin(y(1))
  Omerd=y(4)/a

  d(1) = 1.0/hub_Dlin(a)
  fb=y(3)
  d(2) = fb
  qq=q_DLin(a,qqf)
  d(3) = 3.0*qqf-fb*fb-(1.-qq)*fb
  d(4)=Omnrd*a*d(1)/taud
  d(5)=d(1)/a
  return
end subroutine deriv_Dlin

function fomnrd_Dlin(t)
  USE intreal_types
  use decay
  fomnrd_Dlin=0.
  if(ldecay.eq.0) return
  fomnrd_Dlin=Omnrdi*exp(-t/taud)
  return
end function fomnrd_Dlin

function q_Dlin(a,qf_Dlin)
  USE intreal_types
  use univer
  use univer2
  use h_var
  use fac
  use decay
  !c     This routine gets a couple of q=a''a/(a')^2 and qf ( which is like q without the vacuum term )
  !c     for use in updating dlnD/dlna

  q_Dlin=Omer/a+Omnr+Omvac*a**3+Omcurv*a+Omerd/a+Omnrd
  qs=0.5/q_Dlin
  q_Dlin=qs*(2.0*Omer/a+Omnr-2.0*Omvac*a**3 +2.0*Omerd/a+Omnrd)
  fnr=(Omx+OmB+fhdmclus*Omhdm)/Omnr
  if(Omb.ne.0.and.a.le.1.0e-03) fnr=(Omx+fhdmclus*Omhdm)/Omnr

  qf_Dlin=qs*(2.0*Omer/a*fer+Omnr*fnr+Omnrd*fnr)
  return
end function q_Dlin

function hub_Dlin(a)
  USE intreal_types
  use univer
  use univer2
  use h_var
  use decay
  !c     adot/a  

  hub_Dlin=Omer/a+Omnr+Omvac*a**3+Omcurv*a+Omerd/a+Omnrd
  hub_Dlin=h*sqrt(hub_Dlin/a**3)
  return
end function hub_Dlin

function tin_Dlin(a,tauin_Dlin)
  USE intreal_types
  use univer
  use univer2
  use h_var
  use decay
  use tau_vars

  !c     Start time - assumes EdeS/Radiation dominated as first approx

  aeq=(Omer+Omerd)/(Omnr+Omnrd)
  t1=1./h/sqrt(Omnr)
  t0=1./h/sqrt(Omnr+Omnrd)
  y=sqrt(a+aeq)
  y0=sqrt(aeq)
  t=t0*(0.66666667*(y**3-y0**3)-2.*aeq*(y-y0))
  tau=t0*2.*(y-y0)
  tin_Dlin=t
  tauin_Dlin=tau
  return
end function tin_Dlin


!c **   PRESS ETAL ROUTINE ******************
SUBROUTINE HUNT(XX,N,X,JLO)
  USE intreal_types
  DIMENSION XX(N)
  LOGICAL ASCND
  ASCND=XX(N).GT.XX(1)
  IF(JLO.LE.0.OR.JLO.GT.N)THEN
     JLO=0
     JHI=N+1
     GO TO 3
  ENDIF
  INC=1
  IF(X.GE.XX(JLO).EQV.ASCND)THEN
1    JHI=JLO+INC
     IF(JHI.GT.N)THEN
        JHI=N+1
     ELSE IF(X.GE.XX(JHI).EQV.ASCND)THEN
        JLO=JHI
        INC=INC+INC
        GO TO 1
     ENDIF
  ELSE
     JHI=JLO
2    JLO=JHI-INC
     IF(JLO.LT.1)THEN
        JLO=0
     ELSE IF(X.LT.XX(JLO).EQV.ASCND)THEN
        JHI=JLO
        INC=INC+INC
        GO TO 2
     ENDIF
  ENDIF
3 IF(JHI-JLO.EQ.1)RETURN
  JM=(JHI+JLO)/2
  IF(X.GT.XX(JM).EQV.ASCND)THEN
     JLO=JM
  ELSE
     JHI=JM
  ENDIF
  GO TO 3
END SUBROUTINE HUNT


subroutine lagrint_one(xa,val,est,n,x)
  USE intreal_types
  dimension xa(n),val(n)
  !C     Lagrange 4 pt interpolation or 2 point for 1 scalar

  if(n.eq.4) then
     a1=(x-xa(2))/(xa(1)-xa(2))*(x-xa(3))/(xa(1)-xa(3))*(x-xa(4))/(xa(1)-xa(4))
     a2=(x-xa(1))/(xa(2)-xa(1))*(x-xa(3))/(xa(2)-xa(3))*(x-xa(4))/(xa(2)-xa(4))
     a3=(x-xa(1))/(xa(3)-xa(1))*(x-xa(2))/(xa(3)-xa(2))*(x-xa(4))/(xa(3)-xa(4))
     a4=(x-xa(1))/(xa(4)-xa(1))*(x-xa(2))/(xa(4)-xa(2))*(x-xa(3))/(xa(4)-xa(3))
     est=a1*val(1)+a2*val(2)+a3*val(3)+a4*val(4)
  else 
     a1=(x-xa(2))/(xa(1)-xa(2))
     a2 = 1-a1
     est =a1*val(1)+a2*val(2)
  endif
  return
end subroutine lagrint_one

subroutine lagrint_Dlin(xa,val,est,n,x)
  USE intreal_types
  dimension xa(n),est(3),val(n,3)
  !C     Lagrange 4 pt interpolation or 2 point 

  if(n.eq.4) then
     a1=(x-xa(2))/(xa(1)-xa(2))*(x-xa(3))/(xa(1)-xa(3))*(x-xa(4))/(xa(1)-xa(4))
     a2=(x-xa(1))/(xa(2)-xa(1))*(x-xa(3))/(xa(2)-xa(3))*(x-xa(4))/(xa(2)-xa(4))
     a3=(x-xa(1))/(xa(3)-xa(1))*(x-xa(2))/(xa(3)-xa(2))*(x-xa(4))/(xa(3)-xa(4))
     a4=(x-xa(1))/(xa(4)-xa(1))*(x-xa(2))/(xa(4)-xa(2))*(x-xa(3))/(xa(4)-xa(3))
     do m=1,3
        est(m)=a1*val(1,m)+a2*val(2,m)+a3*val(3,m)+a4*val(4,m)
     enddo
  else 
     a1=(x-xa(2))/(xa(1)-xa(2))
     a2 = 1-a1
     do m=1,3
        est(m) =a1*val(1,m)+a2*val(2,m)
     enddo
  endif
  return
end subroutine lagrint_Dlin


