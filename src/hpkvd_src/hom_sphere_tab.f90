MODULE hom_sphere_tab
END MODULE hom_sphere_tab

!c THESE SUBROUTINES SET UP TABLES FOR vac AND open MODELS - NOTES AT END
!C   MUST BE LINKED WITH psubs_Dlinear.f

!C THE FOLLOWING LINES ARE NEEDED AT THE BEGINNING TO USE THESE ROUTINES
!C EITHER 
!c      common/univer/Omt,Omnr,Omvac,Omcurv,Rvac_nr,Rcurv_nr,Ha_b_nr
!c      common/univer2/OmB,Omx,Omer
!c      common/h0/h0
!c      common/EdeS/iEdeS 
!c      common/univercurv/iamcurved,dcurv
!c      call get_cosmology     
!C OR 
!c      read_cosmology(Omt,Omnr,Omx,OmB,h,Omcurv,Omvac,Omhdm,iamcurved,dcurv)
!c      Dlinear_cosmology(Omt,Omx,OmB,h,Omcurv,Omvac,Omhdm,iamcurved,dcurv)
!C THEN
!c      call Dlinear_setup

!C LINES FOR hom_sphere SETUP
!c      if(iEdeS.ne.1) then
!c         call make_F_zvir_table
!c      else
!c         write(*,*) ' 6 dummy lines '
!c         do i=1,6
!c            read(999,'(a)') dummy
!c         enddo
!c      endif


!C  HUNT AND LAGRINT_ONE ARE IN THE Dlinear_fns PACKAGE

!C  _________FUNCTIONS ON THE TABLE________________

function Frho_tabfn(zvir3)
  USE intreal_types
  use arrays
  use tabvec
  use EdeS
  parameter (fc_EdeS=1.686)
  dimension qx(4),qtab(4)
  if(iEdeS.eq.1) then
     Frho_tabfn=fc_EdeS*zvir3
     return
  endif
  nhunt=ntab
  call hunt(zvir3tab,nhunt,zvir3,ktab)
  a_b=1.0/zvir3
  Dlin=Dlinear(a_b,chi,HD_Ha,D_a)
  Dinv=zvir3/D_a
  IF(ktab.ge.ntab) THEN
     Frho_tabfn=FDctab(ntab)*Dinv
     return
  ELSEIF(ktab.lt.1) THEN
     Frho_tabfn=FDctab(1)*Dinv
     return
  ELSEIF(ktab.lt.2.or.ktab.eq.ntab-1) THEN
     nlint=2
     ktabl=ktab-1
  ELSE 
     nlint=4
     ktabl=ktab-2
  ENDIF
  DO ikk=1,nlint
     qx(ikk)=zvir3tab(ktabl+ikk)
     qtab(ikk)=FDctab(ktabl+ikk)
  ENDDO
  call lagrint_one(qx,qtab,qint,nlint,zvir3)
  Frho_tabfn=qint*Dinv 
  return
end function Frho_tabfn

function zvir3_tabfn(Frho)
  USE intreal_types
  use arrays
  use tabvec
  use EdeS
  parameter (fc_EdeS=1.686)
  dimension qx(4),qtab(4)
  if(iEdeS.eq.1) then
     zvir3_tabfn=Frho/fc_EdeS
     return
  endif
  nhunt=ntab
  call hunt(Ftab,nhunt,Frho,ktab)
  IF(ktab.ge.ntab) THEN
     WRITE(*,*) ' OFF UPPER END OF TABLE '
     WRITE(*,*) ' EXTRAPOLATING VALUES '
     zvir3_tabfn=Frho/Factab(ntab)
     return
  ELSEIF(ktab.lt.1) THEN
     WRITE(*,*) ' OFF LOW END OF TABLE '
     WRITE(*,*) ' EXTRAPOLATING VALUES '
     zvir3_tabfn=Frho/Factab(1)
     return
  ELSEIF(ktab.lt.2.or.ktab.eq.ntab-1) THEN
     nlint=2
     ktabl=ktab-1
  ELSE 
     nlint=4
     ktabl=ktab-2
  ENDIF
  DO ikk=1,nlint
     qx(ikk)=Ftab(ktabl+ikk)
     qtab(ikk)=zvir3tab(ktabl+ikk)
  ENDDO
  call lagrint_one(qx,qtab,qint,nlint,Frho)
  zvir3_tabfn=qint 
  return
end function zvir3_tabfn

function fsc_tabfn(Frho)
  USE intreal_types
  use arrays
  use tabvec
  use EdeS
  parameter (fc_EdeS=1.686)
  dimension qx(4),qtab(4)
  if(iEdeS.eq.1) then
     fsc_tabfn=fc_EdeS
     return
  endif
  nhunt=ntab
  call hunt(Ftab,nhunt,Frho,ktab)
  IF(ktab.ge.ntab) THEN
     fsc_tabfn=FDctab(ntab)
     return
  ELSEIF(ktab.lt.1) THEN
     fsc_tabfn=FDctab(1)
     return
  ELSEIF(ktab.lt.2.or.ktab.eq.ntab-1) THEN
     nlint=2
     ktabl=ktab-1
  ELSE 
     nlint=4
     ktabl=ktab-2
  ENDIF
  DO ikk=1,nlint
     qx(ikk)=Ftab(ktabl+ikk)
     qtab(ikk)=FDctab(ktabl+ikk)
  ENDDO
  call lagrint_one(qx,qtab,qint,nlint,Frho)
  fsc_tabfn=qint 
  return
end function fsc_tabfn

function fsc_tabfn_of_ZZ(zvir3)
  USE intreal_types
  use input_parameters
  use arrays
  use tabvec
  use EdeS

  parameter (fc_EdeS=1.686)
  dimension qx(4),qtab(4)
  if(iEdeS.eq.1) then
     fsc_tabfn_of_ZZ=fc_EdeS
     return
  endif
  nhunt=ntab
  call hunt(zvir3tab,nhunt,zvir3,ktab)
  IF(ktab.ge.ntab) THEN
     fsc_tabfn_of_ZZ=FDctab(ntab)
     return
  ELSEIF(ktab.lt.1) THEN
     fsc_tabfn_of_ZZ=FDctab(1)
     return
  ELSEIF(ktab.lt.2.or.ktab.eq.ntab-1) THEN
     nlint=2
     ktabl=ktab-1
  ELSE 
     nlint=4
     ktabl=ktab-2
  ENDIF
  DO ikk=1,nlint
     qx(ikk)=zvir3tab(ktabl+ikk)
     qtab(ikk)=FDctab(ktabl+ikk)
  ENDDO
  call lagrint_one(qx,qtab,qint,nlint,zvir3)
  fsc_tabfn_of_ZZ=qint 
  return
end function fsc_tabfn_of_ZZ

function Frho_of_fsc_tabfn(fsc)
  USE intreal_types
  use arrays
  use tabvec
  use EdeS
  parameter (fc_EdeS=1.686)
  dimension qx(4),qtab(4)
  if(iEdeS.eq.1) then
     Frho_of_fsc_tabfn=fc_EdeS*10.0
     return
  endif
  nhunt=ntab
  call hunt(FDctab,nhunt,fsc,ktab)
  IF(ktab.ge.ntab) THEN
     Frho_of_fsc_tabfn=Ftab(ntab)
     return
  ELSEIF(ktab.lt.1) THEN
     Frho_of_fsc_tabfn=Ftab(1)
     return
  ELSEIF(ktab.lt.2.or.ktab.eq.ntab-1) THEN
     nlint=2
     ktabl=ktab-1
  ELSE 
     nlint=4
     ktabl=ktab-2
  ENDIF
  DO ikk=1,nlint
     qx(ikk)=FDctab(ktabl+ikk)
     qtab(ikk)=Ftab(ktabl+ikk)
  ENDDO
  call lagrint_one(qx,qtab,qint,nlint,fsc)
  Frho_of_fsc_tabfn=qint 
  return
end function Frho_of_fsc_tabfn

function zvir3_a_3vir_tabfn(a_3vir)
  USE intreal_types
  use arrays
  use tabvec
  use EdeS
  parameter (a_3vir_a_bvir_EdeS=0.3/1.686)
  dimension a_3virtab(ntabmax)
  dimension qx(4),qtab(4)
  data lone/1/
  save lone
  if(iEdeS.eq.1) then
     zvir3_a_3vir_tabfn=a_3vir_a_bvir_EdeS/a_3vir
     return
  endif
  if(lone.eq.1) then
     do j=1,ntab
        a_3virtab(j)=a_3vir_a_bvirtab(j)/zvir3tab(j)
     enddo
     lone=0
  endif

  nhunt=ntab
  call hunt(a_3virtab,nhunt,a_3vir,ktab)
  IF(ktab.ge.ntab) THEN
     WRITE(*,*) ' OFF UPPER END OF TABLE '
     WRITE(*,*) ' EXTRAPOLATING VALUES '
     zvir3_a_3vir_tabfn=a_3vir_a_bvirtab(ntab)/a_3vir
     return
  ELSEIF(ktab.lt.1) THEN
     WRITE(*,*) ' OFF LOW END OF TABLE '
     WRITE(*,*) ' EXTRAPOLATING VALUES '
     zvir3_a_3vir_tabfn=a_3vir_a_bvirtab(1)/a_3vir
     return
  ELSEIF(ktab.lt.2.or.ktab.eq.ntab-1) THEN
     nlint=2
     ktabl=ktab-1
  ELSE 
     nlint=4
     ktabl=ktab-2
  ENDIF
  DO ikk=1,nlint
     qx(ikk)=a_3virtab(ktabl+ikk)
     qtab(ikk)=zvir3tab(ktabl+ikk)
  ENDDO
  call lagrint_one(qx,qtab,qint,nlint,a_3vir)
  zvir3_a_3vir_tabfn=qint 
  return
end function zvir3_a_3vir_tabfn

function a_3vir_tabfn(zvir3)
  USE intreal_types
  use arrays
  use tabvec
  use EdeS

  parameter (a_3vir_a_bvir_EdeS=0.3/1.686)
  dimension qx(4),qtab(4)
  if(iEdeS.eq.1) then
     a_3vir_tabfn=a_3vir_a_bvir_EdeS/zvir3
     return
  endif
  nhunt=ntab
  call hunt(zvir3tab,nhunt,zvir3,ktab)
  IF(ktab.ge.ntab) THEN
     WRITE(*,*) ' OFF UPPER END OF TABLE '
     WRITE(*,*) ' EXTRAPOLATING VALUES '
     a_3vir_tabfn=a_3vir_a_bvirtab(ntab)/zvir3
     return
  ELSEIF(ktab.lt.1) THEN
     WRITE(*,*) ' OFF LOW END OF TABLE '
     WRITE(*,*) ' EXTRAPOLATING VALUES '
     a_3vir_tabfn=a_3vir_a_bvirtab(1)/zvir3
     return
  ELSEIF(ktab.lt.2.or.ktab.eq.ntab-1) THEN
     nlint=2
     ktabl=ktab-1
  ELSE 
     nlint=4
     ktabl=ktab-2
  ENDIF
  DO ikk=1,nlint
     qx(ikk)=zvir3tab(ktabl+ikk)
     qtab(ikk)=a_3vir_a_bvirtab(ktabl+ikk)
  ENDDO
  call lagrint_one(qx,qtab,qint,nlint,zvir3)
  a_3vir_tabfn=qint/zvir3 
  return
end function a_3vir_tabfn

subroutine inithom_sphere_params
  USE intreal_types
  use arrays
  use tabvals_sph
  use params_sph
  use equil_sph
  use univer

  nstepmax=50000
  !c tfac / fraction of local 1-axis Hubble time for dt
  tfac=0.005 
  !c  zinit_fac NOTE zinit=min(zinit_fac*Frho,zinit_fac)
  !c      zinit_fac=101 
  zinit_fac=20
  iforce_strat=0
  !c ivir_strat / 1 a_jeq=fcoll_3 a_b3 / 
  ivir_strat=1
  fcoll_3=0.01

  !c critical overdensity /<= 0 stops/ 
  dcrit=2000
  ntabout=100
  nout=200000
  !c      Frho_max,Frho_min,dlogFrho
  Frho_max=20
  Frho_min=0.8
  dlogFrho=0.025
  itab=0
  Frho=Frho_max
  dFrho_increment=1.0/10.0**dlogFrho
  zinit=max(zinit_fac,zinit_fac*Frho)
  return
end subroutine inithom_sphere_params

subroutine readhom_sphere_params
  USE intreal_types
  use arrays
  use tabvals_sph
  use params_sph
  use equil_sph
  use univer

  nstepmax=50000
  write(*,*) ' tfac / fraction of local 1-axis Hubble time for dt/'
  read(999,*) tfac
  write(*,*) ' zinit_fac NOTE zinit=min(zinit_fac*Frho,zinit_fac)'
  read(999,*) zinit_fac

  iforce_strat=0
  write(*,*) ' ivir_strat / 1 a_jeq=fcoll_3 a_b3 / '
  read(999,*) ivir_strat
  write(*,*) ' fcoll_3 '
  read(999,*) fcoll_3

  write(*,*) 'critical overdensity /<= 0 stops/ '
  read(999,*) dcrit
  ntabout=100
  nout=200000
  write(*,*) 'Frho_max,Frho_min,dlogFrho'
  read(999,*) Frho_max,Frho_min,dlogFrho

  return
end subroutine readhom_sphere_params

subroutine make_F_zvir_table
  USE intreal_types
  use arrays
  use tabvals_sph
  use tabvec
  use evalues_sph
  use params_sph
  use equil_sph
  use univer

  external get_derivs_sph

  call inithom_sphere_params
  itab=0
  Frho=Frho_max
  dFrho_increment=1.0/10.0**dlogFrho

333 continue
  zinit=max(zinit_fac,zinit_fac*Frho)
  call evolve_sphere_tab(zvir3,zturn3,a_3ta,a_3vir_a_bvir,&
       E_Mta_H02rpk2,EK_Mvir_H02rpk2_corr)
  if(zvir3.gt.0.0) then
     itab=itab+1
     Ftab(itab)=Frho
     zvir3tab(itab)=zvir3
     a_b=1.0/zvir3
     Dlin=Dlinear(a_b,chi,HD_Ha,D_a)
     FDctab(itab)=Frho*D_a/zvir3
     Factab(itab)=Frho/zvir3
     a_3vir_a_bvirtab(itab)=a_3vir_a_bvir
     Frho=Frho*dFrho_increment
     if(Frho.ge.Frho_min) goto 333
  endif
  ntab=itab
  return
end subroutine make_F_zvir_table


subroutine evolve_sphere_tab(zvir3,zturn3,a_3ta,a_3vir_a_bvir,E_Mta_H02rpk2,EK_Mvir_H02rpk2_corr)
  USE intreal_types
  use arrays
  use params_sph
  use hard
  use evalues_sph
  use equil_sph
  use univer
  PARAMETER (ny=2)
  PARAMETER (one_third=1.0/3.0)
  dimension y(ny),dy(ny)
  external get_derivs_sph
  nyp=ny
  a_3eq=0.0
  lvirv=0
  t=1.0/zinit
  a3p=0.0            
  tp=0.0
  lvir_3p=0
  v3p2=0.0
  v3p=0.0
  ldcrit=0
  zturn3=-1
  zturn3p=-1
  delta1_ta=0.0
  delta1_tap=0.0
  zvir3=-1
  zdcrit=-1
  pi=4.0*atan(1.0)
  sq3=sqrt(3.0)
  delta1_ta0=9.0*pi*pi/16.0
  !c      delta1_ira_const=(2.0-1.831)/2.0
  delta1_ira_cnst=0.0845
  aLam_3 = Frho/3.0
  !C ROUGH ESTIMATE OF COLLAPSE zzc FROM aLam_1
  zzc_170sph=3.0*(aLam_3)/1.606
  zzc_Fsph=3.0*(aLam_3)/1.686

  call ic_set_sph(zinit,aLam_3,t,nyp,y)

  istep=0
111 continue
  istep=istep+1
  a_3=y(1)
  dtstep=tfac*a_3*sqrt(a_3)*Ha_b_nr
  !C NOT FUNNY RK4 OF PRESS REQUIRES FIRST CALL TO derivs
  call get_derivs_sph(t,y,dy)
  call rk4_sph(y,dy,nyp,t,dtstep,y,get_derivs_sph)
  t=t+dtstep
  delta1_e=t**3/(y(1)**3)
  if(a3p.ne.-1) then
     if(y(1).lt.a3p) then
        !c LINEAR INTERPOLATION: v3_interp=v3p+(v3-v3p)*(a_b_interp-a_bp)/(a_b-a_bp)
        !C    SET THIS TO ZERO
        a3=y(1)
        v3=y(2)
        pfac=-v3p/(v3-v3p)
        a_3ta=a3p+(a3-a3p)*pfac
        a_bta=tp+pfac*(t-tp)
        zturn3=1.0d0/a_bta
        delta1_ta=(a_bta/a_3ta)**3
        if(omvac.eq.0.0) then
           a_3vir=0.5*a_3ta
        else
           call get_cubic_soln(a_3ta,a_3vir)
        endif
        !c            a_3vir_a_3ta=a_3vir/a_3ta
        E_Mta_H02rpk2=-0.5*(Omnr/a_3ta+Omvac*a_3ta**2)*3.0/5.0
        EK_Mvir_H02rpk2_corr=0.5*(3.0*Omvac*a_3vir**2)*3.0/5.0
        !c            write(*,*) 'a_3vir,a_3vir_a_3ta = ',a_3vir,a_3vir_a_3ta
        !c            write(*,*) 'E_Mta_H02rpk2 = ',E_Mta_H02rpk2            

        !C ETA_LLPR=Lambda/(4 pi G rho_ta) = 2OMvac/OMnr/(delta1_ta*zturn3**3)
        eta_LLPR=OMvac/OMnr/delta1_ta/zturn3**3
        delta1_ta_ira=delta1_ta0*((1.0-delta1_ira_cnst*eta_LLPR)/(1.0-eta_LLPR))**sq3
        a3p=-1
        if(eta_LLPR.ne.0) then
           x_ira=2.0/eta_LLPR
        else 
           x_ira=0.0
        endif
     else 
        a3p=y(1)
        v3p2=v3p
        v3p=y(2)
        tp=t
        delta1_tap=delta1_ta
        delta1_ta=delta1_e
        zturn3p=zturn3
        zturn3=1.0d0/t
     endif
  endif
  if(lvirv.ne.lvir_3p) then
     zvir3=1.0/t
     lvir_3p=1
     a_3vir_a_bvir=a_3vir*zvir3
     delta1_e3=(t/y(1))**3
     delta1_vir=(t/a_3vir)**3
     !c         write(*,*) 'delta1_ta,delta1_vir =',delta1_ta,delta1_vir
  endif
  if(lvirv.eq.1.or.istep.ge.nstepmax.or.t.ge.2.0) return
  goto 111
end subroutine evolve_sphere_tab


subroutine ic_set_sph(zinit,aLam_3,t,nyp,y)
  USE intreal_types
  use arrays
  use univer
  use Dlin_params
  parameter (iE_Mout=0)
  dimension y(nyp)
  Rvac_nr=Omvac/Omnr
  Rcurv_nr=Omcurv/Omnr
  zz=zinit
  a_b=1.0/zz
  t=a_b
  !C THIS ASSUMES THAT WE START IN THE nr DOMINATED REGIME

  dlin=Dlinear(a_b,chi,HD_Ha,D_a)
  !C RATIO IS TAKEN HERE wrt H_0 Omnr^{1/2}
  Ha_b_nr=sqrt((1.0+Rvac_nr*a_b**3+Rcurv_nr*a_b)/a_b)
  Ha_b_nrinv=1.0/Ha_b_nr
  y(1)=a_b*(1.0-dlin*aLam_3)
  y(2)=Ha_b_nr*(1.0-(1.0+HD_Ha)*dlin*aLam_3)
  if(iE_Mout.eq.1) then
     Epert_M_H02rpk2=-0.5*D_a*aLam_3*3.0*Omnr*((3.0+2.0*HD_Ha)/5.0+Rvac_nr*a_b**3*2.0*HD_Ha/5.0 &
          +Rcurv_nr*a_b*2.0*(1.0+HD_Ha)/5.0)
     Epert_M_H02rpk2infty=-0.5*D_ainfty*aLam_3*3.0*Omnr
     Eunpert_M_H02rpk2=0.5*Omcurv*3.0/5.0
     Etot_M_H02rpk2=Eunpert_M_H02rpk2+Epert_M_H02rpk2
     Etot_M_H02rpk2infty=Eunpert_M_H02rpk2+Epert_M_H02rpk2infty
     write(*,*) ' Epert_M_H02rpk2 = ',Epert_M_H02rpk2
     write(*,*) ' Epert_M_H02rpk2infty = ',Epert_M_H02rpk2infty
     write(*,*) ' Eunpert_M_H02rpk2 = ',Eunpert_M_H02rpk2
     write(*,*) ' Etot_M_H02rpk2 = ',  Etot_M_H02rpk2
     write(*,*) ' Etot_M_H02rpk2infty = ',  Etot_M_H02rpk2infty
  endif
  return
end subroutine ic_set_sph

subroutine get_derivs_sph(t,y,dy)
  USE intreal_types
  use arrays
  use evalues_sph
  use univer
  use equil_sph
  PARAMETER (one_third=1.0/3.0)
  PARAMETER (ny=2)
  dimension y(ny),dy(ny)
  dimension avec(1)
  dimension bvec_2(1)

  avec(1)=y(1)
  a_b=t
  a_b3=a_b**3
  Ha_b_nr=sqrt((1.0+Rvac_nr*a_b3+Rcurv_nr*a_b)/a_b)
  Ha_b_nrinv=1.0/Ha_b_nr

  if(lvirv.eq.0.and.(avec(1).le.fcoll_3*a_b)) then
     lvirv=1
     a_3eq=fcoll_3*a_b
     if(ivir_strat.eq.1) then
        a_3eq2=a_3eq*1.001
        a_3eq1=a_3eq*1.0003
     endif
  endif

  if(lvirv.eq.1) then
     avec(1)=a_3eq
     y(1)=avec(1)
     y(2)=0.0
     dy(1)=0.0
     dy(2)=0.0
  endif
  delta1_e=a_b3/avec(1)**3

  if(iforce_strat.eq.0) then
     d0=-Rvac_nr*Ha_b_nrinv
     d1_int=1.5*delta1_e/a_b3*Ha_b_nrinv
  else
     d0=(0.5/a_b3-Rvac_nr)*Ha_b_nrinv
     d1_int=1.5*(delta1_e-1.0)/a_b3*Ha_b_nrinv
     if(iforce_strat.eq.4.or.iforce_strat.eq.5.or.iforce_strat.eq.6) then
        dlin=Dlinear(a_b,chi,HD_Ha,D_a)
        d1_ext=1.5/a_b3*Ha_b_nrinv*dlin
        Frho_3=Frho/3.0
        !C I THINK 3 IS FUNDAMENTALLY FLAWED
     elseif(iforce_strat.eq.3) then
        d1_ext=1.5*2.5/a_b3*Ha_b_nrinv
     endif
  endif
  if(iforce_strat.ne.5.and.iforce_strat.ne.6) then
     bvec_2(1)=one_third
  elseif(iforce_strat.eq.6) then
     sum=avec(1)
     sum_3=sum/3.0
     bvec_2(1)=one_third
     !c +0.4*(sum_3-avec(j))/sum_3
  endif
  if(lvirv.eq.0) then
     dy(1)=y(2)*Ha_b_nrinv
     if(iforce_strat.eq.0) then
        dy(2)=-y(1)*d1_int*bvec_2(1)
     elseif(iforce_strat.eq.1) then
        dy(2)=-y(1)*(d1_int*bvec_2(1)+d0)
     elseif(iforce_strat.eq.3) then
        dy(2)=-y(1)*(d1_int*bvec_2(1)+d0+d1_ext*(bvec_2(1)-one_third))
     elseif(iforce_strat.eq.4.or.iforce_strat.eq.6) then
        dy(2)=-y(1)*(d1_int*bvec_2(1)+d0+d1_ext*(aLam_3-Frho_3))
     elseif(iforce_strat.eq.5) then
        dy(2)=-y(1)*d0-d1_ext*a_b*aLam_3
     endif
  endif
  return
end subroutine get_derivs_sph


SUBROUTINE RK4_sph(Y,DYDX,N,X,H,YOUT,DERIVS)
  USE intreal_types
  use arrays
  PARAMETER (NMAX=10)
  DIMENSION Y(N),DYDX(N),YOUT(N),YT(NMAX),DYT(NMAX),DYM(NMAX)
  EXTERNAL DERIVS
  HH=H*0.5
  H6=H/6.
  XH=X+HH
  DO I=1,N
     YT(I)=Y(I)+HH*DYDX(I)
  ENDDO
  CALL DERIVS(XH,YT,DYT)
  DO I=1,N
     YT(I)=Y(I)+HH*DYT(I)
  ENDDO
  CALL DERIVS(XH,YT,DYM)
  DO I=1,N
     YT(I)=Y(I)+H*DYM(I)
     DYM(I)=DYT(I)+DYM(I)
  ENDDO
  CALL DERIVS(X+H,YT,DYT)
  DO I=1,N
     YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
  ENDDO
  RETURN
END SUBROUTINE RK4_sph

subroutine get_cubic_soln(a_3ta,a_3vir)
  USE intreal_types
  use arrays
  use univer
  use cubic_params
  LOGICAL Lzbrac
  external a_3v_cubic
  f_taX2=2.0*(1.0/a_3ta+rvac_nr*a_3ta*a_3ta)
  tol=1.0e-08
  a_3vU=a_3ta
  a_3vL=1.0/f_taX2
  !c      tolL=a_3v_cubic(a_3vL)
  !c      tolU=a_3v_cubic(a_3vU)
  !c      write(*,*) 'a_3vL,tolL,a_3vU,tolU ',a_3vL,tolL,a_3vU,tolU
  call zbrac(a_3v_cubic,a_3vL,a_3vU,Lzbrac)
  if(.not.Lzbrac) then
     write(*,*) ' no bracketing in a_3v_cubic, for vir; stopping '
     stop
  endif
  a_3vir=zbrent(a_3v_cubic,a_3vL,a_3vU,tol)
  !c      tol_vir=a_3v_cubic(a_3vir)
  !c      write(*,*) ' a_3vir,tol_vir ',a_3vir,tol_vir
  return
end subroutine get_cubic_soln

function a_3v_cubic(a_3)
  USE intreal_types
  use arrays
  use univer
  use cubic_params
  !c  OLD         a_3v_cubic=1.0/a_3+Rvac_nr*a_3*a_3-f_taX2
  a_3v_cubic=1.0/a_3+4.0*Rvac_nr*a_3*a_3-f_taX2
  return
end function a_3v_cubic

SUBROUTINE ZBRAC(FUNC,X1,X2,SUCCES)
  USE intreal_types
  use arrays
  PARAMETER (FACTOR=1.6,NTRY=50)
  external FUNC
  LOGICAL SUCCES
  IF(X1.EQ.X2) THEN
     write(*,*) 'You have to guess an initial range for ZBRAC'
     stop
  ENDIF
  F1=FUNC(X1)
  F2=FUNC(X2)
  SUCCES=.TRUE.
  DO J=1,NTRY
     IF(F1*F2.LT.0.)RETURN
     IF(ABS(F1).LT.ABS(F2))THEN
        X1=X1+FACTOR*(X1-X2)
        F1=FUNC(X1)
     ELSE
        X2=X2+FACTOR*(X2-X1)
        F2=FUNC(X2)
     ENDIF
  ENDDO
  SUCCES=.FALSE.
  RETURN
END SUBROUTINE ZBRAC

FUNCTION ZBRENT(FUNC,X1,X2,TOL)
  USE intreal_types
  use arrays
  PARAMETER (ITMAX=100,EPS=3.E-8)
  external FUNC
  A=X1
  B=X2
  FA=FUNC(A)
  FB=FUNC(B)
  IF(FB*FA.GT.0.) THEN
     write(*,*) 'Root must be bracketed for ZBRENT'
     stop
  ENDIF
  FC=FB
  DO ITER=1,ITMAX
     IF(FB*FC.GT.0.) THEN
        C=A
        FC=FA
        D=B-A
        E=D
     ENDIF
     IF(ABS(FC).LT.ABS(FB)) THEN
        A=B
        B=C
        C=A
        FA=FB
        FB=FC
        FC=FA
     ENDIF
     TOL1=2.*EPS*ABS(B)+0.5*TOL
     XM=.5*(C-B)
     IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
        ZBRENT=B
        RETURN
     ENDIF
     IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
        S=FB/FA
        IF(A.EQ.C) THEN
           P=2.*XM*S
           Q=1.-S
        ELSE
           Q=FA/FC
           R=FB/FC
           P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
           Q=(Q-1.)*(R-1.)*(S-1.)
        ENDIF
        IF(P.GT.0.) Q=-Q
        P=ABS(P)
        IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
           E=D
           D=P/Q
        ELSE
           D=XM
           E=D
        ENDIF
     ELSE
        D=XM
        E=D
     ENDIF
     A=B
     FA=FB
     IF(ABS(D) .GT. TOL1) THEN
        B=B+D
     ELSE
        B=B+SIGN(TOL1,XM)
     ENDIF
     FB=FUNC(B)
  ENDDO
  WRITE(*,*) 'ZBRENT exceeding maximum iterations.'
  ZBRENT=B
  RETURN
END FUNCTION ZBRENT

!C  NOTE THAT THE SCALING IS SUCH THAT THE 
!C  VELOCITY IS DIVIDED BY H_0 Omnr^{1/2} 
!C   INTEGRATION IS wrt a HENCE THERE IS ALSO AN adot 
!C   IN THE DENOMINATOR, SO WE HAVE 
!C    dv\tau_0/da  = \tau_0^-2 4pi G \rho_nr X STUFF = (3/2)a_b^-3 X STUFF

!C `TIME' INTEGRATION VARIABLE t=a_background
!c    COSMIC TIME UNIT= (8\pi G/3 rho_bnr*)^{-1/2}

!C  a_3dot=[p_3](Ha)^-1 
!C  p_3dot=-(3/4)[delta1_e /3 ]a_b^-3 a_3 (Ha)^-1
!C  delta1_e = (a_b/a_3)**3
!C  dot wrt a_b  .... 
!C  Ha=sqrt(1+OMvac/Omnra_b^3+Omcurv/Omnr a_b) /a_b^{1/2}
!C INSTEAD WE USE THE BINNEY-TREMAINE FORMULA FOR b_i





