MODULE homeosubs
END MODULE homeosubs

  !c now to run this, one must set up call 

  !c Homogeneous Ellipsoid routines for hpkvd.f p3mhpkvd.f etc.
  !c Version 21 May 1997

  !C `TIME' INTEGRATION VARIABLE t=a_background
  !c    COSMI!C TIME UNIT= (8\pi G/3 rho_b*)^{-1/2}

  !C  a_idot=[p_i](Ha)^-1 
  !C  p_idot=-(3/4)[delta1_e b_i + (2/3 - b_i)]a_b^-3 a_i (Ha)^-1
  !C  delta1_e = (a_b/a_1)(a_b/a_2)(a_b/a_3)
  !C  dot wrt a_b  .... 
  !C  Ha=sqrt(1+OMvac/Omnra_b^3+Omcurv/Omnr a_b) /a_b^{1/2}
  !C  b_idot=[-2p_i/a_i + SUM_{j.ne.i} [b_i-b_j]/[a_i^2-a_j^2] 
  !C     (a_ip_i-a_jp_j) + b_i SUM_j p_j/a_j](Ha)^-1
  !C INSTEAD WE USE THE BINNEY-TREMAINE FORMULA FOR b_i

subroutine read_inithom_ellipse
  USE intreal_types
  use arrays
  use wantevpvtab
  use evpvtable
  use params
  use equil 
  use evalues
  use hard
  use univer
  use rd
  use input_parameters
 
  parameter (ny=6)
  dimension quicktab(0:1,0:1)
  character(LEN=80) zevoutfile
  character(LEN=120) dummy
  nstepmax=100000
  nyp=ny
  write(*,*) ' do you want rd for elliptic?/hardwired 1/ '
  read(999,*) iwant_rd
  write(*,*) 'do you want turnaround & virialization z vs e_v? (1)'
  write(*,*) ' ... (2) makes z vs p_v for fixed e_v '
  write(*,*) ' ... (3) makes table '
  read(999,*) iwant_evmap
  write(*,*) 'want an e_v-p_v table?/dum if iwant_evmap=0/'
  read(999,*) iwant_evpvtab
  write(*,*) ' iforce_strat / 0 no bgnd, 1 std bgnd, '
  write(*,*) '.... 3 bgnd + nonlinear strain, 4 std bgnd ','+ linear strain, '
  write(*,*) '.... 5 linear strain only, 6 SW apprx for b_i /'
  read(999,*) iforce_strat
  write(*,*) 'tfac /fraction of local 1-axis Hubble time for dt/.01'
  read(999,*) tfac
  write(*,*) 'zinit / redshift + 1 /'
  read(999,*) zinit
  write(*,*) ' ivir_strat / 1 a_jeq=fcoll_3 a_b3, ',' 2 a_jeq=fcoll_j a_bj/'
  read(999,*) ivir_strat
  write(*,*) ' fcoll_3 / (=fcoll_1,fcoll_2) /'
  read(999,*) fcoll_3
  fcoll_1=fcoll_3
  fcoll_2=fcoll_3
  write(*,*) 'critical overdensity '
  read(999,*) dcrit
  write(*,*) ' nout / output every nout steps / '
  read(999,*) nout 
  return
end subroutine read_inithom_ellipse

subroutine inithom_ellipse_redshift(Zb)
  USE intreal_types
  use params
  use input_parameters
  use mpivars

  zinit=20.0*Zb
  a_bcrit=1.0/Zb

  return
end subroutine inithom_ellipse_redshift

subroutine makehom_ellipse_tab
  USE intreal_types
  use arrays
  use evalues
  use hard
  use wantevpvtab
  use params
  use evpvtable
  use input_parameters
  use mpivars

  character(LEN=120) fcevpvoutfile

  if(Nfsc.gt.1) then
     fscmin=fsc_tabfn(Fsmin)*1.00001
     fscmax=fsc_tabfn(Fsmax)*0.99999
     if(fscmin.ge.fscmax) then 
        Nfsc=1 
        dfsc=0.0
     else
        dfsc=(fscmax-fscmin)/(Nfsc-1)
     endif
     do ktab=1,Nfsc
        fsctabv(ktab)=fscmin+dfsc*(ktab-1) 
        Ftabv(ktab)=Frho_of_fsc_tabfn(fsctabv(ktab))
     enddo
  else
     fc_EdeS=1.686
     fscmin=fc_EdeS
     fscmax=fc_EdeS
     fsctabv(1)=fc_EdeS
     Ftabv(1)=Fbar 
     dfsc=0.0
  endif
  iwant_evmap=3
  Nev=evtabmax/devtab
  Nevtot=Nev+1
  Npv_ev=pv_evtabmax/dpv_evtab
  Npev_evtot=2*Npv_ev+1
  dpv_evtab=pv_evtabmax/float(Npv_ev)
  do ktab=1,Nfsc
     do jtab=0,Nev
        do itab=-Npv_ev,Npv_ev
           fc1mat(itab,jtab,ktab)=-1.0
        enddo
     enddo
  enddo

  idynax=1
  do ktab=1,Nfsc
     Frho=Ftabv(ktab)
     zinit_fac_ell=20.0
     zinit=max(zinit_fac_ell,zinit_fac_ell*Frho)
     a_bcrit=3.0

     jtab=0
     e_v=devtab*jtab
     p_v=0.0
     call evolve_ellipse_full(idynax,zdynax,Ddynax,fcdynax)
     if(fcdynax.ne.-1.0) then
        do itab=-Npv_ev,Npv_ev
           if(ifc1inv.eq.1) then
              fc1mat(itab,jtab,ktab)=1.0/fcdynax
           else
              fc1mat(itab,jtab,ktab)=fcdynax
           endif
        enddo
     endif

     do jtab=1,Nev
        e_v=devtab*jtab
        if(e_v.gt.evtabmax) goto 119
        do itab=-Npv_ev,Npv_ev
           p_v=e_v*float(itab)*dpv_evtab
           call evolve_ellipse_full(idynax,zdynax,Ddynax,fcdynax)
           if(fcdynax.ne.-1.0) then
              if(ifc1inv.eq.1) then
                 fc1mat(itab,jtab,ktab)=1.0/fcdynax
              else
                 fc1mat(itab,jtab,ktab)=fcdynax
              endif
           endif
        enddo
119     continue
     enddo
  enddo

  if(ihard.eq.1) then
     open(10,file=fcevpvoutfile,form='unformatted',status='unknown')
     write(10) Nfsc,ifc1inv,Nev,Npv_ev,devtab,dpv_evtab,evtabmax,pv_evtabmax,Fsmin,Fsmax,zinit_fac_ell, &
          (fsctabv(ktab),ktab=1,Nfsc),(Ftabv(ktab),ktab=1,Nfsc),(((fc1mat(it,jt,kt),it=-Npv_ev,Npv_ev),jt=0,Nev),kt=1,Nfsc)
     close(10)
  elseif(ihard.eq.2) then
     open(10,file=fcevpvoutfile,form='formatted',status='unknown')
     do ktab=1,Nfsc
        write(10,*) 'fsc, Ftab ',fsctabv(ktab),Ftabv(ktab)
        do jt=0,Nev
           e_v=(jt)*devtab
           write(10,*) 'e_v ',e_v
           write(10,*) (fc1mat(it,jt,ktab),it=-Npv_ev,Npv_ev)
        enddo
     enddo
     close(10)
  endif

  return
end subroutine makehom_ellipse_tab

subroutine inithom_ellipse_params
  USE intreal_types
  use arrays
  use evpvtable
  use params
  use equil
  use evalues
  use univer
  use rd
  use wantevpvtab
  use input_parameters
  
  external get_derivs

  nstepmax=100000
  iwant_rd=1
  tfac=0.01
  iwant_evmap=0
  !c  iforce_strat / 0 no bgnd, 1 std bgnd, 
  !c        3 bgnd + nonlinear strain, 4 std bgnd + linear strain /'
  iforce_strat=4
  !c   ivir_strat / 1 a_jeq=fcoll_3 a_b3,2 a_jeq=fcoll_j a_bj/'

  ivir_strat=2
  !c fcoll_3 / (=fcoll_1,fcoll_2) /'
  fcoll_3=0.171
  !c      fcoll_3=0.01
  fcoll_2=fcoll_3
  fcoll_1=0.01
  !c   critical overdensity 
  dcrit=180.0

  nout=10000

  return
end subroutine inithom_ellipse_params

subroutine readhom_ellipse_tab
  USE intreal_types
  use arrays
  use input_parameters
  use evpvtable

  do ktab=1,Nfscmax
     do jtab=0,Nevmax
        do itab=Npv_evmin,Npv_evmax
           fc1mat(itab,jtab,ktab)=-1.0
        enddo
     enddo
  enddo
  open(10,file=etabfile,form='unformatted',status='old')
  read(10) Nfsc,ifc1inv,Nev,Npv_ev,devtab,dpv_evtab,evtabmax,pv_evtabmax,Fsmin,Fsmax,zinit_fac_ell, &
       (fsctabv(ktab),ktab=1,Nfsc),(Ftabv(ktab),ktab=1,Nfsc),(((fc1mat(it,jt,kt),it=-Npv_ev,Npv_ev),jt=0,Nev),kt=1,Nfsc)
  close(10)

  nout=10000
  return
  do ktab=1,Nfsc
     write(0,*) fsctabv(ktab),Ftabv(ktab)
     do jtab=0,Nev
        write(0,*) (fc1mat(itab,jtab,ktab),itab=-Npv_ev/2,Npv_ev/2)
     enddo
  enddo

  return
end subroutine readhom_ellipse_tab

subroutine gethom_ellipse(fcvir1,Dvir1,zvir1,iflag,iuseinterp)
  USE intreal_types
  use input_parameters
  use evpvtable
  use params
  use equil
  use evalues
  use hard
  use univer
  use rd
  use timing_diagnostics

  !C POSSIBILITIES: INTERPOLATE IN fc1,fc1inv OLD z1 IS NO GOOD 

  !C `TIME' INTEGRATION VARIABLE t=a_background
  !c    COSMIC TIME UNIT= (8\pi G/3 rho_b*)^{-1/2}

  !C  a_idot=[p_i](Ha)^-1 
  !C  p_idot=-(3/4)[delta1_e b_i + (2/3 - b_i)]a_b^-3 a_i (Ha)^-1
  !C  delta1_e = (a_b/a_1)(a_b/a_2)(a_b/a_3)
  !C  dot wrt a_b  .... 
  !C  Ha=sqrt(1+OMvac/Omnra_b^3+Omcurv/Omnr a_b) /a_b^{1/2}
  !C  b_idot=[-2p_i/a_i + SUM_{j.ne.i} [b_i-b_j]/[a_i^2-a_j^2] 
  !C     (a_ip_i-a_jp_j) + b_i SUM_j p_j/a_j](Ha)^-1
  !C INSTEAD WE USE THE BINNEY-TREMAINE FORMULA FOR b_i

  !c  return iflag = -1    ---   bad peak e_v,p_v out of bounds
  !c         iflag =  0    ---   direct calculation
  !c         iflag =  1    ---   table interp


  dimension quicktab(0:1,0:1)

  real(sp) fcvir1v(2)
  !c      character(LEN=80) zevoutfile
  external get_derivs

  iflag=-1
  zvir1=-1
  fcvir1=-1
  if(Frho.lt.0.) return
  !c         Frho = aLam_1 + aLam_2 + aLam_3
  !c         e_v = 0.5*(aLam_3 - aLam_1)/Frho
  !c         P_v = 0.5*(aLam_1 + aLam_3 - 2.0*aLam_2)/Frho
  !c Note these are in INCREASING order
  aLam_3 = Frho/3.0*(1.0+3.0*e_v+p_v)
  aLam_2 = Frho/3.0*(1.0-2.0*p_v)
  aLam_1 = Frho/3.0*(1.0-3.0*e_v+p_v)

  !c flag extreme values
  if(e_v.gt.evtabmax.or.e_v.lt.0.0) return
  if(abs(p_v).gt.e_v) return
  if(igethom_strat.eq.2) then
     call evolve_ellipse(1,zvir1,Dvir1,fcvir1)
     iflag=0
     return
  endif

  iuseinterp=1

  if(Nfsc.eq.0) then     
     iuseinterp=0  
     ifsctab=1
     ifsctab_1=0
  elseif(Nfsc.eq.1) then
     ifsctab=1
     ifsctab_1=0
     nFtabs=1
     nFtabs_1=nFtabs-1
  elseif(Nfsc.gt.1) then
     fsc=fsc_tabfn(Frho)
     nFtabs=2
     if(fsc.lt.fsctabv(1)) then
        nFtabs=1
        ifsctab=1
        iuseinterp=1 
     elseif(fsc.ge.fsctabv(Nfsc)) then
        nFtabs=1
        ifsctab=Nfsc
        iuseinterp=1 
     else
        do ifsc=1,Nfsc-1
           if(fsc.ge.fsctabv(ifsc).and.fsc.lt.fsctabv(ifsc+1)) then
              ifsctab=ifsc
              goto 1199
           endif
        enddo
1199    continue
     endif
     ifsctab_1=ifsctab-1
     fscL=fsctabv(ifsctab)
     if(nFtabs.eq.1) then
        fscU=fscL
     else
        fscU=fsctabv(ifsctab+1)
     endif
  endif
  if(iuseinterp.eq.1) then
     if(e_v.ne.0.0) then
        if(abs(p_v).gt.e_v) return
        ritab=(p_v/e_v)*Npv_ev
        itab_1=ritab
        if(ritab.lt.0.0) itab_1=itab_1-1
        itab=itab_1+1
        diffi=ritab-itab_1
        rjtab=e_v/devtab
        jtab_1=rjtab
        diffj=rjtab-jtab_1
        jtab=jtab_1+1
        if(jtab.gt.Nev) then
           iuseinterp=0
           goto 1120
        endif
        ibadinterp=0
        do jfsc=1,nFtabs
           ibad=0
           jfsctab=ifsctab_1+jfsc
           do jt=0,1
              do it=0,1
                 quicktab(it,jt)=fc1mat(itab_1+it,jtab_1+jt,jfsctab)
                 if(quicktab(it,jt).le.0.0) then
                    ibad=ibad+1
                 endif
              enddo
           enddo
           if(ibad.gt.0) then
              if(ibad.lt.4) then
                 iuseinterp=0
                 goto 1120
              else
                 ibadinterp=ibadinterp+1
                 fcvir1v(jfsc)=-1.0
              endif
           else
              fc1app3=(quicktab(0,0)+(diffi*(quicktab(1,0)-quicktab(0,0))+diffj*(quicktab(0,1)-quicktab(0,0))))
              !C     3rd ORDER INTERPOLATION
              fc1app4=fc1app3+(diffi*diffj*((quicktab(1,1)-quicktab(0,1))-(quicktab(1,0)-quicktab(0,0))))
              !C 4th ORDER INTERPOLATION
              fcvir1v(jfsc)=fc1app4
           endif
        enddo
     else
        ibadinterp=0
        do jfsc=1,nFtabs
           jfsctab=ifsctab_1+jfsc
           fcvir1v(jfsc)=fc1mat(0,0,jfsctab)
           if(fcvir1v(jfsc).eq.-1.0) ibadinterp=ibadinterp+1
        enddo
     endif
     if(ibadinterp.gt.0) then
        iuseinterp=0
        if(ibadinterp.eq.nFtabs) then
           return
        else
           goto 1120
        endif
     endif
     if(nFtabs.gt.1) then
        qq=(fsc-fscL)/(fscU-fscL)
        !c IF ifc1inv.eq.1, it interpolates in fsc^-1
        if(ifc1inv.eq.1) qq=qq*(fscU/fsc)
        fcvir1=fcvir1v(1)+(fcvir1v(2)-fcvir1v(1))*qq
     else 
        fcvir1=fcvir1v(1)
     endif
     if(ifc1inv.eq.1) fcvir1=1.0/fcvir1
     Dvir1=fcvir1/Frho
     zvir1=1.0/afnofD(Dvir1)
!     if(Dvir1.ge.1.0) write(*,*) 'fcvir1,Frho,e_v,p_v ',fcvir1,Frho,e_v,p_v,nFtabs,fsc,fscL,fscU,zvir1
!     if(Dvir1.ge.1.0.and.nFtabs.eq.1) write(*,*) 'nFtabs,fsc,fscL,fscU',nFtabs,fsc,fscL,fscU
     iflag=1
  endif
1120 continue
  if(iuseinterp.eq.0) then
     zinit_fac_ell=20.0
     zinit=max(zinit_fac_ell,zinit_fac_ell*Frho)
     return
     call evolve_ellipse(1,zvir1,Dvir1,fcvir1)
     iflag=0
  endif

  return
end subroutine gethom_ellipse

subroutine evolve_ellipse(idynax,zdynax,Ddynax,fcdynax)
  USE intreal_types
  use params
  use bvalues
  use evalues
  use equil
  use univer

  PARAMETER (ny=6)
  PARAMETER (one_third=1.0/3.0)
  parameter (nzdynmax=7)
  dimension y(ny),dy(ny)

  dimension zdynv(nzdynmax),delta1_ev(nzdynmax),aturnv(3),dcritv(3)
  integer(i4b) ldynv(nzdynmax),ldyn_pv(nzdynmax)
  external get_derivs
  no_vir=3
  no_turn=0
  no_dens=1
  nzdyn=no_vir+no_turn+no_dens
  ivir0=1
  ivir1=no_vir
  iturn0=no_vir+1
  iturn1=no_vir+no_turn
  idens0=no_vir+no_turn+1
  idens1=no_vir+no_turn+no_dens
  dcritv(1)=dcrit

  nyp=ny
  a_3eq=0.0
  a_2eq=0.0
  a_1eq=0.0
  lvirv(1)=0
  lvirv(2)=0
  lvirv(3)=0
  !c      do i=1,3
  !c         aturnv(i)=0.0
  !c     enddo
  zdynax=-1
  Ddynax=-1
  fcdynax=-1
  do ii=1,nzdyn
     ldynv(ii)=0
     delta1_ev(ii)=-1.0
     zdynv(ii)=-1.0
  enddo

  if(Frho.lt.0.) return

  !c      Frho = aLam_1 + aLam_2 + aLam_3
  !c      e_v = 0.5*(aLam_3 - aLam_1)/Frho
  !c      P_v = 0.5*(aLam_1 + aLam_3 - 2.0*aLam_2)/Frho
  !c Note these are in INCREASING order
  aLam_3 = Frho/3.0*(1.0+3.0*e_v+p_v)
  aLam_2 = Frho/3.0*(1.0-2.0*p_v)
  aLam_1 = Frho/3.0*(1.0-3.0*e_v+p_v)

  zzc_180sph=(aLam_1+aLam_2+aLam_3)/1.606
  zzc_Fsph=(aLam_1+aLam_2+aLam_3)/1.686

  call ic_set(zinit,aLam_1,aLam_2,aLam_3,t,nyp,y)

  istep=0
111 continue
  istep=istep+1

  a_3=y(3)
  dtstep=tfac*a_3*sqrt(a_3)*Ha_b_nr
  call get_derivs(t,y,dy)
  call rk4(y,dy,nyp,t,dtstep,y,get_derivs)
  t=t+dtstep
  delta1_e=t**3/(y(1)*y(2)*y(3))

  do ii=ivir0,ivir1
     if(lvirv(ii).ne.ldynv(ii)) then
        zdynv(ii)=1.0/t
        ldynv(ii)=1
        !c         delta1_ev(ii)=t**3/(y(1)*y(2)*y(3))
     endif
  enddo

  do ii=iturn0,iturn1
     if(ldynv(ii).ne.1) then
        iax=ii-ivir1
        if(y(iax).lt.aturnv(iax)) then
           zdynv(ii)=1.0/t
           ldynv(ii)=1
           delta1_ev(ii)=delta1_e
        else
           aturnv(iax)=y(iax)
        endif
     endif
  enddo

  do ii=idens0,idens1
     if(ldynv(ii).ne.1) then
        if(delta1_e.ge.dcritv(ii-iturn1)) then
           zdynv(ii)=1.0/t
           ldynv(ii)=1
        endif
     endif
  enddo

  if(ldynv(idynax).eq.1.or.istep.ge.nstepmax.or.t.ge.a_bcrit) then
     zdynax=zdynv(idynax)
     if(zdynax.eq.-1) then
        Ddynax=-1
        fcdynax=-1
     else
        Ddynax=Dfnofa(1.0/zdynax)
        fcdynax=Frho*Ddynax
     endif
!     write(*,*) zdynax,istep,nstepmax
!     write(*,66) zdynv(3),zdynv(2),zdynv(1),zdynv(idens0),zzc_180sph,zzc_Fsph
66   format('zvir,dcr,zzc_180sph,zzc_Fsph : ',6(f9.3))
     return
  endif

  goto 111

end subroutine evolve_ellipse


subroutine ic_set(zinit,aLam_1,aLam_2,aLam_3,t,nyp,y)
  USE intreal_types
  use univer
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
  y(3)=a_b*(1.0-dlin*aLam_3)
  y(2)=a_b*(1.0-dlin*aLam_2)
  y(1)=a_b*(1.0-dlin*aLam_1)
  y(6)=Ha_b_nr*(1.0-(1.0+HD_Ha)*dlin*aLam_3)
  y(5)=Ha_b_nr*(1.0-(1.0+HD_Ha)*dlin*aLam_2)
  y(4)=Ha_b_nr*(1.0-(1.0+HD_Ha)*dlin*aLam_1)
  return
end subroutine ic_set

subroutine get_derivs(t,y,dy)
  USE intreal_types
  use bvalues
  use evalues
  use univer
  use equil

  PARAMETER (one_third=1.0/3.0)
  PARAMETER (ny=6)
  dimension y(ny),dy(ny)
  dimension avec(3)

  do j=1,3
     avec(j)=y(j)
  enddo
  !C THESE b_i ARE WHITE-SILK alpha_i THAT THEY SOLVE BY ODE
  a_b=t
  a_b3=a_b**3
  Ha_b_nr=sqrt((1.0+Rvac_nr*a_b3+Rcurv_nr*a_b)/a_b)
  Ha_b_nrinv=1.0/Ha_b_nr

  if(lvirv(3).eq.0.and.(avec(3).le.fcoll_3*a_b)) then
     lvirv(3)=1
     a_3eq=fcoll_3*a_b
     if(ivir_strat.eq.1) then
        a_3eq2=a_3eq*1.001
        a_3eq1=a_3eq*1.0003
     endif
  endif
  if(ivir_strat.eq.2.and.lvirv(2).eq.0) a_3eq2=fcoll_2*a_b
  if(ivir_strat.eq.2.and.lvirv(1).eq.0) a_3eq1=fcoll_1*a_b

  if(lvirv(3).eq.1.and.(avec(2).le.a_3eq2.and.lvirv(2).eq.0)) then
     lvirv(2)=1
     a_2eq=a_3eq2
  endif
  if(lvirv(3).eq.1.and.(avec(1).le.a_3eq1.and.lvirv(1).eq.0)) then
     lvirv(1)=1
     a_1eq=a_3eq1
  endif

  if(lvirv(3).eq.1) then
     avec(3)=a_3eq
     y(3)=avec(3)
     y(6)=0.0
     dy(3)=0.0
     dy(6)=0.0
  endif
  if(lvirv(2).eq.1) then
     avec(2)=a_2eq
     y(2)=avec(2)
     y(5)=0.0
     dy(2)=0.0
     dy(5)=0.0
  endif
  if(lvirv(1).eq.1) then
     avec(1)=a_1eq
     y(1)=avec(1)
     y(4)=0.0
     dy(1)=0.0
     dy(4)=0.0
  endif
  delta1_e=a_b3/avec(1)/avec(2)/avec(3)

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
     call get_b_2(avec,bvec_2)
     !C THESE b_i ARE WHITE-SILK alpha_i THAT THEY SOLVE BY ODE
  elseif(iforce_strat.eq.6) then
     sum=0.0
     do j=1,1,3
        sum=sum+avec(j)
     enddo
     sum_3=sum/3.0
     do j=1,3
        bvec_2(j)=1.0/3.0
        !c +0.4*(sum_3-avec(j))/sum_3
     enddo
     bvec_2(1)=bvec_2(1)+0.4*dlin*(aLam_1-Frho_3)
     bvec_2(2)=bvec_2(2)+0.4*dlin*(aLam_2-Frho_3)
     bvec_2(3)=bvec_2(3)+0.4*dlin*(aLam_3-Frho_3)
     !C WS approximation (SUM a^-1) doesn't work, SUM a doesn't work,
     !C     AND linear CORRECTION ISN'T QUITE ASYMMETRIC ENOUGH
     !C     BUT ISN'T BAD
  endif
  if(lvirv(3).eq.0) then
     dy(3)=y(6)*Ha_b_nrinv
     if(iforce_strat.eq.0) then
        dy(6)=-y(3)*d1_int*bvec_2(3)
     elseif(iforce_strat.eq.1) then
        dy(6)=-y(3)*(d1_int*bvec_2(3)+d0)
     elseif(iforce_strat.eq.3) then
        dy(6)=-y(3)*(d1_int*bvec_2(3)+d0+d1_ext*(bvec_2(3)-one_third))
     elseif(iforce_strat.eq.4.or.iforce_strat.eq.6) then
        dy(6)=-y(3)*(d1_int*bvec_2(3)+d0+d1_ext*(aLam_3-Frho_3))
     elseif(iforce_strat.eq.5) then
        dy(6)=-y(3)*d0-d1_ext*a_b*aLam_3
     endif
  endif
  if(lvirv(2).eq.0) then
     dy(2)=y(5)*Ha_b_nrinv
     if(iforce_strat.eq.0) then
        dy(5)=-y(2)*d1_int*bvec_2(2)
     elseif(iforce_strat.eq.1) then
        dy(5)=-y(2)*(d1_int*bvec_2(2)+d0)
     elseif(iforce_strat.eq.3) then
        dy(5)=-y(2)*(d1_int*bvec_2(2)+d0+d1_ext*(bvec_2(2)-one_third))
     elseif(iforce_strat.eq.4.or.iforce_strat.eq.6) then
        dy(5)=-y(2)*(d1_int*bvec_2(2)+d0+d1_ext*(aLam_2-Frho_3))
     elseif(iforce_strat.eq.5) then
        dy(5)=-y(2)*d0-d1_ext*a_b*aLam_2
     endif
  endif
  if(lvirv(1).eq.0) then
     dy(1)=y(4)*Ha_b_nrinv
     if(iforce_strat.eq.0) then
        dy(4)=-y(1)*d1_int*bvec_2(1)
     elseif(iforce_strat.eq.1) then
        dy(4)=-y(1)*(d1_int*bvec_2(1)+d0)
     elseif(iforce_strat.eq.3) then
        dy(4)=-y(1)*(d1_int*bvec_2(1)+d0+d1_ext*(bvec_2(1)-one_third))
     elseif(iforce_strat.eq.4.or.iforce_strat.eq.6) then
        dy(4)=-y(1)*(d1_int*bvec_2(1)+d0+d1_ext*(aLam_1-Frho_3))
     elseif(iforce_strat.eq.5) then
        dy(4)=-y(1)*d0-d1_ext*a_b*aLam_1
     endif
  endif
  return
end subroutine get_derivs

subroutine get_b_2(avec,bvec_2)
  USE intreal_types
  use rd
  PARAMETER (one_third=1.0/3.0)
  dimension avec(3),bvec_2(3)
  !C THIS ROUTINE RETURNS 
  !C   b_i/2=0.5*a_1a_2a_3 \int du [a_i^2+u] Prod_j [a_j^2+u]^{1/2}
  !C   SO THE ELLIPSOID GRAV POT IS 4PI G RHO_E (b_i/2.0) X_I^2/2  
  !C   IN THE PRINCIPAL AXIS SYSTEM
  r_2=avec(2)/avec(1)
  r_3=avec(3)/avec(1)
  if(r_3.ge.r_2) then
     if(r_3.ge.0.999) then
        bvec_2(1)=one_third
        bvec_2(2)=one_third
        bvec_2(3)=one_third
        return
     elseif(r_3.le.0.001) then
        ecc2=1.0-r_3**2
        one_ecc2=r_3**2
        ecc=sqrt(ecc2)
        one_ecc=0.5*r_3**2
        bvec_2(3)=0.5*((one_ecc2)/ecc2*(1.0/(one_ecc2)-log((1.0+ecc)/(one_ecc))*0.5/ecc))
        bvec_2(2)=bvec_2(3)
        bvec_2(1)=1.0-(bvec_2(2)+bvec_2(3))
        return
     else
        ecc2=1.0-r_3**2
        ecc=sqrt(ecc2)
        bvec_2(3)=0.5*((1.0-ecc2)/ecc2*(1.0/(1.0-ecc2)-log((1.0+ecc)/(1.0-ecc))*0.5/ecc))
        bvec_2(2)=bvec_2(3)
        bvec_2(1)=1.0-(bvec_2(2)+bvec_2(3))
        return
     endif
  endif
  if(iwant_rd.eq.1) then
     r_22=r_2*r_2
     r_32=r_3*r_3
     qrat=one_third*r_2*r_3
     bvec_2(3)=qrat*elliptic_rd(r_22,1.0,r_32)
     bvec_2(2)=qrat*elliptic_rd(r_32,1.0,r_22)
     bvec_2(1)=1.0-(bvec_2(2)+bvec_2(3))
  else
     sinth=sqrt(1.0-r_3**2)
     x = sinth/r_3
     rk2=(1.0-r_2**2)/(1.0-r_3**2)
     rkappa2=1.0-rk2
     rkappa=sqrt(rkappa2)
     FLegell_1 = el2(x,rkappa,1.0,1.0)
     ELegell_2 = el2(x,rkappa,1.0,rkappa2)
     qrat=r_2*r_3/sinth**3
     bvec_2(1)=qrat*(FLegell_1-ELegell_2)/rk2
     bvec_2(3)=qrat*((r_2/r_3)*sinth-ELegell_2)/rkappa2
     bvec_2(2)=1.0-(bvec_2(1)+bvec_2(3))
     !c     b_2=qrat*(ELegell_2-rkappa2*FLegell_1-(r_3/r_2)*rk2*sinth)/rk2/rkappa2
  endif
  return
end subroutine get_b_2

FUNCTION elliptic_rd(x,y,z)
  USE intreal_types
  real(sp)  elliptic_rd,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6
  PARAMETER (ERRTOL=0.05,TINY=1.0e-25,BIG=4.5e21,C1=3.0/14.0,C2=1.0/6.0,&
       C3=9.0/22.0,C4=3.0/26.0,C5=0.25*C3,C6=1.5*C4)
  REAL(SP) alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,&
       sqrtz,sum,xt,yt,zt
  !c     COMPUTES CARLSON'S ELLIPTIC INTEGRAL 
  !C     R_D(x,y,z)=1.5*int_0^infty dt/[(t+z)^3/2 (t+x)^1/2 (t+y)^1/2]
  !C     x,y >= 0, AND AT MOST ONE ZERO
  if(min(x,y).lt.0.0.or.min(x+y,z).lt.TINY.or.max(x,y,z).gt.BIG) then
     write(*,*) 'invalid arguments in rd '
  endif
  xt=x
  yt=y
  zt=z
  sum=0.0
  fac=1.0
1 continue
  sqrtx=sqrt(xt)
  sqrty=sqrt(yt)
  sqrtz=sqrt(zt)
  alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
  sum=sum+fac/(sqrtz*(zt+alamb))
  fac=0.25*fac
  xt=0.25*(xt+alamb)
  yt=0.25*(yt+alamb)
  zt=0.25*(zt+alamb)
  ave=0.2*(xt+yt+3.0*zt)
  delx=(ave-xt)/ave
  dely=(ave-yt)/ave
  delz=(ave-zt)/ave
  if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL) go to 1
  ea=delx*dely
  eb=delz*delz
  ec=ea-eb
  ed=ea-6.0*eb
  ee=ed+ec+ec
  elliptic_rd=3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
  return
END FUNCTION elliptic_rd

FUNCTION EL2(X,QQC,AA,BB)
  USE intreal_types
  PARAMETER(PI=3.14159265, CA=.0003, CB=1.E-9)
  IF(X.EQ.0.)THEN
     EL2=0.
  ELSE IF(QQC.NE.0.)THEN
     QC=QQC
     A=AA
     B=BB
     C=X**2
     D=1.+C
     P=SQRT((1.+QC**2*C)/D)
     D=X/D
     C=D/(2.*P)
     Z=A-B
     EYE=A
     A=0.5*(B+A)
     Y=ABS(1./X)
     F=0.
     L=0
     EM=1.
     QC=ABS(QC)
1    B=EYE*QC+B
     E=EM*QC
     G=E/P
     D=F*G+D
     F=C
     EYE=A
     P=G+P
     C=0.5*(D/P+C)
     G=EM
     EM=QC+EM
     A=0.5*(B/EM+A)
     Y=-E/Y+Y
     IF(Y.EQ.0.)Y=SQRT(E)*CB
     IF(ABS(G-QC).GT.CA*G)THEN
        QC=SQRT(E)*2.
        L=L+L
        IF(Y.LT.0.)L=L+1
        GO TO 1
     ENDIF
     IF(Y.LT.0.)L=L+1
     E=(ATAN(EM/Y)+PI*L)*A/EM
     IF(X.LT.0.)E=-E
     EL2=E+C*Z
  ELSE
     write(*,*) 'failure in EL2'
  ENDIF
  RETURN
END FUNCTION EL2


SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT,DERIVS)
  USE intreal_types
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
END SUBROUTINE RK4

subroutine evolve_ellipse_full(idynax,zdynax,Ddynax,fcdynax)
  USE intreal_types
  use arrays
  use evpvtable
  use params
  use bvalues
  use hard
  use evalues 
  use equil
  use univer
  use input_parameters
 
  PARAMETER (ny=6)
  PARAMETER (one_third=1.0/3.0)
  parameter (nzdynmax=7)
  dimension y(ny),dy(ny)

  dimension zdynv(nzdynmax),delta1_ev(nzdynmax),aturnv(3),dcritv(3)
  integer(i4b) ldynv(nzdynmax),ldyn_pv(nzdynmax)
  external get_derivs
  no_vir=3
  no_turn=3
  no_dens=1
  nzdyn=no_vir+no_turn+no_dens
  ivir0=1
  ivir1=no_vir
  iturn0=no_vir+1
  iturn1=no_vir+no_turn
  idens0=no_vir+no_turn+1
  idens1=no_vir+no_turn+no_dens
  do i=1,3
     aturnv(i)=0.0
  enddo
  dcritv(1)=dcrit
  !c            zvir3=zdynv(3)
  !c            zvir2=zdynv(2)
  !c            zvir1=zdynv(1)
  !c            zturn3=zdynv(6)
  !c            zturn2=zdynv(5)
  !c            zturn1=zdynv(4)
  !c            zdcrit=zdynv(7)
  do ii=1,nzdyn
     ldynv(ii)=0
     delta1_ev(ii)=-1.0
     zdynv(ii)=-1.0
  enddo

  nyp=ny
  a_3eq=0.0
  a_2eq=0.0
  a_1eq=0.0
  lvirv(1)=0
  lvirv(2)=0
  lvirv(3)=0
  zdynax=-1
  Ddynax=-1
  fcdynax=-1
  if(iwant_evmap.eq.1) then
     Frho=Fbar
     e_v=e_v+de_v
     if(e_v.gt.e_vmax) then
        stop
     endif
     p_v=p_vbar
  elseif(iwant_evmap.eq.2) then
     Frho=Fbar
     e_v=e_vbar
     p_v=p_v+dp_v
     if(p_v.gt.e_v) then
        stop
     endif
  elseif(iwant_evmap.eq.0) then
     if(e_v.lt.0.0) then
        stop
     endif
  elseif(iwant_evmap.eq.-1) then
     iwant_evmap=0
  endif


  aLam_3 = Frho/3.0*(1.0+3.0*e_v+p_v)
  aLam_2 = Frho/3.0*(1.0-2.0*p_v)
  aLam_1 = Frho/3.0*(1.0-3.0*e_v+p_v)
  if(aLam_3.lt.aLam_2.or.aLam_3.lt.aLam_1.or.aLam_2.lt.aLam_1) then
     write(*,*) 'problem with aLams, exiting...',aLam_1,aLam_2,aLam_3
     call mpi_finalize(ierr)
     stop
  endif

  !C ROUGH ESTIMATE OF COLLAPSE zzc FROM aLam_1
  zzc_180sph=(aLam_1+aLam_2+aLam_3)/1.606
  zzc_Fsph=(aLam_1+aLam_2+aLam_3)/1.686
  !c      zzc_1pan=aLam_1
  !c      zzc_3pan=aLam_3

  call ic_set(zinit,aLam_1,aLam_2,aLam_3,t,nyp,y)

  istep=0
111 continue
  istep=istep+1
  a_3=y(3)
  dtstep=tfac*a_3*sqrt(a_3)*Ha_b_nr
  !C NOTE FUNNY RK4 OF PRESS REQUIRES FIRST CALL TO derivs
  call get_derivs(t,y,dy)
  call rk4(y,dy,nyp,t,dtstep,y,get_derivs)
  t=t+dtstep
  delta1_e=t**3/(y(1)*y(2)*y(3))

  if(ivir1.ge.ivir0) then
     do ii=ivir0,ivir1
        if(lvirv(ii).ne.ldynv(ii)) then
           zdynv(ii)=1.0/t
           ldynv(ii)=1
           !c               delta1_ev(ii)=t**3/(y(1)*y(2)*y(3))
           delta1_ev(ii)=delta1_e
        endif
     enddo
  endif
  if(iturn1.ge.iturn0) then
     do ii=iturn0,iturn1
        if(ldynv(ii).ne.1) then
           iax=ii-ivir1
           if(y(iax).lt.aturnv(iax)) then
              zdynv(ii)=1.0/t
              ldynv(ii)=1
              delta1_ev(ii)=delta1_e
           else
              aturnv(iax)=y(iax)
           endif
        endif
     enddo
  endif

  if(idens1.ge.idens0) then
     do ii=idens0,idens1
        if(ldynv(ii).ne.1) then
           if(delta1_e.ge.dcritv(ii-iturn1)) then
              zdynv(ii)=1.0/t
              ldynv(ii)=1
              delta1_ev(ii)=delta1_e
           endif
        endif
     enddo
  endif
  if(iwant_evmap.eq.0) then
     if((istep.eq.(istep/nout)*nout).or.lvirv(1).eq.1.or.istep.ge.nstepmax) then
        zz=1.0/t
        a_b=1.0/zz
        delta1_e=t**3/(y(1)*y(2)*y(3))
        tau_1=(bvec_2(1)-one_third) 
        tau_2=(bvec_2(2)-one_third) 
        tau_3=(bvec_2(3)-one_third) 
11      format(i5,6(1pe11.3),2x,2i2)
115     format(' tau321=(b_i/2-1/3) ',3(1pe11.3))
     endif
  endif
  if(lvirv(1).eq.1.or.istep.ge.nstepmax.or.t.ge.2.0) then
     if(iwant_evmap.eq.0) write(*,*) 'last axis (1) is virialized '
67   format(' F,e_v,p_v,Lam123: ',6f9.3)
     !c         write(*,66) zturn3,zvir3,zturn2,zvir2,zturn1,zvir1,zdcrit
66   format('zvir123,zturn,dcr: ',7(f9.3))
68   format('zc180sph,zcF: ',2(F9.3),' dvir,dcr: ',4(F9.3))
     if(ihard.eq.1.and.iwant_evmap.ne.3) then
        write(10,67) Frho,e_v,p_v,aLam_1,aLam_2,aLam_3
        write(10,66) (zdynv(ii),ii=1,nzdyn)
        write(10,68) zzc_3sph,zzc_Fsph,delta1_ev(3),delta1_ev(2),delta1_ev(1),delta1_ev(7)
        write(10,169) (delta1_ev(ii),ii=1,nzdyn)
     elseif(ihard.eq.2.and.iwant_evmap.ne.3) then
        write(10,69) Frho,e_v,p_v,zdynv(3),zdynv(2),zdynv(1),zdynv(7),zzc_Fsph
69      format(8F9.3)
169     format('densvir,turn,cut: ',7(f9.3))
     endif
     zdynax=zdynv(idynax)
     if(zdynax.eq.-1) then
        Ddynax=-1
        fcdynax=-1
     else
        Ddynax=Dfnofa(1.0/zdynax)
        fcdynax=Frho*Ddynax
     endif
     return
  endif
  goto 111
end subroutine evolve_ellipse_full



