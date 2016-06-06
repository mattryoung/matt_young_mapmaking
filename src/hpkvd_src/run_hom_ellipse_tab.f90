PROGRAM run_hom_ellipse
  call run_hom_ellipse_tab
END PROGRAM run_hom_ellipse

!C THE PROGRAM IS SET UP FOR vac AND open MODELS
!C  NOTE THAT THE SCALING IS SUCH THAT THE 
!C  VELOCITY IS DIVIDED BY H_0 Omnr^{1/2} 
!C   INTEGRATION IS wrt a HENCE THERE IS ALSO AN adot 
!C   IN THE DENOMINATOR, SO WE HAVE 
!C    dv\tau_0/da  = \tau_0^-2 4pi G \rho_nr X STUFF = (3/2)a_b^-3 X STUFF

subroutine run_hom_ellipse_tab
  USE intreal_types
  use wantevpvtab
  use evpvtable
  use params
  use equil
  use evalues
  use hard
  use rd

  !C `TIME' INTEGRATION VARIABLE t=a_background
  !c    COSMIC TIME UNIT= (8\pi G/3 rho_bnr*)^{-1/2}

  !C  a_idot=[p_i](Ha)^-1 
  !C  p_idot=-(3/4)[delta1_e b_i + (2/3 - b_i)]a_b^-3 a_i (Ha)^-1
  !C  delta1_e = (a_b/a_1)(a_b/a_2)(a_b/a_3)
  !C  dot wrt a_b  .... 
  !C  Ha=sqrt(1+OMvac/Omnra_b^3+Omcurv/Omnr a_b) /a_b^{1/2}
  !C  b_idot=[-2p_i/a_i + SUM_{j.ne.i} [b_i-b_j]/[a_i^2-a_j^2] 
  !C     (a_ip_i-a_jp_j) + b_i SUM_j p_j/a_j](Ha)^-1
  !C INSTEAD WE USE THE BINNEY-TREMAINE FORMULA FOR b_i

  parameter(Npv_evmax=10,Npv_evmin=-10,Nevmax=19,Nfscmax=10)
  parameter (ny=6)

  dimension quicktab(0:1,0:1)

  character(LEN=120) fcvoutfile
  character(LEN=120) dummy
  external get_derivs
  external get_derivs_sph

  call read_cosmology(Omt,Omnr,Omx,OmB,h,Omcurv,Omvac,Omhdm,fhdmclus,iamcurved,dcurv)
  call Dlinear_cosmology(Omt,Omnr,Omx,OmB,h,Omcurv,Omvac,Omhdm,fhdmclus,iamcurved,dcurv)

  igethom_strat=1
  if(omcurv.ne.0.or.omvac.ne.0) igethom_strat=2
  igethom_strat=1
  !C FORCE THIS TO CHECK TABLE USE
  write(*,*) ' iwant_evpvtab /1 READ, 2 MAKE, 3 MAKE-READ/'
  read(999,*) iwant_evpvtab
  if(iwant_evpvtab.eq.3) ihard=1
  if(iwant_evpvtab.eq.2) ihard=2
  write(*,*) ' table check '
  read(999,*) itabcheck
  if(itabcheck.eq.1) then
     call inithom_ellipse_params
     write(*,*) 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx', dcrit
     if(iwant_evpvtab.eq.1) then
        call make_F_zvir_table
        call readhom_ellipse_tab
     else
        call make_F_zvir_table
        write(*,*) ' going into makehom_ellipse_tab '
        call makehom_ellipse_tab
     endif
     Zb=1.0
     call inithom_ellipse_redshift(Zb)
1915 continue
     write(*,*) ' Frho,e_v,p_v /e_v< 0 stop /'
     read(999,*) Frho,e_v,p_v 
     if(e_v.lt.0.0) stop
     call gethom_ellipse(fcvir1tab,Dvir1tab,zvir1tab,iflag,igethom_strat,iuseinterp)    
     write(*,*) 'XXXXXXXXXXX iuseinterp ',iuseinterp
     if(zvir1tab.gt.0.0) then
        aB=1.0/zvir1tab
        Dtab1=Dlinear(aB,chi,HD_Ha,D_aB)
        zvir1check=1.0/afnofD(Dvir1tab)
        !c            aB=1.0/(zvir1tab*D_a)
        !c            Dtab2=Dlinear(aB,chi,HD_Ha,D_a)
        !c            aB=1.0/(zvir1tab*D_a)
        !c            Dtab3=Dlinear(aB,chi,HD_Ha,D_a)
     endif
     idynax=1
     zinit_fac_ell=20.0
     zinit=max(zinit_fac_ell,zinit_fac_ell*Frho)
     call evolve_ellipse(idynax,zvir1f,Dvir1f,fcvir1f)
     if(zvir1f.gt.0.0) then
        aB=1.0/zvir1f
        Df=Dlinear(aB,chi,HD_Ha,D_af)
        !c            Dest=Dlinear(1.0/(zvir1tab*D_a),chi,HD_Ha,D_a)
     endif

     iwant_evmap=-1
     zinit_fac_ell=20.0
     zinit=max(zinit_fac_ell,zinit_fac_ell*Frho)
     call evolve_ellipse_full(idynax,zvir1,Dvir1,fcvir1)
     if(zvir1.gt.0.0) then
        aB=1.0/zvir1
        D1=Dlinear(aB,chi,HD_Ha,D_a1)
     endif
     write(*,*) '*** Frho,e_v,p_v ',Frho,e_v,p_v
     write(*,*) 'FFFF fcvir1tab,fcvir1f,fcvir1 ',fcvir1tab,fcvir1f,fcvir1
     write(*,*) 'ZZZZ zvir1 ',zvir1tab,zvir1check,zvir1f,zvir1
     write(*,*) 'DDDD Dtab: ',Dtab1,Dvir1tab,Dvir1f,Dvir1
     write(*,*) 'DaDaDa : ',D_aB,D_af,D_a1
     !c         write(*,*) '*** Dtab: ',Dtab1,Dtab2,Dtab3
     !c         write(*,*) '......!D: ',Dest,Df,D1
     goto 1915
  else
     write(*,*) 'dummy '
     read(999,'(a)') dummy
     write(*,*) ' Frho,e_v,p_v /DUMMY /'
     read(999,*) Frho,e_v,pv 
     write(*,*) ' Frho,e_v,p_v /DUMMY /'
     read(999,*) Frho,e_v,p_v 
  endif

  nstepmax=100000
  nyp=ny
  call read_inithom_ellipse
  call inithom_ellipse_redshift(Zb)
  if(iwant_evmap.ne.0) then
     if(iwant_evpvtab.eq.1) then
        write(*,*) 'Nfsc,ifc1inv,Fsmin,Fsmax,zinit_fac_ell,','devtab,dpv_evtab,evtabmax,pv_evtabmax'
        read(999,*) Nfsc,ifc1inv,Fsmin,Fsmax,zinit_fac_ell,devtab,dpv_evtab,evtabmax,pv_evtabmax
        if(Nfsc.gt.1) then
           fscmin=fsc_tabfn(Fsmin)*1.0001
           fscmax=fsc_tabfn(Fsmax)*0.9999
           if(fscmin.ge.fscmax) then 
              Nfsc=1 
              dfsc=0.0
           else
              dfsc=(fscmax-fscmin)/(Nfsc-1)
           endif
           do ktab=1,Nfsc
              fsctabv(ktab)=fscmin+dfsc*(ktab-1) 
              Ftabv(ktab)= Frho_of_fsc_tabfn(fsctabv(ktab))
           enddo
           write(*,*) 'fsctab  ',(fsctabv(ktab),ktab=1,Nfsc)
        else
           fc_EdeS=1.686
           fscmin=fc_EdeS
           fscmax=fc_EdeS
           fsctabv(1)=fc_EdeS
           Ftabv(1)=Fsmax 
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
     endif
     if(iwant_evmap.eq.1) then
        write(*,*) ' Fbar,e_vmax,de_v,p_vbar /p_vbar is fixed /'
        read(999,*) Fbar,e_vmax,de_v,p_vbar
        e_v=-de_v
     elseif(iwant_evmap.eq.2) then
        write(*,*) ' Fbar,e_vbar,dp_v /e_vbar is fixed /'
        read(999,*) Fbar,e_vbar,dp_v
        p_v=-e_vbar-dp_v
     endif
  else
     write(*,*) 'DUMMY Nfsc,ifc1inv,Fsmin,Fsmax,zinit_fac_ell,','devtab,dpv_evtab,', 'evtabmax,pv_evtabmax'
     read(999,*) Nfsc,ifc1inv,Fsmin,Fsmax,zinit_fac_ell,devtab,dpv_evtab,evtabmax,pv_evtabmax

  endif
  write(*,*) 'hardcopy? /1 usual output, 2 for plot or tab/'
  read(999,*) ihard
  write(*,*) ' FILENAME: /(1) zevout file OR (2) zpvout file '
  write(*,*) ' ... OR (3) tabfile/'
  write(*,*) ' ... OR DUMMY'
  read(999,'(a)') zevoutfile
  if(ihard.ne.0) then
     if(iwant_evpvtab.eq.0) then
        open(10,file=trim(fcvoutfile),form='formatted',status='unknown')
     elseif(iwant_evpvtab.eq.1) then
        open(10,file=trim(fcvoutfile),form='unformatted',status='unknown')
     endif
  endif

110 continue
  if(iwant_evpvtab.eq.1) then
     zinit_fac_ell=20.0
     idynax=1
     do ktab=1,Nfsc
        Frho=Ftabv(ktab)
        zinit=max(zinit_fac_ell,zinit_fac_ell*Frho)
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
119        continue
        enddo
     enddo
     write(10) Nfsc,Nev,Npv_ev,devtab,dpv_evtab,evtabmax,pv_evtabmax,Fsmin,Fsmax,zinit_fac_ell, & 
          (fsctabv(ktab),ktab=1,Nfsc),(Ftabv(ktab),ktab=1,Nfsc),(((fc1mat(it,jt,kt),it=-Npv_ev,Npv_ev),jt=0,Nev),kt=1,Nfsc)
     close(10)
889  continue
     write(*,*) 'input ev,pv /ev < 0 stop /'
     read(999,*) ev,pv
     if(ev.lt.0.) stop
     if(ev.ne.0.0) then
        if(abs(pv).gt.ev) then
           fc1app3=0.0
           fc1app4=0.0
        else
           ritab=(pv/ev)*Npv_ev
           itab_1=ritab
           if(ritab.lt.0.0) itab_1=itab_1-1
           itab=itab_1+1
           diffi=ritab-itab_1
           rjtab=ev/devtab
           jtab_1=rjtab
           diffj=rjtab-jtab_1
           jtab=jtab_1+1
           write(*,*) fc1mat(itab_1,jtab_1,1),fc1mat(itab,jtab_1,1),fc1mat(itab_1,jtab,1),fc1mat(itab,jtab,1)
           fc1app3=fc1mat(itab_1,jtab_1,1) +(diffi*(fc1mat(itab,jtab_1,1) -fc1mat(itab_1,jtab_1,1)) &
                +diffj*(fc1mat(itab_1,jtab,1)-fc1mat(itab_1,jtab_1,1)))
           fc1app4=fc1app3+(diffi*diffj*((fc1mat(itab,jtab,1)-fc1mat(itab_1,jtab,1)) &
                -(fc1mat(itab,jtab_1,1)-fc1mat(itab_1,jtab_1,1))))
        endif
     else
        fc1app3=fc1mat(0,0,1) 
        fc1app4=fc1app3
     endif
     write(*,*) ' fc1app3,fc1app4  ',fc1app3,fc1app4
     !c we could get rid of the diffi*diffj step without losing accuracy
     goto 889
  endif

  idynax=1
  write(*,*) 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX', dcrit
  call evolve_ellipse_full(idynax,zdynax,Ddynax,fcdynax)
  goto 110
end subroutine run_hom_ellipse_tab

subroutine evolve_ellipse_full_1(idynax,zdynax,Ddynax,fcdynax)
  USE intreal_types
  use evpvtable
  use params
  use bvalues
  use hard
  use evalues
  use equil
  use univer

  PARAMETER (ny=6)
  PARAMETER (one_third=1.0/3.0)
  parameter (nzdynmax=7)
  parameter(Npv_evmax=10,Npv_evmin=-10,Nevmax=19,Nfscmax=10)
  dimension y(ny),dy(ny)

  dimension zdynv(nzdynmax),delta1_ev(nzdynmax),aturnv(3),dcritv(3)
  integer ldynv(nzdynmax),ldyn_pv(nzdynmax)
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
     write(*,*) 'Frho,e_v,p_v /if e_v negative stop/'
     read(999,*) Frho,e_v,p_v
     !c      write(*,*) 'aLam_1,aLam_2,aLam_3 '
     !c      read(999,*) aLam_1,aLam_2,aLam_3
     if(e_v.lt.0.0) then
        stop
     endif
  elseif(iwant_evmap.eq.-1) then
     write(*,*) ' uses default Frho,e_v,p_v '
     iwant_evmap=0
  endif
  aLam_3 = Frho/3.0*(1.0+3.0*e_v+p_v)
  aLam_2 = Frho/3.0*(1.0-2.0*p_v)
  aLam_1 = Frho/3.0*(1.0-3.0*e_v+p_v)
  if(aLam_3.lt.aLam_2.or.aLam_3.lt.aLam_1.or.aLam_2.lt.aLam_1) return
  !C ROUGH ESTIMATE OF COLLAPSE zzc FROM aLam_1
  zzc_180sph=(aLam_1+aLam_2+aLam_3)/1.606
  zzc_Fsph=(aLam_1+aLam_2+aLam_3)/1.686
  !c      zzc_1pan=aLam_1
  !c      zzc_3pan=aLam_3

  call ic_set(zinit,aLam_1,aLam_2,aLam_3,t,nyp,y)

  istep=0
  if(iwant_evmap.eq.0) then
     write(*,*)'         zz        a_b        a_3        a_2   ','     a_1      1+delta vir:3,2'
  endif
111 continue
  istep=istep+1
  a_3=y(3)
  dtstep=tfac*a_3*sqrt(a_3)*Ha_b_nr
  !C NOT FUNNY RK4 OF PRESS REQUIRES FIRST CALL TO derivs
  call get_derivs(t,y,dy)
  call rk4(y,dy,nyp,t,dtstep,y,get_derivs)
  t=t+dtstep
  delta1_e=t**3/(y(1)*y(2)*y(3))

  if(ivir1.ge.ivir0) then
     do ii=ivir0,ivir1
        if(lvirv(ii).ne.ldynv(ii)) then
           zdynv(ii)=1.0/t
           ldynv(ii)=1
           !c               write(*,*) 'virializing ',ii
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
        write(*,11) istep,zz,a_b,y(3),y(2),y(1),delta1_e,lvirv(3),lvirv(2)
        tau_1=(bvec_2(1)-one_third) 
        tau_2=(bvec_2(2)-one_third) 
        tau_3=(bvec_2(3)-one_third) 
        write(*,115) tau_3,tau_2,tau_1 
11      format(i5,6(1pe11.3),2x,2i2)
115     format(' tau321=(b_i/2-1/3) ',3(1pe11.3))
     endif
  endif
  if(lvirv(1).eq.1.or.istep.ge.nstepmax.or.t.ge.2.0) then
     if(iwant_evmap.eq.0) write(*,*) 'last axis (1) is virialized '
     write(*,*) '  '
     write(*,67) Frho,e_v,p_v,aLam_1,aLam_2,aLam_3
67   format(' F,e_v,p_v,Lam123: ',6f9.3)
     !c         write(*,66) zturn3,zvir3,zturn2,zvir2,zturn1,zvir1,zdcrit
     write(*,66) (zdynv(ii),ii=1,nzdyn)
66   format('zvir123,zturn,dcr: ',7(f9.3))
     write(*,68) zzc_180sph,zzc_Fsph,delta1_ev(3),delta1_ev(2),delta1_ev(1),delta1_ev(7)
68   format('zc180sph,zcF: ',2(F9.3),' dvir,dcr: ',4(F9.3))
     write(*,169) (delta1_ev(ii),ii=1,nzdyn)
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
end subroutine evolve_ellipse_full_1



