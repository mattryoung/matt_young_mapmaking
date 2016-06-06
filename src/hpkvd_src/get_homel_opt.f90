subroutine get_homel_opt(npart,ipp,alatt,hlatt,ir2min,F,&
     etavec,Sbar,RTHL,Srb,ZZon,DDon,Fnupk,Fevpk,Fpvpk,&
     strain_mat,Sbar2)

  USE intreal_types
  use bound

  !  Radial integration of F to find Fbar=f_v crossing point
  !  closest to R^2=ir2min*alatt fiducial. Return RTHL and Srb.

  !  Use homogeneous ellipsoid model to truncate integration

  !  Input:
  !              DDon = D_on correction to F at peak
  !                      (assumed uniform over peak-patch, use to
  !                       multiply F,etave!c values as DDon)

  !  Get back shear eigenvalues F_nu,F_ev,F_pv where
  !              lam1 = F_nu/3( 1 - 3*F_ev + F_pv )
  !              lam2 = F_nu/3( 1 - 2*F_pv )
  !              lam3 = F_nu/3( 1 + 3*F_ev + F_pv )

  !  So          F_nu = lam1 + lam2 + lam3 = F
  !              F_ev = ( lam3 - lam1 )/ 2*F_nu
  !              F_pv = ( lam1 + lam3 - 2*lam2 )/ 2*F_nu

  !  These are calculated over a top-hat filtering
  !  The call is      :      call evolve_ellipse(zvir1)
  !  IF zvir=-1 IT DID NOT COLLAPSE ALONG AXIS 1 (last axis)

  
  use evalues
  use particle
  use table  
  use input_parameters
  use arrays

  implicit none
  integer(i4b) iuseinterp
  real(sp) fcvir1p,Dvir1p

  integer(i4b) nn(3),npart,ipp
  real(sp) F(n1,n2,n3)
  real(sp) etavec(n1,n2,n3,3)
  integer(i4b) ipk(3),iv(3),nshell(np)
  real(sp) Fnupk,Fevpk,Fpvpk
  real(sp) Sshell(3,np),SRshell(3,3,np)
  real(sp) Sbar(3),Ebar(3,3),strain_mat(3,3)
  real(sp) Fbar(np),rad(np)

  ! ------------------------------------------------------------
  ! Sbar2 added by M.A.A. to do local 2LPT measurement
  integer  cell(3),naverage
  real(sp) flocal
  real(sp) Sbar2(3)
  real(sp) s2local(3),slocal(3)
  real(sp) strlocal(6)
  !
  ! ------------------------------------------------------------

  integer(i4b) iblack
  integer(i4b) ir2min,n1xn2,ir2upp,nsh,m,m0,m1,ir2p,ifcrit,&
               iflag,iflagel
  integer(i4b) L,K,jp,ir2,j0,icon,i2c,nSbar
  integer(i4b) mupp,muppnew,mupp_p,mlow,mlownew,mstart,mp
  real(sp) fsc_tabfn_of_ZZ
  real(sp) RTHL,Srb,ZZon,DDon,alatt,fcrit,fcritx,zvir1,zvir1p
  real(sp) Frhoh,Frhpk,hlatt,hlatt2,hlatt_1,hlatt_2,diff,con,aww
  real(sp) pi,fourpi,anor,aRnor,wnor,wRnor,one_third
  real(sp) rupp,Fshell,dFbar,dFbarp,rad3,rad3p,rlow,u0,u02,u12
  real(sp) wt,dZvir,RTHL3,RTHL5,drad3,Fbarx,rad5,rad5p

  real(sp) Lam(3)
  real(sp) hRinteg

  RTHL=0.0
  Srb=1.0
  if (ir2min.le.0.or.hlatt.le.0.0) return

  fcrit=fsc_tabfn_of_ZZ(ZZon)

  pi=4.0*atan(1.0)
  fourpi=4.0*pi
  anor=1.0/hlatt**3
  aRnor=3.0*anor/alatt
  wnor=anor/fourpi
  wRnor=3.0*wnor/alatt

  one_third=1.0/3.0
  nn(1)=n1
  nn(2)=n2
  nn(3)=n3
  n1xn2=nn(1)*nn(2)
  call get_ijk(n1xn2,nn(1),ipp,ipk(1),ipk(2),ipk(3))

  hlatt2=hlatt*hlatt
  hlatt_1=1.0/hlatt
  hlatt_2=hlatt_1*hlatt_1

  rupp=sqrt(float(ir2min))+2.0*hlatt
  ir2upp=int(rupp*rupp)+1

  Fshell=0.0
  nsh=0
  !c  First point = first shell
  m=1
  m0=1
  rad(1)=0.0
  Fbar(1)=F(ipk(1),ipk(2),ipk(3))
  nshell(1)=1
  ir2p=irs2(2)
  dFbarp=Fbar(1)
  do L=1,3
     Sshell(L,1)=etavec(ipk(1),ipk(2),ipk(3),L)
     Sshell(L,2)=0.0
     do K=1,3
        SRshell(L,K,1)=0.0
        SRshell(L,K,2)=0.0
     enddo
  enddo
  !c  Fill array out to ir2min (need to go out to ir2max for Ebar)
  ifcrit=1
  do jp=2,npart
     ir2=irs2(jp)
     if(ir2.ne.ir2p) then
        m=m+1
        rad(m)=sqrt(float(ir2p))
        nshell(m)=nsh
        dFbar=0.0
        if(nsh.gt.0) dFbar=Fshell/float(nsh)
        rad3p=rad(m-1)**3
        rad3=rad(m)**3
        Fbar(m)=(rad3p*Fbar(m-1)+0.5*(dFbarp+dFbar)*(rad3-rad3p))/rad3
        dFbarp=dFbar

        if(ifcrit.eq.1) then
           m0=m
           if(ir2p.gt.ir2min) then
              rupp=rad(m0)+2.0*hlatt
              ir2upp=int(rupp*rupp)+1
              ifcrit=0
              mupp=m
           endif
        else
           if(rad(m).lt.rupp) mupp=m
        endif

        ir2p=ir2
        Fshell=0.0
        nsh=0
        do L=1,3
           Sshell(L,m+1)=0.0
           do K=1,3
              SRshell(L,K,m)=SRshell(L,K,m)/rad(m)
              SRshell(L,K,m+1)=0.0
           enddo
        enddo

     endif
     iblack=0
     do L=1,3
        iv(L)=ipk(L)+ixsvec(jp,L)
        !c  Implicit wrap.
        if(iwrap.eq.1) then
           if(iv(L).lt.1) iv(L)=iv(L)+nn(L)
           if(iv(L).gt.nn(L)) iv(L)=iv(L)-nn(L)
        else
           if((iv(L).lt.1).or.(iv(L).gt.nn(L))) iblack=1
        endif
     enddo
     if(iblack.eq.0) then
        nsh=nsh+1
        Fshell=Fshell+F(iv(1),iv(2),iv(3))
        do L=1,3
           Sshell(L,m+1)=Sshell(L,m+1)+etavec(iv(1),iv(2),iv(3),L)
           do K=1,3
              SRshell(L,K,m+1)=SRshell(L,K,m+1)+etavec(iv(1),iv(2),iv(3),L)*ixsvec(jp,K)
           enddo
        enddo
     endif
     !c  Now bounce out if at ir2max
     if (ir2.gt.ir2upp) goto 252
  enddo
  !c  Check at ir2min (m0) for Fbar=fcrit
252 continue

  if(Fbar(m0).ge.fcrit) then
     !c  Go outward to Fbar<fcrit
     ifcrit=1
     j0=jp+1
     do jp=j0,npart
        ir2=irs2(jp)
        if(ir2.ne.ir2p) then
           m=m+1
           rad(m)=sqrt(float(ir2p))
           nshell(m)=nsh
           dFbar=0.0
           if(nsh.gt.0) dFbar=Fshell/float(nsh)
           rad3p=rad(m-1)**3
           rad3=rad(m)**3
           Fbar(m)=(rad3p*Fbar(m-1)+0.5*(dFbarp+dFbar)*(rad3-rad3p))/rad3
           dFbarp=dFbar

           !c  check for Fbar<fcrit : make sure go out to r+2h
           if(ifcrit.eq.1) then
              m0=m
              rupp=rad(m0)+2.0*hlatt
              if(Fbar(m).lt.fcrit) ifcrit=0
           else
              if(rad(m).ge.rupp) then
                 mupp=m-1
                 goto 254
              endif
           endif

           ir2p=ir2
           Fshell=0.0
           nsh=0
           do L=1,3
              Sshell(L,m+1)=0.0
              do K=1,3
                 SRshell(L,K,m)=SRshell(L,K,m)/rad(m)
                 SRshell(L,K,m+1)=0.0
              enddo
           enddo
        endif
        iblack=0
        do L=1,3
           iv(L)=ipk(L)+ixsvec(jp,L)
           if(iwrap.eq.1) then
              if(iv(L).lt.1) iv(L)=iv(L)+nn(L)
              if(iv(L).gt.nn(L)) iv(L)=iv(L)-nn(L)
           else
              if((iv(L).lt.1).or.(iv(L).gt.nn(L))) iblack=1
           endif
        enddo
        if(iblack.eq.0) then
           nsh=nsh+1
           Fshell=Fshell+F(iv(1),iv(2),iv(3))
           do L=1,3
              Sshell(L,m+1)=Sshell(L,m+1)+etavec(iv(1),iv(2),iv(3),L)
              do K=1,3
                 SRshell(L,K,m+1)=SRshell(L,K,m+1) + ixsvec(jp,K)*etavec(iv(1),iv(2),iv(3),L)
              enddo
           enddo
        endif
     enddo

     !c  To reach here, end of array was encountered without fcrit
     
     if(ifcrit.eq.1) then
        mupp=m     
     else
        return
     endif
254  continue

  else
     !c  Back down through shells to find fcrit

     mstart=max(m0-1,1)
     do mp=mstart,1,-1
        if(Fbar(mp).ge.fcrit) goto 256
        m0=mp
     enddo

     !c  To reach here, first shell was reached without downcross

     RTHL=-1.0
     Srb=0.0
     return

256  continue

     !c  Find mupp
     rupp=rad(m0)+2.0*hlatt
     muppnew=m0
     do m1=m0+1,mupp
        if(rad(m1).lt.rupp) muppnew=m1
     enddo
     mupp=muppnew

  endif

  !c  Now have Fbar(m0-1).ge.fcrit.gt.Fbar(m0)
  !c  Check inward for homeoellipse virialization
  !c  Find mlow
  rlow=max(rad(m0)-2.0*hlatt,0.0)
  mlow=m0
  if(m0.gt.1) then
     do m1=m0-1,1,-1
        if(rad(m1).gt.rlow) mlow=m1
     enddo
  endif
  mupp_p=mupp

  !c  Zvir at m0
  !c  Need to sum over particles on shells mlow to mupp
  do L=1,3
     do K=1,3
        Ebar(L,K)=0.0
     enddo
  enddo
  do m1=mlow,mupp
     u0=hlatt_1*(rad(m0)-rad(m1))
     u02=u0*u0
     if(u02.lt.4.0) then
        u12=hlatt_2*(rad(m0)**2+rad(m1)**2)
        wt=hRinteg(u0,u12)/(u12-u02)**2
        do L=1,3
           do K=1,3
              Ebar(L,K)=Ebar(L,K)+0.5*wt*(SRshell(L,K,m1)+SRshell(K,L,m1))
           enddo
        enddo
     endif
  enddo

  do L=1,3
     do K=1,3
        Ebar(L,K)=wRnor*Ebar(L,K)/rad(m0)
     enddo
  enddo

  !c  Solve for eigenvalues
  call get_evals(3,Ebar,Lam,iflag)
  !c  These should be in INCREASING order. 
  Frho = Lam(1) + Lam(2) + Lam(3)
  if(Frho.gt.0.0) then
     e_v = 0.5*(Lam(3) - Lam(1))/Frho
     p_v = 0.5*(Lam(3) + Lam(1) - 2.0*Lam(2))/Frho
  else
     e_v = 0.0
     p_v = 0.0
  endif
  !c  Renormalize to Fbar
  Frhoh= Frho
  Frho = Fbar(m0)
  do L=1,3
     do K=1,3
        strain_mat(L,K)=Ebar(L,K)
     enddo
  enddo
  if(iflag.eq.0) then
     call gethom_ellipse(fcvir1p,Dvir1p,zvir1p,iflagel,iuseinterp)
  else
     zvir1p=-1
  endif
  !c  Check to see that things havent virialized yet
  if(zvir1p.ge.ZZon) then
     zvir1=-1.0
     goto 300
  endif
  zvir1=zvir1p
  Frhpk=Frhoh
  Fnupk=Frho
  Fevpk=e_v
  Fpvpk=p_v

  if(m0.eq.1) goto 299
  !c  Step down through shells
  mstart=m0-1
  do mp=mstart,1,-1
     rupp=rad(mp)+2.0*hlatt
     rlow=max(rad(mp)-2.0*hlatt,0.0)
     if(rad(mupp).gt.rupp) then
        !c  New mupp
        muppnew=mp
        do m1=mp+1,mupp
           if(rad(m1).lt.rupp) muppnew=m1
        enddo
        mupp=muppnew
     endif

     if(mlow.gt.1) then
        if(rad(mlow-1).gt.rlow) then
           !c  New mlow
           mlownew=mlow
           do m1=mlow-1,1,-1
              if(rad(m1).gt.rlow) mlownew=m1
           enddo
           mlow=mlownew
        endif
     endif
     !c  Check ellipsoid
     do L=1,3
        do K=1,3
           Ebar(L,K)=0.0
        enddo
     enddo

     if(mp.gt.1) then
        do m1=mlow,mupp
           u0=hlatt_1*(rad(mp)-rad(m1))
           u02=u0*u0
           if(u02.lt.4.0) then
              u12=hlatt_2*(rad(mp)**2+rad(m1)**2)
              wt=hRinteg(u0,u12)/(u12-u02)**2
              do L=1,3
                 do K=1,3
                    Ebar(L,K)=Ebar(L,K)+0.5*wt*(SRshell(L,K,m1)+SRshell(K,L,m1))
                 enddo
              enddo
           endif
        enddo

        do L=1,3
           do K=1,3
              Ebar(L,K)=wRnor*Ebar(L,K)/rad(mp)
           enddo
        enddo

     else
        do m1=mlow,mupp
           u0=hlatt_1*rad(m1)
           con=100.0*u0*u0
           icon=int(con)
           diff=con-icon
           i2c=icon+1
           !c  linear interp
           aww=akk(i2c) + diff*(akk(i2c+1) - akk(i2c))
           wt=0.5*aww/rad(m1)
           do L=1,3
              do K=1,3
                 Ebar(L,K)=Ebar(L,K)+0.5*wt*(SRshell(L,K,m1)+SRshell(K,L,m1))
              enddo
           enddo
        enddo

        do L=1,3
           do K=1,3
              Ebar(L,K)=aRnor*Ebar(L,K)
           enddo
        enddo
     endif

     do L=1,3
        do K=1,3
           strain_mat(L,K)=Ebar(L,K)
        enddo
     enddo
     call get_evals(3,Ebar,Lam,iflag)
     Frho = Lam(1) + Lam(2) + Lam(3)
     if(Frho.gt.0.0) then
        e_v = 0.5*(Lam(3) - Lam(1))/Frho
        p_v = 0.5*(Lam(3) + Lam(1) - 2.0*Lam(2))/Frho
     else
        e_v = 0.0
        p_v = 0.0
     endif

     Frhoh= Frho
     Frho = Fbar(mp)
     if(iflag.eq.0) then
        call gethom_ellipse(fcvir1p,Dvir1p,zvir1p,iflagel,iuseinterp)
     else
        zvir1p=-1
     endif
     m0=mp+1
     if(zvir1p.ge.ZZon) goto 300

     mupp_p=mupp
     zvir1=zvir1p
     Frhpk=Frhoh
     Fnupk=Frho
     Fevpk=e_v
     Fpvpk=p_v
  enddo

299 continue

  RTHL=-1.0
  Srb=0.0
  return

300 continue

  !c  Now have localized - zvir1p.ge.ZZon.gt.zvir1 for m0-1,m0
  dZvir=zvir1p-zvir1
  if(zvir1.gt.0.0.and.dZvir.ne.0.0) then
     RTHL=rad(m0-1)+(rad(m0)-rad(m0-1))*(zvir1p-ZZon)/dZvir
  else
     RTHL=rad(m0-1)
  endif
  if(RTHL.le.0.0) then
     RTHL=-1.0
     Srb=0.0
     return
  endif
  RTHL3=RTHL**3
  RTHL5=RTHL3*RTHL*RTHL

  rad3p=rad(m0-1)**3
  rad3=rad(m0)**3
  dFbar=Fbar(m0)-Fbar(m0-1)
  drad3=rad3-rad3p
  if(zvir1.gt.0.0.and.drad3.ne.0.0) then
     Fbarx=Fbar(m0-1)+(RTHL3-rad3p)*dFbar/drad3
     Frhpk=Frhoh+(RTHL3-rad3p)*(Frhpk-Frhoh)/drad3
     !c         Fevpk=e_v+(RTHL3-rad3p)*(Fevpk-e_v)/drad3
     !c         Fpvpk=p_v+(RTHL3-rad3p)*(Fpvpk-p_v)/drad3
  else
     Fbarx=Fbar(m0-1)
  endif
  do K=1,3
     do L=1,3
        strain_mat(L,K)=strain_mat(L,K)*Fbarx/Frhoh
     enddo
  enddo
  Fnupk=Fbarx
  Fevpk=e_v
  Fpvpk=p_v
  !c volume average S (inside RTHL)
  do L=1,3
     Sbar(L)=0.0
  enddo
  nSbar=0
  do m1=1,m0-1
     if(nshell(m1).gt.0) then
        do L=1,3
           Sbar(L)=Sbar(L)+Sshell(L,m1)
        enddo
        nSbar=nSbar+nshell(m1)
     endif
  enddo
  do L=1,3
     Sbar(L)=Sbar(L)/float(nSbar)
  enddo

  !c average S on shell RTHL using SPH kernel
  !c RTHL between m-1,m so integrate between mlow,mupp_p is safe
  !c      do L=1,3
  !c         Sbarh(L)=0.0
  !c      enddo
  !c      do m1=mlow,mupp_p
  !c         u0=hlatt_1*(RTHL-rad(m1))
  !c         u02=u0*u0
  !c         if(u02.lt.4.0) then
  !c            wt=0.5*wnor*hinteg(u0)*hlatt2/(RTHL*rad(m1))
  !c            do L=1,3
  !c               Sbarh(L)=Sbarh(L)+wt*Sshell(L,m1)
  !c            enddo
  !c         endif
  !c      enddo

  !c  Integrate Fbar to get energy factor     
  Srb=0.0
  rad5p=0.0
  if(m0.gt.2) then
     do mp=2,m0-1
        rad5=rad(mp)**5
        Srb=Srb+0.5*(Fbar(mp-1)+Fbar(mp))*(rad5-rad5p)
        rad5p=rad5
     enddo
  endif
  if(zvir1.gt.0.0.and.dZvir.ne.0.0) Srb=Srb+0.5*(Fbar(m0-1)+Fbarx)*(RTHL5-rad5p)

  Srb=Srb/(Fbarx*RTHL5)

  ! -------------------------------------------------------
  ! Now make local measurement of Sbar2 by looping over all
  ! cells such that r<RTHL
  ! 
  naverage = 0
  sbar2=0
  do jp=2,npart
     if(sqrt(float(irs2(jp)))>RTHL) goto 88

     cell(:) = ipk(:)+ixsvec(jp,:) 

     slocal = etavec(cell(1),cell(2),cell(3),:)

     call get_strain(cell(1),cell(2),cell(3),n1,n2,n3,&
          etavec,strlocal,8)
     strlocal = strlocal / alatt

     ! Normalize strain
     flocal = F(cell(1),cell(2),cell(3)) / &
          (strlocal(1)+strlocal(2)+strlocal(3))
     strlocal(:) = - strlocal(:) * flocal 
     flocal = F(cell(1),cell(2),cell(3)) 

     ! Get local 2LPT 

     s2local(1) = flocal*slocal(1) + &
          (strlocal(1)*slocal(1)+strlocal(5)*slocal(3)+&
          strlocal(6)*slocal(2))
     s2local(2) = flocal*slocal(2) + &
          (strlocal(2)*slocal(2)+strlocal(4)*slocal(3)+&
          strlocal(6)*slocal(1))
     s2local(3) = flocal*slocal(3) + &
          (strlocal(3)*slocal(3)+strlocal(4)*slocal(2)+&
          strlocal(5)*slocal(1))

     s2local  = s2local  * 3. / 14.
     naverage = naverage + 1
     sbar2    = sbar2    + s2local

  enddo

88 continue

  sbar2 = sbar2 / naverage

  !
  ! Added by M.A.A. 19/06/14
  ! -------------------------------------------------------

  return
end subroutine get_homel
