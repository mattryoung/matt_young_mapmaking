module peakvoidsubs
end module peakvoidsubs

subroutine get_pks(f_v,Non,xbx,ybx,zbx,alattv,ired, Rsmooth)
  use ccoordstrat
  use input_parameters
  use mpivars
  use arrays

  real(sp) alattv(3)

  !c Non is the old value coming in, new one comes out

  cen1=0.5*(n1+1)
  cen2=0.5*(n2+1)
  cen3=0.5*(n3+1)

  xmin = xbx - dcore_box / 2
  xmax = xbx + dcore_box / 2
  ymin = ybx - dcore_box / 2
  ymax = ybx + dcore_box / 2
  zmin = zbx - dcore_box / 2
  zmax = zbx + dcore_box / 2
  do k=2,n3-1
     zc=zbx+alattv(3)*(k-cen3)
     if(zc<zmin.or.zc>zmax) cycle
     do j=2,n2-1
        yc=ybx+alattv(2)*(j-cen2)
        if(yc<ymin.or.yc>ymax) cycle
        do i=2,n1-1
           xc=xbx+alattv(1)*(i-cen1)
           if(xc<xmin.or.xc>xmax) cycle

           if(mask(i,j,k)==1) cycle
           ff=delta(i,j,k)
           if(ff.lt.f_v) cycle

	   do kk=-1,1
	     kkk=k+kk
	     do jj=-1,1
	       jjj=j+jj
	       do ii=-1,1
	         iii=i+ii

	         if(ff<delta(iii,jjj,kkk)) goto 99

	       enddo
	     enddo
	   enddo

           rc = sqrt(xc**2+yc**2+zc**2)
           ac = afn(rc)
           rdsc = 1 / ac - 1

           if(rdsc > maximum_redshift.and.ievol==1) cycle

           Non=Non+1
           xpk(Non,ired)=xc
           ypk(Non,ired)=yc
           zpk(Non,ired)=zc
           Fcollv(Non,ired)=ff
           vzpk(Non,ired)=Rsmooth
           vypk(Non,ired)=etay(i,j,k)*fac
           vxpk(Non,ired)=etax(i,j,k)*fac
           lagrange(Non,ired)=i+(j-1)*n1+(k-1)*n1*n2
           mask(i,j,k)=1
           if (Non.ge.Npkmaxl) return

99         continue
        enddo
     enddo
  enddo

  return

end subroutine get_pks

subroutine get_homel(npart,ipp,alatt,ir2min, &
     Sbar,RTHL,Srb,ZZon,Fnupk,Fevpk,Fpvpk,strain_mat,Sbar2,Rfclvi)
  USE intreal_types
  use memory_management
  use bound
  use mpivars

  !c  Radial integration of F profile to find Fbar=f_v crossing point
  !c  closest to R^2=ir2min*alatt fiducial. Return RTHL and Srb.

  !c  Use homogeneous ellipsoid model to truncate integration

  !c  Get back shear eigenvalues F_nu,F_ev,F_pv where
  !c              lam1 = F_nu/3( 1 - 3*F_ev + F_pv )
  !c              lam2 = F_nu/3( 1 - 2*F_pv )
  !c              lam3 = F_nu/3( 1 + 3*F_ev + F_pv )

  !c  So          F_nu = lam1 + lam2 + lam3 = F
  !c              F_ev = ( lam3 - lam1 )/ 2*F_nu
  !c              F_pv = ( lam1 + lam3 - 2*lam2 )/ 2*F_nu

  !c  These are calculated over a top-hat filtering
  !c  IF zvir=-1 IT DID NOT COLLAPSE ALONG AXIS 1 (last axis)

  use evalues
  use particle
  use table  
  use input_parameters
  use arrays
  use timing_diagnostics

  implicit none
  integer(i4b) iuseinterp
  real(sp) fcvir1p,Dvir1p

  integer(i4b) nn(3),npart,ipp
  real(sp) Fnupk,Fevpk,Fpvpk
  real(sp) etavec(3)
  integer ipk(3),iv(3)
  integer, allocatable :: nshell(:)
  real,    allocatable :: Sshell(:,:),SRshell(:,:,:),Fbar(:),rad(:)
  real(sp) Sbar(3),Ebar(3,3),strain_mat(3,3)

  ! ------------------------------------------------------------
  ! Sbar2 added by M.A.A. to do local 2LPT measurement
  integer  cell(3),naverage
  real(sp) flocal
  real(sp) Sbar2(3)
  real(sp) s2local(3),slocal(3)
  real(sp) strlocal(6)
  real time1,time2,time3,time4,time5
  !
  ! ------------------------------------------------------------

  integer(i4b) iblack
  integer(i4b) ir2min,n1xn2,ir2upp,nsh,m,m0,m1,ir2p,ifcrit,&
               iflag,iflagel
  integer(i4b) L,K,jp,ir2,j0,icon,i2c,nSbar
  integer(i4b) mupp,muppnew,mupp_p,mlow,mlownew,mstart,mp
  real(sp) fsc_tabfn_of_ZZ
  real(sp) RTHL,Srb,ZZon,alatt,fcrit,fcritx,zvir1,zvir1p
  real(sp) Frhoh,Frhpk,hlatt,hlatt2,hlatt_1,hlatt_2,diff,con,aww
  real(sp) pi,fourpi,anor,aRnor,wnor,wRnor,one_third
  real(sp) rupp,Fshell,dFbar,dFbarp,rad3,rad3p,rlow,u0,u02,u12
  real(sp) wt,dZvir,RTHL3,RTHL5,drad3,Fbarx,rad5,rad5p

  !c      INTEGER NOW(3),LAPSED_TIME
  !C FOR USE WITH IRIS TO GET TIMING OF VARIOUS OPERATIONS

  real(sp) Lam(3)

  !c Functions
  real(sp) hRinteg

  !G.S filterscale peak found at
  real(sp) Rfclvi

  allocate(nshell(npart+1))
  allocate(Sshell(3,npart+1))
  allocate(SRshell(3,3,npart+1))
  allocate(Fbar(npart+1))
  allocate(rad(npart+1))

  if(report_memory) then
     do m0=1,npart
        nshell(m0)=m0*ipp
        Fbar(m0)=m0*ipp**2
        rad(m0)=m0*ipp**3
        Sshell(1,m0)=m0*ipp**4
        Sshell(2,m0)=m0*ipp**5
        Sshell(3,m0)=m0*ipp**6
        do m=1,3
           SRshell(1,m0,m)=m0*ipp**7+m*npart**2
           SRshell(2,m0,m)=m0*ipp**8+m*npart**3
           SRshell(3,m0,m)=m0*ipp**9+m*npart**4
        enddo
     enddo
     call ReportMemory(C_CHAR_"High Water Mark"//C_NULL_CHAR)
     nshell=0
     Sshell=0
     SRshell=0
     Fbar=0
     rad=0
  endif
        

  RTHL=0.0
  Srb=1.0
  hlatt=1.0

  fcrit=1.686
!  fcrit=fsc_tabfn_of_ZZ(ZZon)

  !c In-box correction factor D_on
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

  !G.S 22/09/2015
  !if masked by peak larger than max search radius for filterscale then quit get_homel
  if (rmax2rs > 0.0) then
     if ( mask(ipk(1),ipk(2),ipk(3)) > int((Rfclvi/alatt)*rmax2rs) ) then
        RTHL=-1.0
        Srb=0.0
        deallocate(nshell,Sshell,SRshell,Fbar,rad)
        return
     endif
  endif
  !if masked by peak larger than the minimum search radius but less than
  !max search radius set minimum search radius to masked value
!  if ( int((0.8*mask(ipk(1),ipk(2),ipk(3)))**2) > ir2min) then
!     ir2min = int((0.8*mask(ipk(1),ipk(2),ipk(3)))**2)
!  endif

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
  call checkiv(ipk,nn,0)
  Fbar(1)=delta(ipk(1),ipk(2),ipk(3))
  nshell(1)=1
  ir2p=irs2(2)
  dFbarp=Fbar(1)
  etavec(1)=etax(ipk(1),ipk(2),ipk(3))
  etavec(2)=etay(ipk(1),ipk(2),ipk(3))
  etavec(3)=etaz(ipk(1),ipk(2),ipk(3))
  do L=1,3
     Sshell(L,1)=etavec(L)
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
        call checkm(m,npart,0)
        call checkm(m-1,npart,1)
        rad(m)=sqrt(float(ir2p))
        nshell(m)=nsh
        dFbar=0.0
        if(nsh.gt.0) dFbar=Fshell/float(nsh)
        rad3p=rad(m-1)**3
        rad3=rad(m)**3
        Fbar(m)=(rad3p*Fbar(m-1)+&
             0.5*(dFbarp+dFbar)*(rad3-rad3p))/rad3
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
        call checkm(m,npart,2)
        call checkm(m+1,npart,3)
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
        call checkiv(iv,nn,1)
        call checkm(m+1,npart,4)
        Fshell=Fshell+delta(iv(1),iv(2),iv(3))
  	etavec(1)=etax(iv(1),iv(2),iv(3))
        etavec(2)=etay(iv(1),iv(2),iv(3))
        etavec(3)=etaz(iv(1),iv(2),iv(3))
        do L=1,3
           Sshell(L,m+1)=Sshell(L,m+1)+etavec(L)
           do K=1,3
              SRshell(L,K,m+1)=SRshell(L,K,m+1)+&
                   etavec(L)*ixsvec(jp,K)
           enddo
        enddo
     endif
     !c  Now bounce out if at ir2upp
     if (ir2.gt.ir2upp) goto 252
  enddo
  !c  Check at ir2min (m0) for Fbar=fcrit
  call checkm(mupp,npart,45)
252 continue

  if(Fbar(m0).ge.fcrit) then
     !c  Go outward to Fbar<fcrit
     ifcrit=1
     j0=jp+1
     do jp=j0,npart
        ir2=irs2(jp)
        if(ir2.ne.ir2p) then
           m=m+1
           call checkm(m,npart,6)
           call checkm(m-1,npart,7)
           rad(m)=sqrt(float(ir2p))
           nshell(m)=nsh
           dFbar=0.0
           if(nsh.gt.0) dFbar=Fshell/float(nsh)
           rad3p=rad(m-1)**3
           rad3=rad(m)**3
           Fbar(m)=(rad3p*Fbar(m-1)+&
                0.5*(dFbarp+dFbar)*(rad3-rad3p))/rad3
           dFbarp=dFbar

           !c  check for Fbar<fcrit : make sure go out to r+2h
           if(ifcrit.eq.1) then
              m0=m
              call checkm(m0,npart,8)
              rupp=rad(m0)+2.0*hlatt
              call checkm(m,npart,9)
              if(Fbar(m).lt.fcrit) ifcrit=0
           else
              call checkm(m,npart,10)
              if(rad(m).ge.rupp) then
                 mupp=m-1
                 goto 254
              endif
           endif

           ir2p=ir2
           Fshell=0.0
           nsh=0
           call checkm(m,npart,11)
           call checkm(m+1,npart,12)
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
           call checkiv(iv,nn,2)
           Fshell=Fshell+delta(iv(1),iv(2),iv(3))
  	   etavec(1)=etax(iv(1),iv(2),iv(3))
  	   etavec(2)=etay(iv(1),iv(2),iv(3))
  	   etavec(3)=etaz(iv(1),iv(2),iv(3))
           call checkm(m+1,npart,13)
           do L=1,3
              Sshell(L,m+1)=Sshell(L,m+1)+etavec(L)
              do K=1,3
                 SRshell(L,K,m+1)=SRshell(L,K,m+1)+&
                      ixsvec(jp,K)*etavec(L)
              enddo
           enddo
        endif
     enddo

     !c  To reach here, end of array was encountered without fcrit
     
     if(ifcrit.eq.1) then
        mupp=m     
     else
        deallocate(nshell,Sshell,SRshell,Fbar,rad)
        return
     endif
254  continue

  else
     !c  Back down through shells to find fcrit

     mstart=max(m0-1,1)
     do mp=mstart,1,-1
        call checkm(mp,npart,14)
        if(Fbar(mp).ge.fcrit) goto 256
        m0=mp
     enddo

     !c  To reach here, first shell was reached without downcross

     RTHL=-1.0
     Srb=0.0
     deallocate(nshell,Sshell,SRshell,Fbar,rad)
     return

256  continue

     !c  Find mupp
     call checkm(m0,npart,15)
     rupp=rad(m0)+2.0*hlatt
     muppnew=m0
     do m1=m0+1,mupp
        call checkm(mupp,npart,16)
        if(rad(m1).lt.rupp) muppnew=m1
     enddo
     mupp=muppnew

  endif
  !c  Now have Fbar(m0-1).ge.fcrit.gt.Fbar(m0)
  !c  Check inward for homeoellipse virialization
  !c  Find mlow
  call checkm(m0,npart,17)
  rlow=max(rad(m0)-2.0*hlatt,0.0)
  mlow=m0
  if(m0.gt.1) then
     do m1=m0-1,1,-1
        call checkm(m1,npart,18)
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
     call checkm(m0,npart,19)
     call checkm(m1,npart,20)
     u0=hlatt_1*(rad(m0)-rad(m1))
     u02=u0*u0
     if(u02.lt.4.0) then
        u12=hlatt_2*(rad(m0)**2+rad(m1)**2)
        wt=hRinteg(u0,u12)/(u12-u02)**2
        do L=1,3
           do K=1,3
              Ebar(L,K)=Ebar(L,K)+&
                   0.5*wt*(SRshell(L,K,m1)+SRshell(K,L,m1))
           enddo
        enddo
     endif
  enddo


  call checkm(m0,npart,21)
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
  call checkm(m0,npart,22)
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
     call checkm(mp,npart,22)
     rupp=rad(mp)+2.0*hlatt
     rlow=max(rad(mp)-2.0*hlatt,0.0)
     call checkm(mupp,npart,23)
     if(rad(mupp).gt.rupp) then
        !c  New mupp
        muppnew=mp
        do m1=mp+1,mupp
           call checkm(m1,npart,24)
           if(rad(m1).lt.rupp) muppnew=m1
        enddo
        mupp=muppnew
     endif

     if(mlow.gt.1) then
        call checkm(mlow-1,npart,25)
        if(rad(mlow-1).gt.rlow) then
           !c  New mlow
           mlownew=mlow
           do m1=mlow-1,1,-1
              call checkm(m1,npart,26)
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
           call checkm(m1,npart,27)
           call checkm(mp,npart,28)
           u0=hlatt_1*(rad(mp)-rad(m1))
           u02=u0*u0
           if(u02.lt.4.0) then
              u12=hlatt_2*(rad(mp)**2+rad(m1)**2)
              wt=hRinteg(u0,u12)/(u12-u02)**2
              do L=1,3
                 do K=1,3
                    Ebar(L,K)=Ebar(L,K)+&
                         0.5*wt*(SRshell(L,K,m1)+SRshell(K,L,m1))
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
           call checkm(m1,npart,1)
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
                 Ebar(L,K)=Ebar(L,K)+&
                      0.5*wt*(SRshell(L,K,m1)+SRshell(K,L,m1))
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
     call checkm(mp,npart,0)
     Frho = Fbar(mp)
     if(iflag.eq.0) then
        call gethom_ellipse(fcvir1p,Dvir1p,zvir1p,&
             iflagel,iuseinterp)
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
  deallocate(nshell,Sshell,SRshell,Fbar,rad)
  return

300 continue

  !c  Now have localized - zvir1p.ge.ZZon.gt.zvir1 for m0-1,m0
  dZvir=zvir1p-zvir1
  if(zvir1.gt.0.0.and.dZvir.ne.0.0) then
     call checkm(m0,npart,1)
     call checkm(m0-1,npart,1)
     RTHL=rad(m0-1)+(rad(m0)-rad(m0-1))*(zvir1p-ZZon)/dZvir
  else
     call checkm(m0-1,npart,1)
     RTHL=rad(m0-1)
  endif
  if(RTHL.le.0.0) then
     RTHL=-1.0
     Srb=0.0
     deallocate(nshell,Sshell,SRshell,Fbar,rad)
     return
  endif
  RTHL3=RTHL**3
  RTHL5=RTHL3*RTHL*RTHL

  call checkm(m0,npart,1)
  call checkm(m0-1,npart,1)
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
     call checkm(m1,npart,3)
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

  !c  Integrate Fbar to get energy factor     
  Srb=0.0
  rad5p=0.0
  if(m0.gt.2) then
     do mp=2,m0-1
        call checkm(mp,npart,0)
        rad5=rad(mp)**5
        call checkm(mp-1,npart,0)
        Srb=Srb+0.5*(Fbar(mp-1)+Fbar(mp))*(rad5-rad5p)
        rad5p=rad5
     enddo
  endif
  call checkm(m0-1,npart,0)
  if(zvir1.gt.0.0.and.dZvir.ne.0.0) &
       Srb=Srb+0.5*(Fbar(m0-1)+Fbarx)*(RTHL5-rad5p)

  Srb=Srb/(Fbarx*RTHL5)

  ! ----------------------------------------------------------
  ! Now make local measurement of Sbar2 by looping over all 
  ! cells such that r<RTHL
  ! 
  naverage = 0
  sbar2=0
  do jp=2,npart
     if(sqrt(float(irs2(jp)))>RTHL) cycle
     do L=1,3
        iv(L)=ipk(L)+ixsvec(jp,L)
        !c  Implicit wrap.
        if(iwrap.eq.1) then
           if(iv(L).lt.1) iv(L)=iv(L)+nn(L)
           if(iv(L).gt.nn(L)) iv(L)=iv(L)-nn(L)
        else
           if(iv(L)<1.or.iv(L)>nn(L)) goto 88
        endif
     enddo

     slocal(1) = etax(iv(1),iv(2),iv(3))
     slocal(2) = etay(iv(1),iv(2),iv(3))
     slocal(3) = etaz(iv(1),iv(2),iv(3))

     call get_strain(iv(1),iv(2),iv(3),n1,n2,n3,&
          etax,etay,etaz,strlocal,8)
     strlocal = strlocal / alatt

     ! Normalize strain
     call checkiv(iv,nn,3)
     flocal = delta(iv(1),iv(2),iv(3)) / &
          (strlocal(1)+strlocal(2)+strlocal(3))
     strlocal(:) = - strlocal(:) * flocal 
     flocal = delta(iv(1),iv(2),iv(3)) 

     ! Get local 2LPT 

     s2local(1) = flocal*slocal(1) + &
      (strlocal(1)*slocal(1)+strlocal(5)*slocal(3)+strlocal(6)*slocal(2))
     s2local(2) = flocal*slocal(2) + &
      (strlocal(2)*slocal(2)+strlocal(4)*slocal(3)+strlocal(6)*slocal(1))
     s2local(3) = flocal*slocal(3) + &
      (strlocal(3)*slocal(3)+strlocal(4)*slocal(2)+strlocal(5)*slocal(1))

     s2local = s2local * 3. / 14.

     naverage = naverage + 1

     sbar2 = sbar2 + s2local
   
     ! G.S. 22.09.15 UPDATE MASK OF ALL VALUES WITHIN RADIUS RTHL
     mask(iv(1),iv(2),iv(3)) = max(mask(iv(1),iv(2),iv(3)),int(RTHL))

  enddo
88 continue

  !c Added by PB 2015.06.18 since a few peaks have RTHL < 1 which produced Nans otherwise
  if(naverage.eq.0) then
     sbar2 = 0.0
  else
     sbar2 = sbar2 / naverage
  endif
  
  !
  ! Added by M.A.A. 19/06/14
  ! --------------------------------------------------------

  
  deallocate(nshell,Sshell,SRshell,Fbar,rad)
  return

end subroutine get_homel

function hRinteg(u0,u12)
  USE intreal_types

  !c Integral of SPH kernel 4pi*W(x)dx*x(x^2-u12) from x0 to 2

  real(sp) hRinteg,u0,u12,x0,x02,x04,h1,h2,one7

  hRinteg=0.0
  x0=abs(u0)
  if(x0.ge.2.0) return
  hRinteg=1.4*u12-31.0/35.0
  if(x0.eq.0.0) return
  x02=x0*x0
  x04=x02*x02
  one7=1.0/7.0

  !c value at upper boundary x=2
  hRinteg=1.6*u12-6.4*one7
  if(x0.ge.1.0) then
     h1=x04*(2.0-x0*(2.4-x0*(1.0-x0*one7)))
     h2=x02*(4.0-x0*(4.0-x0*(1.5-x0*0.2)))
     hRinteg=hRinteg-h2*u12+h1
     return
  else
     hRinteg=hRinteg-0.2*(u12-one7)
     if(x0.gt.0.0) then
        h1=x04*(1.0-x02*(1.0-3.0*x0*one7))
        h2=x02*(2.0-x02*(1.5-0.6*x0))
        hRinteg=hRinteg-h2*u12+h1
     endif
  endif

  return
end function hRinteg


function hinteg(u0)
  USE intreal_types

  !c Integral of SPH kernel 4pi*W(x)dx*x from x0 to 2

  real(sp) hinteg,u0,x0,x02,x04,h1

  hinteg=0.0
  x0=abs(u0)
  if(x0.ge.2.0) return
  hinteg=1.4
  if(x0.eq.0.0) return
  x02=x0*x0
  x04=x02*x02

  !c value at upper boundary x=2
  hinteg=1.6
  if(x0.ge.1.0) then
     h1=x02*(4.0-x0*(4.0-x0*(1.5-x0*0.2)))
     hinteg=hinteg-h1
     return
  else
     hinteg=hinteg-0.2
     if(x0.gt.0.0) then
        h1=x02*(2.0-x02*(1.5-0.6*x0))
        hinteg=hinteg-h1
     endif
  endif

  return
end function hinteg

subroutine get_ijk(n1xn2,n1,j,j1,j2,j3)
  USE intreal_types
  !c Corrected version stm 3 Jun 1993
  jj=j-1
  !c this used to be just j
  j3=1+jj/n1xn2
  je=jj-n1xn2*(j3-1)
  j2=1+je/n1
  j1=je-n1*(j2-1)+1
  !c this used to not have the +1 (since it was j not j-1)
  return
end subroutine get_ijk

FUNCTION RAN1(IDUM)
  USE intreal_types
  DIMENSION R(97)
  PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
  PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
  PARAMETER (M3=243000,IA3=4561,IC3=51349)
  DATA IFF /0/
  save
  IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
     IFF=1
     IX1=MOD(IC1-IDUM,M1)
     IX1=MOD(IA1*IX1+IC1,M1)
     IX2=MOD(IX1,M2)
     IX1=MOD(IA1*IX1+IC1,M1)
     IX3=MOD(IX1,M3)
     DO J=1,97
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IA2*IX2+IC2,M2)
        R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
     ENDDO
     IDUM=1
  ENDIF
  IX1=MOD(IA1*IX1+IC1,M1)
  IX2=MOD(IA2*IX2+IC2,M2)
  IX3=MOD(IA3*IX3+IC3,M3)
  J=1+(97*IX3)/M3
  IF(J.GT.97.OR.J.LT.1) THEN
     write(*,*) 'ran1 problem '
     STOP
  ENDIF
  RAN1=R(J)
  R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
  RETURN
END FUNCTION RAN1

subroutine xran3_init(idum)
  USE intreal_types
  use cran3
  parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
  mj=mseed-iabs(idum)
  mj=mod(mj,mbig)
  ma(55)=mj
  mk=1
  do i=1,54
     ii=mod(21*i,55)
     ma(ii)=mk
     mk=mj-mk
     if(mk.lt.mz)mk=mk+mbig
     mj=ma(ii)
  enddo
  do k=1,4
     do i=1,55
        ma(i)=ma(i)-ma(1+mod(i+30,55))
        if(ma(i).lt.mz)ma(i)=ma(i)+mbig
     enddo
  enddo
  inext=0
  inextp=31
  if(idum.lt.0) idum=1
  return
end subroutine xran3_init

function xran3(idum)
  USE intreal_types
  use cran3
  !c     implicit real*4(m)
  !c     parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=2.5e-7)
  parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
  inext=inext+1
  if(inext.eq.56)inext=1
  inextp=inextp+1
  if(inextp.eq.56)inextp=1
  mj=ma(inext)-ma(inextp)
  if(mj.lt.mz)mj=mj+mbig
  ma(inext)=mj
  xran3=mj*fac
  return
end function xran3

subroutine xran3b_init(idum)
  USE intreal_types
  use cran3
  !c     implicit real*4(m)
  !c     parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=2.5e-7)
  parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
  mj=mseed-iabs(idum)
  mj=mod(mj,mbig)
  ma(55)=mj
  mk=1
  do i=1,54
     ii=mod(21*i,55)
     ma(ii)=mk
     mk=mj-mk
     if(mk.lt.mz)mk=mk+mbig
     mj=ma(ii)
  enddo
  do k=1,4
     do i=1,55
        ma(i)=ma(i)-ma(1+mod(i+30,55))
        if(ma(i).lt.mz)ma(i)=ma(i)+mbig
     enddo
  enddo
  inext=0
  inextp=31
  if(idum.lt.0) idum=1
  return
end subroutine xran3b_init

function xran3b(idum)
  USE intreal_types
  use cran3b
  !c     implicit real*4(m)
  !c     parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=2.5e-7)
  parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
  inext=inext+1
  if(inext.eq.56)inext=1
  inextp=inextp+1
  if(inextp.eq.56)inextp=1
  mj=ma(inext)-ma(inextp)
  if(mj.lt.mz)mj=mj+mbig
  ma(inext)=mj
  xran3b=mj*fac
  return
end function xran3b

!c** PRESS ETAL ROUTINE ******************
SUBROUTINE HUNT1(XX,N,X,JLO)
  USE intreal_types
  DIMENSION XX(N)
  LOGICAL ASCND
  save
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
END SUBROUTINE HUNT1

subroutine lagrint(xa,val,est,n,mdim,x)
  USE intreal_types
  dimension xa(n),est(mdim),val(n,mdim)
  !C Lagrange 4 pt interpolation or 2 point 

  if(n.eq.4) then
     a1=(x-xa(2))/(xa(1)-xa(2))*(x-xa(3))/(xa(1)-xa(3))*(x-xa(4))/(xa(1)-xa(4))
     a2=(x-xa(1))/(xa(2)-xa(1))*(x-xa(3))/(xa(2)-xa(3))*(x-xa(4))/(xa(2)-xa(4))
     a3=(x-xa(1))/(xa(3)-xa(1))*(x-xa(2))/(xa(3)-xa(2))*(x-xa(4))/(xa(3)-xa(4))
     a4=(x-xa(1))/(xa(4)-xa(1))*(x-xa(2))/(xa(4)-xa(2))*(x-xa(3))/(xa(4)-xa(3))
     do m=1,mdim
        est(m)=a1*val(1,m)+a2*val(2,m)+a3*val(3,m)+a4*val(4,m)
     enddo
  else 
     a1=(x-xa(2))/(xa(1)-xa(2))
     a2 = 1-a1
     do m=1,mdim
        est(m) =a1*val(1,m)+a2*val(2,m)
     enddo
  endif
  return
end subroutine lagrint

subroutine get_strain(l,m,n,n1,n2,n3,etax,etay,etaz,strain,order)

  implicit none

  integer i,j,k,n1,n2,n3,l,m,n,o,c,order
  real etax(n1,n2,n3),etay(n1,n2,n3),etaz(n1,n2,n3),strain(6)

  real e(-4:4,3,3),f(-4:4),etavec(3)

  real nth_order_derivative

  ! e(i,j,k) is jth component of ith neighbor in k direction


  do o=-4,4 ! o is for order

     i = l + o
     j = m + o
     k = n + o

     if(i>n1) i=i-n1
     if(j>n2) j=j-n2
     if(k>n3) k=k-n3
     if(i<1) i=i+n1
     if(j<1) j=j+n2
     if(k<1) k=k+n3

     ! direction 1 : increment index 1
     e(o,1,1) = etax(i,m,n) 
     e(o,2,1) = etay(i,m,n) 
     e(o,3,1) = etaz(i,m,n) 
     ! direction 2 : increment index 2
     e(o,1,2) = etax(l,j,n) 
     e(o,2,2) = etay(l,j,n) 
     e(o,3,2) = etaz(l,j,n) 
     ! direction 3 : increment index 3
     e(o,1,3) = etax(l,m,k) 
     e(o,2,3) = etay(l,m,k) 
     e(o,3,3) = etaz(l,m,k) 

  enddo

  f = e(:,1,1)
  strain(1) = nth_order_derivative(f,order)

  f = e(:,2,2)
  strain(2) = nth_order_derivative(f,order)

  f = e(:,3,3)
  strain(3) = nth_order_derivative(f,order)

  f = ( e(:,2,3) + e(:,3,2) ) /2.
  strain(4) = nth_order_derivative(f,order)

  f = ( e(:,1,3) + e(:,3,1) ) /2.
  strain(5) = nth_order_derivative(f,order)

  f = ( e(:,1,2) + e(:,2,1) ) /2.
  strain(6) = nth_order_derivative(f,order)

end subroutine get_strain

subroutine checkiv(iv,nn,i)

  integer iv(3),nn(3)

  do L=1,3
     if(iv(L)<1.or.iv(L)>nn(L)) then
        write(*,*) 'call = ',i
        write(*,*) 'iv out of bounds:',iv(L),L,nn(L)
        call mpi_finalize(ierr)
        stop
     endif
  enddo
  
end subroutine checkiv

subroutine checkm(m,n,i)

  if(m<1.or.m>n) then
     write(*,*) 'call = ',i
     write(*,*) 'out of bounds:',m,n
     call mpi_finalize(ierr)
     stop
  endif
  
end subroutine checkm

real function nth_order_derivative(f,order)

  implicit none

  integer order
  real f(-4:4),derivative

  if(order==8) then
     derivative = &
          1./280. * ( f(-4) - f(4) ) + &
          4./105. * (-f(-3) + f(3) ) + &
          1./5.   * ( f(-2) - f(2) ) + &
          4./5.   * (-f(-1) + f(1) ) 
  else
     derivative = &
          1./2. * (-f(-1) + f(1) ) 
  endif

  nth_order_derivative = derivative

  return 
  
end function nth_order_derivative

subroutine icloud(npart,nhunt)
  USE intreal_types
  use arrays  
  use mpivars

  integer, allocatable :: indx(:)

  ! THIS SUBROUTINE GENERATES Lattice POSITIONS 
  ! irs2 is the (r/alatt)^2 ixs,iys,izs are the positions in alatt units     

  n11=n1+1
  n21=n2+1
  n31=n3+1

  n1_21=n1/2+1
  n2_21=n2/2+1
  n3_21=n3/2+1

  irad2=nhunt**2 ! Here we set largest halo
  if(nhunt>nhalomax) then
    if(myid==0) write(*,*) 'nhunt > nhalomax, exiting',nhunt,nhalomax
    call mpi_finalize(ierr)
    stop
  endif
  m=0
  do  i = 1,n31
     izz=(i-n3_21)
     do  j = 1,n21
        iyy=(j-n2_21)
        do  k = 1,n11
           ixx=(k-n1_21)
           irr=ixx*ixx+iyy*iyy+izz*izz
           if(irr.le.irad2) then
              m = m + 1
           endif 
        enddo
     enddo
  enddo
  npart=m

  allocate(indx(npart+20))
  allocate(ixsvec(npart+20,3))
  allocate(irs2(npart+20))

  m=0
  do  i = 1,n31
     izz=(i-n3_21)
     do  j = 1,n21
        iyy=(j-n2_21)
        do  k = 1,n11
           ixx=(k-n1_21)
           irr=ixx*ixx+iyy*iyy+izz*izz
           if(irr.le.irad2) then
              m = m + 1
              irs2(m) = ixx*ixx+iyy*iyy+izz*izz
              indx(m)  = m
           endif 
        enddo
     enddo
  enddo

  call sort2(npart,irs2,indx)
  do i=1,npart
     irs2(indx(i))=i
  enddo

  mold=0
  do  i = 1,n31
    izz=(i-n3_21)
    do  j = 1,n21
       iyy=(j-n2_21)
       do  k = 1,n11
	  ixx=(k-n1_21)
           irr=ixx*ixx+iyy*iyy+izz*izz
           if(irr.le.irad2) then
              mold = mold + 1
              m = irs2(mold)

              ixsvec(m,1) = ixx
              ixsvec(m,2) = iyy
              ixsvec(m,3) = izz
           endif
        enddo
     enddo
  enddo

  do m=1,npart
    irs2(m) = sum(ixsvec(m,:)**2)
  enddo

  deallocate(indx)

  return
end subroutine icloud


subroutine atab4
  USE intreal_types
  use table
  parameter (two_thirds=2.0/3.0)
  pi=4.0*atan(1.0)

  !C  This subroutine provides the kernel akk and the derivative of the kernel,
  !C  dkk. The latter is divided by r. The grid density assignment kernel is akg.
  !C  The kernel akk is a generalization of schoenbergs m4 kernel with continuous 
  !C  derivatives up to the second. 

  ca=3./(2.*pi)
  cb=1./(4.*pi)

  do i=1,1100
     u=(i-1.)*0.01
     ru=sqrt(u)
     akk(i)=0.
     dkk(i)=0.
     if(ru.lt.1.)then
        akk(i)=ca*(two_thirds-u+0.5*u*ru)
        dkk(i)=ca*(-2.0+1.5*ru)
     else if(ru.ge.1.0.and.ru.le.2.0)then
        akk(i)=cb*(2.0-ru)**3
        dkk(i)=-3.0*cb*(2.0-ru)**2/ru
     end if
  enddo
  return
end subroutine atab4


subroutine get_evals(Nu,Qin,rdiag,iflag)
  USE intreal_types

  PARAMETER (nkrp=3)
  real(sp) Qin(nkrp,nkrp),rdiag(nkrp)
  real(dp) evec(nkrp,nkrp)
  real(dp) diag(nkrp),offdiag(nkrp)
  save nuphys
  data nuphys/3/

  iflag=0

  do ib=1,Nu
     do ia=1,Nu
        evec(ia,ib)=dble(Qin(ia,ib))
     enddo
  enddo

  call TRED2(evec,Nu,nuphys,diag,offdiag)
  call TQLInm(diag,offdiag,Nu,nuphys,evec,iflag)

  !C Kth COLUMN OF evec IS THE NORMALIZED EVECTOR OF EVALUE diag(K)
  !C   THAT IS   evec(ia,K),ia=1,3 
  !c IFLAG=0 normally, IFLAG=-1 for bad solution

  if (iflag.eq.0) then
     do ib=1,Nu
        rdiag(ib)=sngl(diag(ib))
        do ia=1,Nu
           Qin(ia,ib)=sngl(evec(ia,ib))
        enddo
     enddo
     call eval_sort(Nu,nuphys,rdiag,Qin)
  else
     !c TQLI failed
     do ib=1,Nu
        rdiag(ib)=-1
     enddo
  endif

  return
end subroutine get_evals


SUBROUTINE TRED2(A,N,NP,D,E)
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  REAL(DP) A(NP,NP),D(NP),E(NP)
  IF(N.GT.1)THEN
     DO I=N,2,-1  
        L=I-1
        H=0.0d0
        SCALE=0.0d0
        IF(L.GT.1)THEN
           DO K=1,L
              SCALE=SCALE+ABS(A(I,K))
           enddo
           IF(SCALE.EQ.0.)THEN
              E(I)=A(I,L)
           ELSE
              DO K=1,L
                 A(I,K)=A(I,K)/SCALE
                 H=H+A(I,K)**2
              enddo
              F=A(I,L)
              G=-SIGN(SQRT(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=0.0d0
              DO J=1,L
                 A(J,I)=A(I,J)/H
                 G=0.0d0
                 DO K=1,J
                    G=G+A(J,K)*A(I,K)
                 enddo
                 IF(L.GT.J)THEN
                    DO K=J+1,L
                       G=G+A(K,J)*A(I,K)
                    enddo
                 ENDIF
                 E(J)=G/H
                 F=F+E(J)*A(I,J)
              enddo

              HH=F/(H+H)
              DO J=1,L
                 F=A(I,J)
                 G=E(J)-HH*F
                 E(J)=G
                 DO K=1,J
                    A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
                 enddo
              enddo
           ENDIF
        ELSE
           E(I)=A(I,L)
        ENDIF
        D(I)=H
     enddo
  ENDIF
  D(1)=0.0d0
  E(1)=0.0d0
  DO I=1,N
     L=I-1
     IF(D(I).NE.0.)THEN
        DO J=1,L
           G=0.0d0
           DO K=1,L
              G=G+A(I,K)*A(K,J)
           enddo
           DO K=1,L
              A(K,J)=A(K,J)-G*A(K,I)
           enddo
        enddo
     ENDIF
     D(I)=A(I,I)
     A(I,I)=1.0d0
     IF(L.GE.1)THEN
        DO J=1,L
           A(I,J)=0.0d0
           A(J,I)=0.0d0
        enddo
     ENDIF
  enddo
  RETURN
END SUBROUTINE TRED2

SUBROUTINE TQLInm(D,E,N,NP,Z,iflag)
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  REAL(DP) D(NP),E(NP),Z(NP,NP)
  real(dp), parameter :: EPS=epsilon(x)
  iflag=0
  IF (N.GT.1) THEN
     DO I=2,N
        E(I-1)=E(I)
     enddo
     E(N)=0.0d0
     DO L=1,N
        ITER=0
1       DO M=L,N-1
           DD=ABS(D(M))+ABS(D(M+1))
           ! SEE http://www.nr.com/forum/showthread.php?p=5251
           ! IF (ABS(E(M))+DD.EQ.DD) GO TO 2 ! OLD WAY
           if(ABS(E(M))<=EPS*DD) GO TO 2     ! NEW WAY 
        enddo
        M=N
2       IF(M.NE.L)THEN
           !c            IF(ITER.EQ.30)PAUSE 'too many iterations'
           IF(ITER.EQ.60) THEN
              write(*,*) 'TQLI : too many iterations = ',ITER
              iflag=-1
              return
           ENDIF
           ITER=ITER+1
           G=(D(L+1)-D(L))/(2.0d0*E(L))
           R=SQRT(G**2+1.0d0)
           G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
           S=1.0d0
           C=1.0d0
           P=0.0d0
           DO I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(ABS(F).GE.ABS(G))THEN
                 C=G/F
                 R=SQRT(C**2+1.0d0)
                 E(I+1)=F*R
                 S=1.0d0/R
                 C=C*S
              ELSE
                 S=F/G
                 R=SQRT(S**2+1.0d0)
                 E(I+1)=G*R
                 C=1.0d0/R  
                 S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.0d0*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
              DO K=1,N
                 F=Z(K,I+1)
                 Z(K,I+1)=S*Z(K,I)+C*F
                 Z(K,I)=C*Z(K,I)-S*F
              enddo
           enddo
           D(L)=D(L)-P
           E(L)=G
           E(M)=0.0d0
           GO TO 1
        ENDIF
     enddo
  ENDIF
  RETURN
END SUBROUTINE TQLInm

subroutine eval_sort(Nu,nuphys,diag,evec)
  USE intreal_types
  !C SORTS THE EVALUES SO BIGGEST IN 3, SMALLEST IN 1
  !C SORTS THE EVECTORS AS WELL 
  parameter (n=3)
  real(sp) diag(nuphys),evec(nuphys,nuphys)
  real(sp) wksp(n)
  integer(i4b) iwksp(n)

  call rindexx(Nu,diag,iwksp)

  !C THE NEXT STEP CHANGES TO DECREASING ORDER FROM INCREASING ORDER
  !c      Nu_2=Nu/2
  !c      Nu1=Nu+1
  !c      do 10 j=1,Nu_2
  !c         idum=iwksp(j)
  !c         iwksp(j)=iwksp(Nu1-j)
  !c         iwksp(Nu1-j)=idum
  !c 10   continue

  do j=1,Nu
     wksp(j)=diag(j)
  enddo
  do j=1,Nu
     diag(j)=wksp(iwksp(j))
  enddo
  do j=1,Nu
     wksp(j)=evec(1,j)
  enddo
  do j=1,Nu
     evec(1,j)=wksp(iwksp(j))
  enddo
  do j=1,Nu
     wksp(j)=evec(2,j)
  enddo
  do j=1,Nu
     evec(2,j)=wksp(iwksp(j))
  enddo
  do j=1,Nu
     wksp(j)=evec(3,j)
  enddo
  do j=1,Nu
     evec(3,j)=wksp(iwksp(j))
  enddo
  return
end subroutine eval_sort


SUBROUTINE RINDEXX(N,ARRIN,INDX)
  USE intreal_types
  REAL(SP) ARRIN(*),Q
  INTEGER(I4B) INDX(*)
  DO J=1,N
     INDX(J)=J
  enddo
  L=N/2+1
  IR=N
10 CONTINUE
  IF(L.GT.1)THEN
     L=L-1
     INDXT=INDX(L)
     Q=ARRIN(INDXT)
  ELSE
     INDXT=INDX(IR)
     Q=ARRIN(INDXT)
     INDX(IR)=INDX(1)
     IR=IR-1
     IF(IR.EQ.1)THEN
        INDX(1)=INDXT
        RETURN
     ENDIF
  ENDIF
  I=L
  J=L+L
20 IF(J.LE.IR)THEN
     IF(J.LT.IR)THEN
        IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
     ENDIF
     IF(Q.LT.ARRIN(INDX(J)))THEN
        INDX(I)=INDX(J)
        I=J
        J=J+J
     ELSE
        J=IR+1
     ENDIF
     GO TO 20
  ENDIF
  INDX(I)=INDXT
  GO TO 10
END SUBROUTINE RINDEXX


SUBROUTINE sort2(n,arr,brr)
      INTEGER n,M,NSTACK
      INTEGER(kind=4) arr(n),brr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,b,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
          do 11 i=j-1,l,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
11        continue
          i=l-1
2         arr(i+1)=a
          brr(i+1)=b
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        temp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
          temp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
          temp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
          temp=brr(l)
          brr(l)=brr(l+1)
          brr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
        b=brr(l+1)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        temp=brr(i)
        brr(i)=brr(j)
        brr(j)=temp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        brr(l+1)=brr(j)
        brr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort2'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
END SUBROUTINE SORT2

SUBROUTINE ISORT4(N,RA,RB,RC,RD,WKSP,IWKSP)
  USE intreal_types
  INTEGER RA(N),RB(N),RC(N),RD(N),WKSP(N),IWKSP(N)
  CALL INDEXX(N,RA,IWKSP)
  DO J=1,N
     WKSP(J)=RA(J)
  enddo
  DO J=1,N
     RA(J)=WKSP(IWKSP(J))
  enddo
  DO J=1,N
     WKSP(J)=RB(J)
  enddo
  DO J=1,N
     RB(J)=WKSP(IWKSP(J))
  enddo
  DO J=1,N
     WKSP(J)=RC(J)
  enddo
  DO J=1,N
     RC(J)=WKSP(IWKSP(J))
  enddo
  DO J=1,N
     WKSP(J)=RD(J)
  enddo
  DO J=1,N
     RD(J)=WKSP(IWKSP(J))
  enddo
  RETURN
END SUBROUTINE ISORT4

SUBROUTINE INDEXX(N,ARRIN,INDX)
  USE intreal_types
  INTEGER ARRIN(N),INDX(N),Q

  DO J=1,N
     INDX(J)=J
  enddo
  L=N/2+1
  IR=N
10 CONTINUE
  IF(L.GT.1)THEN
     L=L-1
     INDXT=INDX(L)
     Q=ARRIN(INDXT)
  ELSE
     INDXT=INDX(IR)
     Q=ARRIN(INDXT)
     INDX(IR)=INDX(1)
     IR=IR-1
     IF(IR.EQ.1)THEN
        INDX(1)=INDXT
        RETURN
     ENDIF
  ENDIF
  I=L
  J=L+L
20 IF(J.LE.IR)THEN
     IF(J.LT.IR)THEN
        IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
     ENDIF
     IF(Q.LT.ARRIN(INDX(J)))THEN
        INDX(I)=INDX(J)
        I=J
        J=J+J
     ELSE
        J=IR+1
     ENDIF
     GO TO 20
  ENDIF
  INDX(I)=INDXT
  GO TO 10
END SUBROUTINE INDEXX
