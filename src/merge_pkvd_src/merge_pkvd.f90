program merge

  use mpivars

  implicit none 

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,ntasks,ierr)

  call  merge_pkvd

  call mpi_barrier(mpi_comm_world,ierr)

  call mpi_finalize(ierr)

end program merge

!C need to modify for the no-eovlution case, ie fixed redshift in the box?? 

subroutine merge_pkvd
  USE intreal_types
  use input_parameters
  use textlib
 !c modules needed for for Tree routines
  USE mparam
  use rpart
  USE ipart
  USE mnodes
  USE mbound
  USE mlists
  use mpivars
  !c Program to perform final box-seam trimming and application of
  !c Zeldovich dynamics to output catalog from SZ_PK.

  !c cita  stm     7 Jan 1991  - version new cluster format
  !c       stm     9 Jan 1991  - add cluster id
  !c       stm    11 Jan 1991  - version SZ_PK_MRG for merging
  !c       stm    14 Jan 1991  - box edge trim
  !c       stm    16 Feb 1991  - trim by RTH only option
  !c       stm    20 Feb 1991  - revised
  !c       stm    12 Mar 1991  - RTH trim one radius option
  !c       stm    13 Mar 1991  - revise RTH adjust
  !c       stm    16 Apr 1991  - bug fix, no boxfile option
  !c       stm    18 Apr 1991  - bug fix, adjust vTHvir with RTH change
  !c       stm    19 Jun 1991  - add negv, no outfile and .001 options
  !c ***   stm    25 Jun 1991  - version to merge in final config.
  !c       stm     7 Aug 1991  - correct velocity after zeldovich
  !c       stm     8 Aug 1991  - output in index order
  !c       stm    13 Aug 1991  - link-merging option, replaces RTH adjust
  !c       stm    21 Aug 1991  - independent index, trim, merge
  !c       stm     8 Nov 1991  - correct apportion of mass in merge
  !c       stm     9 Nov 1991  - trim,merge,link sequence revised/corrected
  !c       stm    17 Nov 1991  - new formats
  !c       stm     8 Jan 1992  - bug fixes for trim/merge (imrgtr=0)
  !c       stm    17 Jan 1992  - use IMRG vector to track links
  !c       stm    22 Jan 1992  - more options on reduction
  !c       stm     5 May 1992  - add Tree utilities
  !c       stm    12 Jun 1992  - heavily revised logic and options
  !c ***   stm    15 Jun 1992  - NEW WORKING VERSION, ran7
  !c       stm    16 Jun 1992  - new sorting
  !c       stm    15 Feb 1993  - version DEC-Alpha compatible
  !c                             streamline input, defaults
  !c       stm    17 Feb 1993  - min,max RTHL selection

  !C PROGRAM TRIMS CLUSTERS THAT OVERLAP
 
  !c declarations for merge
  parameter (nclmax=100)
  !c      parameter (nmax=200000)
  parameter (nmax=6000000)

  ! MEMORY ALLOCATIONS ESTIMATE
  ! num_arrays(nmax) = (9 + 6 + 6 + 2) + (9 + 4 + 3 + 6 + 6 + 2) + (7)
  !                  = 60
  ! So if ~1.5Gb per processor then nmax < 6.7e6 should work 
  real(sp) xon(nmax),yon(nmax),zon(nmax)
  real(sp) vxon(nmax),vyon(nmax),vzon(nmax)
  real(sp) vx2on(nmax),vy2on(nmax),vz2on(nmax)
  real(sp) strain_bar(6,nmax)
  real(sp) Fcollv(nmax),RTHLv(nmax),vTHvir(nmax)
  real(sp) hp(nmax),F_ev(nmax),F_pv(nmax)
  integer(i4b) iwas(nmax),lagrange(nmax)

  real(sp) xon2(nmax),yon2(nmax),zon2(nmax)
  real(sp) vxon2(nmax),vyon2(nmax),vzon2(nmax)
  real(sp) xout(nmax),yout(nmax),zout(nmax)
  real(sp) vxout(nmax),vyout(nmax),vzout(nmax),rout(nmax)
  real(sp) vx2on2(nmax),vy2on2(nmax),vz2on2(nmax)
  real(sp) strain_bar2(6,nmax)
  real(sp) Fcollv2(nmax),RTHLv2(nmax),vTHvir2(nmax)
  real(sp) hp2(nmax),F_ev2(nmax),F_pv2(nmax)
  integer(i4b) iwas2(nmax),lagrange2(nmax)

  integer(i4b) ibad(nmax),indxr(nmax)
  integer(i4b) imrg(nmax),lnk(nmax)
  real(sp) RTHfv(nmax),RTHL3(nmax)
  real(sp) aiwas(nmax)
  real(sp) s1LTHv(3),s2LTHv(3)

  integer(i4b) icindx(nclmax),nmrgv(nclmax),nbadv(nclmax)
  integer(i4b) idclus(nclmax),ic_cutv(nclmax),nicv(nclmax)
  real(sp) fcritv(nclmax),Rfclv(nclmax)
  real(sp) dLmaxv(nclmax),dFmaxv(nclmax),vtemp(nclmax),szeld(nclmax)
  real debugmeanx,debugmeany,debugmeanz,debugmeanvx,debugmeanvy,debugmeanvz,debugmeanRTH
  logical selectRTH
  character*128 outcode

  ! DOMAIN DECOMPOSITION VARIABLES
  integer, allocatable :: tilelist(:)

  real boxsize,Rbuff  
  integer itile,jtile,ktile,ntile,numtiles,tile,neach,tile_index
  real fxpk, fypk, fzpk
  real xpk, ypk, zpk, vxpk, vypk, vzpk, rpk, sbar1,sbar2,sbar3 
  integer Npk, Non2
  integer, allocatable :: Npk_eachtask(:), Npk_eachtaskl(:), Npk_begin(:)
  real    RTHLmin,  RTHLmax, RTHLmax_in
  real    RTHLminl, RTHLmaxl
  real, parameter :: epsilon=1e-6

  integer :: npk_tot
  integer :: outnum

  ! DATE AND TIME
  character(8)          :: dnt_date
  character(10)         :: dnt_time
  character(5)          :: dnt_zone
  integer, dimension(8) :: dnt_values
  integer initial_time, final_time, count_rate
  integer elapsed_hours, elapsed_minutes, elapsed_seconds
  double precision elapsed_time
  integer seedin
  integer(kind=8) pos_offset
  character*128 seedinstr


  ! Usage check
!  if(iargc()/=1) then
!     write(*,11)
!11   format('Usage: merge_pkvd <seedin>')
!     stop
!  endif

  ! REPORT DATE AND TIME
  if(myid==0) then
     call date_and_time(dnt_date,dnt_time,dnt_zone,dnt_values)
     write(*,12) dnt_values(3),dnt_values(2),&
          dnt_values(1),dnt_values(5),dnt_values(6),&
          dnt_values(7),dnt_zone
     call system_clock(initial_time,count_rate)
  endif

  seedin = i4arg(1,13579)
  ! Read parameters
  open(unit=1,file="merge_params.txt")
  read(1,*) boxsize,ntile,Omt,Omvac,h,iZeld,iLindx,&
       iLexc,iLmrg,iFindx,iFexc,iFmrg,ioutshear
  read(1,'(a)') filein
  read(1,'(a)') fileout
  close(1)

  numtiles  = ntile**3
  dcore_box = boxsize / ntile
  
  write (seedinstr, *) seedin
  seedinstr = adjustl(seedinstr)

  filein = trim(filein)//'.'//trim(seedinstr)
  fileout = trim(fileout)//'.'//trim(seedinstr)

  ! Hardwires
  selectRTH   = .false. ! don't exlcude based on RTHL value
  iwrap_merge = 0       ! periodic merging off
  nic         = 0       ! No cluster file
  iLagrp      = 0       ! keep Lagrange
  iredvth     = 1       ! change vTHvir on merge
  ikeeplnk    = 1       ! put link parent ID in Lagrange
  ncut_ic     = 0       ! don't cut out clusters
  f_disp      = 1       ! zel'dovich displacement factor = 1
  idovoid     = 0       ! don't do voids 
  fexc        = 1       ! final state exclusion factor
  fmrg        = 1       ! final state merging factor

  ! Initialization
  amtorad=atan(1.0)/45.0/60.0
  onethd=1.0/3.0
  do ic=1,nclmax
     nicv(ic)=0
     nmrgv(ic)=0
     nbadv(ic)=0
     icindx(ic)=0
     idclus(ic)=0
  enddo

  ! Random tile list
  iseed=13580
  allocate(tilelist(numtiles))
  do i=1,numtiles
    tilelist(i)=i
  enddo
  do i=numtiles,2,-1
    x = ran(iseed)
    j = int(x * i) + 1
    if(j<1.or.j>i) then
      write(*,*) 'j out of bounds',x,j
      call mpi_finalize(ierr)
      stop
    endif
    m = tilelist(j)
    tilelist(j) = tilelist(i)
    tilelist(i) = m    
  enddo
  
  ! Merging flags
  !  i(L,F)indx | 0-None, 1-RTHL, 2-RTHf, 3-Rf, 4-RTHL*
  !  i(L,F)exc  | 0-None, 1-(D<f*(R1+R2)) 
  !  i(L,F)mrg  | 0=None, 1-Link, 2-AddMass+mrg 

  ! Zeldovich flags
  !  iZeld      ! 0-None, 1-1LPT, 2-2LPT, >3-SEEJRB


  !open(unit=1,file=parameterfile)
  !read(1,*) Omt,Omvac,h        ! cosmology
  !read(1,*) iZeld              ! Zeldovich flag
  !read(1,*) iLindx,iLexc,iLmrg ! merging flags (Lagrangian)
  !read(1,*) iFindx,iFexc,iFmrg ! merging flags (Eulerian)
  !read(1,*) ioutshear          ! are Sbar2s in raw file
  !close(1)

  ! Set cosmology
  call read_cosmology(   Omt,h,Omvac,iamcurved,dcurv)
  call Dlinear_cosmology(Omt,h,Omvac,iamcurved,dcurv)
  
  ! THIS IS TO MAKE SURE THERE IS NO EXISTING FILE WITH THAT NAME
  if(myid==0) open(1,file=fileout,status='new',access='stream')
  call mpi_barrier(mpi_comm_world,ierr)
  if(myid.ne.0) open(1,file=fileout,status='unknown',access='stream')

  ! -----------------------------------------------------------------------
  ! LOOP OVER TILES
  ! -----------------------------------------------------------------------
  neach  = int(numtiles/ntasks)+1
  Npk_totl = 0
  RTHLmaxl = -1e10
  RTHLminl = 1e10

  do tile_index=myid+1,ntasks*neach,ntasks
  
  if(tile_index>numtiles) goto 201
  tile = tilelist(tile_index)

  ktile = (tile - 1 )                     / ntile**2 + 1
  jtile = (tile - ntile**2*(ktile-1) - 1) / ntile    + 1
  itile =  tile - ntile**2*(ktile-1) - (jtile-1)*ntile     

  open(4,file=filein,status='old',access='stream')
  read(4) Npk,RTHLmax_in

  Rbuff  = 2 * RTHLmax_in
  dL_box = dcore_box + 2 * Rbuff
  xtile = (itile - 0.5) * dcore_box - boxsize/2 
  ytile = (jtile - 0.5) * dcore_box - boxsize/2 
  ztile = (ktile - 0.5) * dcore_box - boxsize/2 

  ! DOMAIN DECOMPOSITION FILTERING OF INPUT DATA
  Non2 = 0
  if(ioutshear.eq.0) then
     do i=1,Npk
        
        read(4) xpk,ypk,zpk,vxpk,vypk,vzpk,rpk
        fxpk = abs(2*(xpk-xtile)/dL_box)
        if(fxpk<=1) then
           
           fypk = abs(2*(ypk-ytile)/dL_box)
           if(fypk<=1) then
              
              fzpk = abs(2*(zpk-ztile)/dL_box)
              if(fzpk<=1) then
                 
                 Non2         = Non2 + 1
                 xon2(Non2)   = xpk
                 yon2(Non2)   = ypk
                 zon2(Non2)   = zpk
                 vxon2(Non2)  = vxpk
                 vyon2(Non2)  = vypk
                 vzon2(Non2)  = vzpk
                 RTHLv2(Non2) = rpk
                 
              endif
           endif
        endif
     
     enddo
  elseif(ioutshear.eq.1) then
     do i=1,Npk
        
        read(4) xpk,ypk,zpk,vxpk,vypk,vzpk,rpk,sbar1,sbar2,sbar3
        fxpk = abs(2*(xpk-xtile)/dL_box)
        if(fxpk<=1) then
           
           fypk = abs(2*(ypk-ytile)/dL_box)
           if(fypk<=1) then
              
              fzpk = abs(2*(zpk-ztile)/dL_box)
              if(fzpk<=1) then
                 
                 Non2         = Non2 + 1
                 xon2(Non2)   = xpk
                 yon2(Non2)   = ypk
                 zon2(Non2)   = zpk
                 vxon2(Non2)  = vxpk
                 vyon2(Non2)  = vypk
                 vzon2(Non2)  = vzpk
                 RTHLv2(Non2) = rpk
                 strain_bar2(1,Non2) = sbar1
                 strain_bar2(2,Non2) = sbar2
                 strain_bar2(3,Non2) = sbar3
                 
              endif
           endif
        endif
     
     enddo
  endif
  close(4)

  if(Non2 > nmax) then
     write(*,*) 'Tiles too large, exiting... Non2, nmax, myid = ',Non2,nmax,myid
     call mpi_finalize(err)
     stop
  endif
     
  ! Read raw catalog
  !if(ioutshear.eq.0) then
  !   write(*,*) 'reading raw catalog, N = ',Non2,' RTHLmax = ',RTHLmax
  !   read(4) ( xon2(i), yon2(i), zon2(i),&
  !        vxon2(i),vyon2(i),vzon2(i),&
  !        RTHLv2(i),i=1,Non2) 
  !   close(4)
  !endif
  !if(ioutshear.eq.1) then
  !   write(*,*) 'reading raw catalog, N = ',Non2,' RTHLmax = ',RTHLmax
  !   read(4) ( xon2(i), yon2(i), zon2(i),&
  !        vxon2(i),vyon2(i),vzon2(i),&
  !        RTHLv2(i), (strain_bar2(j,i),j=1,3),i=1,Non2) 
  !   close(4)
  !endif


  ! Cull raw catalog
  Non=0
  jp=0
  nrejc=0
  nmrg=0
  do jp2=1,Non2
     if(iwas2(jp2).le.-10000.or.RTHLv2(jp2)<=0) then
        nrejc=nrejc+1
        goto 100
     endif
     if(iwas2(jp2).lt.0) nmrg=nmrg+1
     jp=jp+1
     xon(jp)=xon2(jp2)
     yon(jp)=yon2(jp2)
     zon(jp)=zon2(jp2)
     vxon(jp)=vxon2(jp2)
     vyon(jp)=vyon2(jp2)
     vzon(jp)=vzon2(jp2)
     vx2on(jp)=vx2on2(jp2)
     vy2on(jp)=vy2on2(jp2)
     vz2on(jp)=vz2on2(jp2)
     if(ioutshear.eq.1) then
        do kk=1,6
           strain_bar(kk,jp)=strain_bar2(kk,jp2)
        enddo
     endif
     Fcollv(jp)=Fcollv2(jp2)
     RTHLv(jp)=RTHLv2(jp2)
     vTHvir(jp)=vTHvir2(jp2)
     iwas(jp)=iwas2(jp2)
     lagrange(jp)=lagrange2(jp2)
     if(Fcollv(jp).gt.0.0001) then
        RTHfv(jp)=0.3*RTHLv(jp)/Fcollv(jp)
     else
        RTHfv(jp)=0.3*RTHLv(jp)/1.606
     endif
     if(iLmrg.eq.1.or.iFmrg.eq.2) RTHL3(jp)=RTHLv(jp)**3
     hp(jp)=hp2(jp2)
     F_ev(jp)=F_ev2(jp2)
     F_pv(jp)=F_pv2(jp2)
     id=abs(iwas(jp))
     if(nic.gt.0) then
        aiwas(jp)=Rfclv(id)
     else
        aiwas(jp)=-abs(iwas(jp))
     endif
     if(nic.gt.0.and.id.gt.0.and.id.le.nclmax) then
        ic=icindx(id)
        nicv(ic)=nicv(ic)+1
        if(iwas(jp).lt.0) nmrgv(ic)=nmrgv(ic)+1
     endif
100  continue
  enddo
  Non=jp

  npx=0
  npy=0
  npz=0
  do i=1,Non
     if(xon(i)>0) npx=npx+1
     if(yon(i)>0) npy=npy+1
     if(zon(i)>0) npz=npz+1
  end do

  ! MERGE
  ! --------------------------------------------------
  ! Here we are reversing the sign of strain_bar2
  ! so strain_bar_ij = -von_i,j
  !
  ! Change by M.A.A. 18/06/14
  ! --------------------------------------------------
  strain_bar = - strain_bar

  !c  Files read into merged arrays

  if(nic.gt.0) then
     do ic=1,nic
        nicv(ic)=0
        nmrgv(ic)=0
        nbadv(ic)=0
        dLmaxv(ic)=0.0
        dFmaxv(ic)=0.0
        vtemp(ic)=0.0
        szeld(ic)=0.0
     enddo
  endif

  !c  Get extrema of some useful quantities
  iquit=0
  radFmax=0.0
  radLmax=0.0
  vthrms=0.0
  szzrms=0.0
  nrms=0
  iwasbad=0
  irthbad=0
  do jp=1,Non
     if(iwas(jp).le.-10000.or.RTHLv(jp).le.0.0) then
        if(iwas(jp).le.-10000) iwasbad=iwasbad+1
        if(RTHLv(jp).le.0.0) irthbad=irthbad+1
        iquit=iquit+1
        ibad(jp)=1
     else
        ibad(jp)=0
        if(iwas(jp).lt.0) ibad(jp)=-1
        if(radFmax.lt.RTHfv(jp)) radFmax=RTHfv(jp)
        if(radLmax.lt.RTHLv(jp)) radLmax=RTHLv(jp)

        nrms=nrms+1
        vthrms=vthrms+vTHvir(jp)*vTHvir(jp)
        szzrms=szzrms+vxon(jp)**2+vyon(jp)**2+vzon(jp)**2

        id=abs(iwas(jp))
        if(nic.gt.0.and.id.gt.0.and.id.le.nclmax) then
           ic=icindx(id)
           if(iwas(jp).gt.0) then
              nicv(ic)=nicv(ic)+1
           else
              nmrgv(ic)=nmrgv(ic)+1
           endif
           dFmaxv(ic)=dFmaxv(ic)+RTHfv(jp)
           dLmaxv(ic)=dLmaxv(ic)+RTHLv(jp)

           vtemp(ic)=vtemp(ic)+vTHvir(jp)
           szeld(ic)=szeld(ic)+vxon(jp)**2+vyon(jp)**2+vzon(jp)**2
        endif
     endif
  enddo
  vthrms=sqrt(vthrms/float(nrms))
  szzrms=sqrt(szzrms/float(nrms))

  dLmax=2.0*radLmax
  dFmax=2.0*radFmax

  if(nic.gt.0) then
     do ic=1,iabs(nic)
        id=idclus(ic)
        if(nicv(ic).gt.0.or.nmrgv(ic).gt.0) then
           dLmaxv(ic)=dLmaxv(ic)/(nicv(ic)+nmrgv(ic))
           dFmaxv(ic)=dFmaxv(ic)/(nicv(ic)+nmrgv(ic))
           vtemp(ic)=vtemp(ic)/(nicv(ic)+nmrgv(ic))
           szeld(ic)=sqrt(szeld(ic)/(nicv(ic)+nmrgv(ic)))
        endif
!        write(*,5001) id,nicv(ic),nmrgv(ic),dLmaxv(ic),vtemp(ic),szeld(ic)
5001    format(1x,i3,2x,i6,2x,i6,3(2x,f8.2))
     enddo
  endif

  !c Index peaks (Lagrangian)
!  if(iLindx.eq.3) write(*,*) 'Lagrangian indexing by Rf rank '     
!  if(iLindx.eq.2) write(*,*) 'Lagrangian indexing by RTHf rank '
!  if(iLindx.eq.1) write(*,*) 'Lagrangian indexing by RTHL rank '
!  if(iLindx>3.or.iLindx<1) write(*,*) 'Lagrangian indexing by particle number '

  !c  Fill array INDXR with labels of points to check
  nindx=0
  do jp=1,Non
     if(ibad(jp).eq.1) go to 110
     if(selectRTH.and.(RTHLv(jp).lt.Rmgmin.or.RTHLv(jp).gt.Rmgmax)) then
        ibad(jp)=1
        go to 110
     endif
     if(iwas(jp).lt.0) ibad(jp)=-1
     nindx=nindx+1
     indxr(nindx)=jp
110  continue
  enddo
  ngood=nindx

  if(iLindx.ne.0) then

     !c  Sort INDXR entries in decreasing order of ...
     if(iLindx.eq.3) then
        call indexd3(nindx,aiwas,RTHLv,Fcollv,indxr)
     elseif(iLindx.eq.2) then
        call indexd3(nindx,RTHfv,RTHLv,Fcollv,indxr)
     else
        call indexd3(nindx,RTHLv,Fcollv,aiwas,indxr)
     endif

  endif

  !c  Lagrangian exclusion
!  if(iLexc.eq.4) then
!     write(*,*) 'Lagrangian exclusion for D < R1+R2 '
!  elseif(iLexc.eq.3) then
!     write(*,*) 'Lagrangian exclusion for D < R1 '
!  elseif(iLexc.eq.2) then
!     write(*,*) 'Lagrangian exclusion for D < R1-R2 '
!  elseif(iLexc.eq.1) then
!     write(*,*) 'Lagrangian exclusion for D = 0 '
!  else
!     write(*,*) 'No Lagrangian exclusion '
!  endif

  if(iLexc.ne.0) then

     !c  Build tree for trimming
     do jndx=1,nindx
        jp=indxr(jndx)
        !c  Store positions in r-array for Tree
        r(1,jndx)=xon(jp)
        r(2,jndx)=yon(jp)
        r(3,jndx)=zon(jp)
        
     enddo
     nobj=nindx
     call Bounds(nobj,r)
     call XToIx
     call BuildTree
     ! write(*,*) 'Lagrangian Oct-tree completed' 

     ntrim=0
     do jndx2=1,nindx-1
        !c  get list of all neighbors within dmax  
        jp2=indxr(jndx2)
        if(ibad(jp2).eq.1) goto 217
        if(iLindx.eq.1) then
           rnbr2=RTHLv(jp2)*RTHLv(jp2)
        else
           rnbr2=radLmax*radLmax
        endif
        if(iLexc.eq.4) rnbr2=4.0*rnbr2
        !c  call TreeSearch
        call Neighbors(jndx2,rnbr2)
        do jlst=1,nlist
           !c  check those ranked below test particle
           jndx1=ilist(jlst)
           if(jndx1.le.jndx2) goto 216
           jp1=indxr(jndx1)
           if(ibad(jp1).eq.1) goto 216
           !c Test:
           !c               write(18,*) jndx2,jndx1,sqrt(rlist(jlst))
           if(iLexc.eq.4) then
              test=(RTHLv(jp1)+RTHLv(jp2))**2
              if(rlist(jlst).gt.test) go to 216
           elseif(iLexc.eq.3) then
              test=(max(RTHLv(jp1),RTHLv(jp2)))**2
              if(rlist(jlst).gt.test) go to 216
           elseif(iLexc.eq.2) then
              test=(RTHLv(jp1)-RTHLv(jp2))**2
              if(rlist(jlst).gt.test) go to 216
           elseif(iLexc.eq.1) then
              if(rlist(jlst).gt.0.0) goto 216
           endif
           !c  Flag it --
           ibad(jp1)=1
           ntrim=ntrim+1
216        continue
        enddo
217     continue
     enddo
     ngood=ngood-ntrim
     if(myid==0) then
        write(*,*) 'Lagrangian exclusion completed for box ',tile_index
        write(*,*) ' number   trimmed = ',ntrim
        write(*,*) ' number remaining = ',ngood
     endif
  endif

  !c  Lagrangian merging
!  if(iLmrg.eq.1) then
!     write(*,*) 'Binary Lagrangian merging with reducued mass.'
!  else
!     write(*,*) 'No Lagrangian merging '
!  endif

  if(iLmrg.gt.0) then

     !c  Build tree for trimming (if not already)
     if(iLexc.eq.0) then
        do jndx=1,nindx
           jp=indxr(jndx)
           !c  Store positions in r-array for Tree
           r(1,jndx)=xon(jp)
           r(2,jndx)=yon(jp)
           r(3,jndx)=zon(jp)
        enddo
        nobj=nindx

        call Bounds(nobj,r)
        call XToIx
        call BuildTree
!        write(*,*) 'Lagrangian Oct-tree completed' 
     endif

     nmerg=0
     do jndx2=1,nindx-1
        !c  get list of all neighbors within dmax  
        jp2=indxr(jndx2)
        if(ibad(jp2).eq.1) goto 219
        if(iLindx.eq.1) then
           rnbr2=RTHLv(jp2)*RTHLv(jp2)
        else
           rnbr2=radLmax*radLmax
        endif
        rnbr2=4.0*rnbr2
        !c  call TreeSearch
        call Neighbors(jndx2,rnbr2)
        do jlst=1,nlist
           !c  check those ranked below test particle
           jndx1=ilist(jlst)
           if(jndx1.le.jndx2) goto 218
           jp1=indxr(jndx1)
           if(ibad(jp1).eq.1) goto 218
           !c  reduce mass by overlaps                        
           test=(RTHLv(jp1)+RTHLv(jp2))**2
           rad2=(xon(jp2)-xon(jp1))**2+(yon(jp2)-yon(jp1))**2+(zon(jp2)-zon(jp1))**2
           if(rad2.gt.test) goto 218
           rad1=sqrt(rad2)
           if(RTHLv(jp1).ge.RTHLv(jp2)) then
              jpA=jp1
              jpB=jp2
           else
              jpA=jp2
              jpB=jp1
           endif
           rA=RTHLv(jpA)
           rB=RTHLv(jpB)
           rA2=rA*rA
           rB2=rB*rB
           if(rad2.ge.(rA2-rB2))then
              d_A=0.5*(rad2+rA2-rB2)/rad1
              d_B=rad1-d_A
              h_A=d_A/rA
              h_B=d_B/rB
              dmmA=0.5-0.75*h_A+0.25*h_A**3
              dmmB=0.5-0.75*h_B+0.25*h_B**3
              RTHL3(jpA)=RTHL3(jpA)-dmmA*rA**3
              RTHL3(jpB)=RTHL3(jpB)-dmmB*rB**3
              if(ibad(jpA).eq.0) nmerg=nmerg+1
              ibad(jpA)=-1
              if(ibad(jpB).eq.0) nmerg=nmerg+1
              ibad(jpB)=-1
           elseif(rad1.gt.(rA-rB)) then
              d_A=0.5*(rad2+rA2-rB2)/rad1
              d_B=d_A-rad1
              h_A=d_A/rA
              h_B=d_B/rB
              dmmA=0.5-0.75*h_A+0.25*h_A**3
              !c  Now dmmB=1-dmmB(above)
              dmmB=0.5+0.75*h_B-0.25*h_B**3
              RTHL3(jpA)=RTHL3(jpA)-dmmA*rA**3
              RTHL3(jpB)=RTHL3(jpB)-dmmB*rB**3
              if(ibad(jpA).eq.0) nmerg=nmerg+1
              ibad(jpA)=-1
              if(ibad(jpB).eq.0) nmerg=nmerg+1
              ibad(jpB)=-1
           else
              !c  Smaller wholly inside larger
              RTHL3(jpB)=0.0
              if(ibad(jpB).eq.0) nmerg=nmerg+1
              ibad(jpB)=-1
           endif
218        continue
        enddo
219     continue
     enddo

     !c  Clean out the destroyed clusters

     ndest=0
     do jndx=1,nindx
        jp=indxr(jndx)
        if(ibad(jp).ge.0) goto 230
        if(RTHL3(jp).le.0.0) then
           ibad(jp)=1
           ndest=ndest+1
        endif
230     continue
     enddo

     nunmg=ngood-nmerg
     ngood=ngood-ndest
     if(myid==0) then
        write(*,*) 'Lagrangian merging completed for box ', tile_index
        write(*,*) ' number    merged = ',nmerg
        write(*,*) ' number destroyed = ',ndest
        write(*,*) ' number remaining = ',ngood
        write(*,*) ' number  unmerged = ',nunmg
     endif
  endif

  !c Report Lagrangian results.

  if(iLexc.ne.0.or.iLmrg.ne.0) then

     ngoodt=0
     nbadt=0
     nchg=0
     if(nic.gt.0) then
        do ic=1,nic
           nicv(ic)=0
           nmrgv(ic)=0
           nbadv(ic)=0
        enddo
     endif

     do jp=1,Non
        if(ibad(jp).eq.0) then
           ngoodt=ngoodt+1
        else if(ibad(jp).eq.-1) then
           nmrgt=nmrgt+1
        else
           nbadt=nbadt+1
        endif

        id=abs(iwas(jp))
        if(nic.gt.0.and.id.gt.0.and.id.le.nclmax) then
           ic=icindx(id)
           if(ibad(jp).eq.0) then
              nicv(ic)=nicv(ic)+1
           else if(ibad(jp).eq.-1) then
              nmrgv(ic)=nmrgv(ic)+1
           else if(ibad(jp).eq.1) then
              nbadv(ic)=nbadv(ic)+1
           endif
        endif
     enddo
     if(myid==0) then
        write(*,*) ' there are ', nbadt, ' total trimmed peaks '
        write(*,*) ' there are ', nmrgt, ' total merged peaks '
        write(*,*) ' there are ', ngoodt, ' total original peaks '
     endif
     if(nic.gt.0) then
!        write (*,*) ' id, nicv, nmrgv, nbadv '
        do ic=1,iabs(nic)
           id=idclus(ic)
!           write(*,*) id,nicv(ic),nmrgv(ic),nbadv(ic)
           nicv(ic)=0
           nmrgv(ic)=0
           nbadv(ic)=0
        enddo
     endif

  endif

  if(iFexc.ne.0.or.iFmrg.ne.0) then
     !c  Final-state exclusion and merging
     !c  Index peaks (Final state)
!     if(iFindx.eq.4) then
!        write(*,*) 'Final state indexing by current RTHL rank '
!     elseif(iFindx.eq.3) then
!        write(*,*) 'Final state indexing by Rf rank '
!     elseif(iFindx.eq.2) then
!        write(*,*) 'Final state indexing by RTHf rank '
!     elseif(iFindx.eq.1) then
!        write(*,*) 'Final state indexing by original RTHL rank '
!     else
!        write(*,*) 'Final state indexing by particle number '
!     endif

     !c  Fill array INDXR with labels of points to check
     nindx=0
     do jp=1,Non
        if(ibad(jp).eq.1) go to 313
        if(selectRTH.and.(RTHLv(jp).lt.Rmgmin.or.RTHLv(jp).gt.Rmgmax)) then
           ibad(jp)=1
           go to 313
        endif
        nindx=nindx+1
        indxr(nindx)=jp
        if(iLagrp.eq.1) lagrange(jp)=nindx
313     continue
     enddo
!     write(*,*) 'Final-state clusters selected = ',nindx
     ngood=nindx

     if(iFindx.ne.0) then

        !c  Sort INDXR entries in decreasing order of ...
        if(iFindx.eq.4.and.iLmrg.eq.1) then
           call indexd3(nindx,RTHL3,RTHLv,Fcollv,indxr)
        elseif(iFindx.eq.3) then
           call indexd3(nindx,aiwas,RTHLv,Fcollv,indxr)
        elseif(iFindx.eq.2) then
           call indexd3(nindx,RTHfv,RTHLv,Fcollv,indxr)
        else
           call indexd3(nindx,RTHLv,Fcollv,aiwas,indxr)
        endif

!        write(*,*) 'Final-state Index-ordering completed.'
     endif
   endif

  !c  Final-state (Eulerian) exclusion
 ! write(*,*)
 ! if(iFexc.eq.1) then
 !    write(*,*) 'Final-state exclusion for D < fR*(R1+R2) '
 ! else
 !    write(*,*) 'No Final-state exclusion '
 ! endif

  if(iFexc.gt.0) then

     !c  Build tree for trimming
     do jndx=1,nindx
        jp=indxr(jndx)
        !c  Store positions in r-array for Tree
        r(1,jndx)=xon(jp)-vxon(jp)
        r(2,jndx)=yon(jp)-vyon(jp)
        r(3,jndx)=zon(jp)-vzon(jp)
     enddo
     nobj=nindx

     call Bounds(nobj,r)
     call XToIx
     call BuildTree
  !   write(*,*) 'Final-state Oct-tree completed' 

     ntrim=0
     do jndx2=1,nindx-1
        !c  get list of all neighbors within dmax  
        jp2=indxr(jndx2)
        if(ibad(jp2).eq.1) goto 317
        if(iFindx.eq.2) then
           rnbr2=RTHfv(jp2)*RTHfv(jp2)
        else
           rnbr2=radFmax*radFmax
        endif
        rnbr2=4.0*fexc*fexc*rnbr2
        !c  call TreeSearch
        call Neighbors(jndx2,rnbr2)
        do jlst=1,nlist
           !c  check those ranked below test particle
           jndx1=ilist(jlst)
           if(jndx1.le.jndx2) goto 316
           jp1=indxr(jndx1)

           test=(fexc*(RTHfv(jp1)+RTHfv(jp2)))**2
           if(rlist(jlst).gt.test) go to 316
           !c  Flag it --
           ibad(jp1)=1
           ntrim=ntrim+1
316        continue
        enddo
317     continue
     enddo
     ngood=ngood-ntrim
     !write(*,*) 'Final-state exclusion completed'
     !write(*,*) ' number   trimmed = ',ntrim
     !write(*,*) ' number remaining = ',ngood
  endif

  !c  Final-state (Eulerian) merging
!  write(*,*)
!  if(iFmrg.eq.2) then
!     write(*,*) 'Final-state full merging '
!  elseif(iFmrg.eq.1) then
!     write(*,*) 'Final-state linking '
!  else
!     write(*,*) 'No Final-state merging '
!  endif

  if(iFmrg.gt.0) then

     if(iFexc.eq.0) then
        !c  Build tree for trimming if none already
        do jndx=1,nindx
           jp=indxr(jndx)
           !c  Store positions in r-array for Tree
           r(1,jndx)=xon(jp)-vxon(jp)
           r(2,jndx)=yon(jp)-vyon(jp)
           r(3,jndx)=zon(jp)-vzon(jp)
        enddo
        nobj=nindx

        call Bounds(nobj,r)
        call XToIx
        call BuildTree
!        write(*,*) 'Final-state Oct-tree completed' 
     endif

     !c  Make lists for merge
     do jndx=1,nindx
        jp=indxr(jndx)
        imrg(jp)=0
        lnk(jndx)=0
     enddo
     nobj=nindx

     nlnk=0
     nbadmg=0
     nthrow=0
     nlinks=0
     do jndx=1,nindx
        jp=indxr(jndx)
        !c  Link new cluster
        if(imrg(jp).ne.0) goto 322
        !c  Pop first nolnk to lnk. 
        lnk(1)=jndx
        nlnk=1
        nlnkp=0
        npop=1
        !c  npop counts number popped this round
        do while (npop.gt.0) 
           !c  if something popped, check again
           !c  check nolink list 
           npop=0
           !c  move through particles linked previous round
           ilnlow=nlnkp+1
           ilnupp=nlnk
           do ilnk=ilnlow,ilnupp
              jndx2=lnk(ilnk)
              jp2=indxr(jndx2)
              !c  get list of all neighbors within distance
              if(ibad(jp2).eq.1) goto 320
              if(iFindx.eq.2) then
                 rnbr2=RTHfv(jp2)*RTHfv(jp2)
              else
                 rnbr2=radFmax*radFmax
              endif
              rnbr2=4.0*fmrg*fmrg*rnbr2
              !c  call TreeSearch
              call Neighbors(jndx2,rnbr2)
              !c  add those not already linked
              ipop=0
              do ilst=1,nlist
                 jndx1=ilist(ilst)
                 jp1=indxr(jndx1)
                 if(imrg(jp1).ne.0) goto 319
                 rad2=rlist(ilst)
                 test=(fmrg*(RTHfv(jp1)+RTHfv(jp2)))**2
                 if(rad2.gt.test) goto 319
                 !c  found a merge
                 ipop=1
                 npop=npop+1
                 lnk(nlnk+npop)=jndx1
                 imrg(jp1)=1
319              continue
              enddo
320           continue
           enddo

           nlnkp=nlnk
           nlnk=nlnk+npop
        end do
        !c  no more linked
        !c  collapse lnk to lagrange designation of lnk(1)
        !c  sum linked masses if iFmrg=2
        if(nlnk.gt.1) then
           jndx1=lnk(1)
           jp1=indxr(jndx1)
           lagrp1=lagrange(jp1)
           do ilnk=2,nlnk
              jndx2=lnk(ilnk)
              jp2=indxr(jndx2)
              if(ikeeplnk.eq.1) lagrange(jp2)=lagrp1
              if(iFmrg.eq.2) then
                 !c  sum masses and delete mergers
                 RTHL3(jp1)=RTHL3(jp1)+RTHL3(jp2)
                 ibad(jp2)=1
              else
                 !c  just flag as altered
                 ibad(jp2)=-1
              endif
              nbadmg=nbadmg+1
           enddo
        endif
        nthrow=nthrow+nlnk-1
        nlinks=nlinks+1
322     continue
     enddo
     if(myid==0) then
        write(*,*) 'Final-state merging completed'
        write(*,*) ' nthrow,nlinks = ',nthrow,nlinks
        write(*,*) ' number clusters merged = ',nbadmg
        if(iFmrg.eq.2) write(*,*) ' Summed cluster masses on link'
     endif
  endif

  if(iFexc.ne.0.or.iFmrg.ne.0) then

     ngoodt=0
     nbadt=0
     nmrg=0
     if(nic.gt.0) then
        do ic=1,nic
           nicv(ic)=0
           nmrgv(ic)=0
           nbadv(ic)=0
        enddo
     endif

     do jp=1,Non
        if(ibad(jp).eq.-1) then
           nmrg=nmrg+1
        elseif(ibad(jp).eq.0) then
           ngoodt=ngoodt+1
        else
           nbadt=nbadt+1
        endif

        id=abs(iwas(jp))
        if(nic.gt.0.and.id.gt.0.and.id.le.nclmax) then
           ic=icindx(id)
           if(ibad(jp).eq.-1) then
              nmrgv(ic)=nmrgv(ic)+1
           else if(ibad(jp).eq.0) then
              nicv(ic)=nicv(ic)+1
           else if(ibad(jp).eq.1) then
              nbadv(ic)=nbadv(ic)+1
           endif
        endif
     enddo
     if(myid==0) then
        write(*,*) ' there are ', nbadt, ' trimmed peaks '
        write(*,*) ' there are ', nmrg, ' merged peaks '
        write(*,*) ' there are ', ngoodt, ' unmerged peaks '
     endif
     if(nic.gt.0) then
!        write (*,*) ' id, nomrgv, nmrgv, nbadv '
        do ic=1,nic
           id=idclus(ic)
!           write(*,*) id,nicv(ic),nmrgv(ic),nbadv(ic)
           nicv(ic)=0
           nmrgv(ic)=0
           nbadv(ic)=0
        enddo
     endif

  endif

  !c Write out in indexed order
  chin=6000.0/h
  if(iZeld.ne.0) then
!     write(*,*) 
!     write(*,*) 'Applying Zeldovich displacements, iZeld 1 zeld, 2 2nd order Lagrangian, 3,4,5 variants '
!     write(*,*) 'Factor = ',f_disp
  
     vrms=0.0
     jp2=0
     do jndx=1,nindx
        jp=indxr(jndx)
        if(ibad(jp).eq.1) goto 120
        if(iFmrg.eq.2.and.ibad(jp).eq.-1) goto 120
        !G.S. 13/11/2015
        !Added to toss peaks before move, else have a chance of peaks moving
        !further than a buffersize and being lost
        fxpk = abs(2*(xon(jp)-xtile)/dcore_box)
        fypk = abs(2*(yon(jp)-ytile)/dcore_box)
        fzpk = abs(2*(zon(jp)-ztile)/dcore_box)
        if(fxpk>=1) goto 120
        if(fypk>=1) goto 120
        if(fzpk>=1) goto 120
        xp2=(xon(jp)-vxon(jp))**2
        yp2=(yon(jp)-vyon(jp))**2
        zp2=(zon(jp)-vzon(jp))**2
        rr=xp2+yp2+zp2
        jp2=jp2+1
        !C THIS SHOULD ONLY BE DONE WITH  EVOLUTION
        !C FOR NO EVOLUTION ONE z OF THE BOX WITH ONE HDa IS BETTER
        chib=chifn_Sigma(sqrt(rr),iamcurved,dcurv)
        ab=afn(chib)
        Zb=1.0/ab
        !c  (1.0-chib/chin)**2
        Db=Dlinear(ab,chib,HD_Ha,D_a)
        hubb=100.0*hub_Dlin(ab)
        !c            HD=100.0*h0*Zb*sqrt(Zb)
        HD=HD_Ha*hubb
        HDa=HD/Zb
        s1Lthv(1)=f_disp*vxon(jp)
        s1Lthv(2)=f_disp*vyon(jp)
        s1Lthv(3)=f_disp*vzon(jp)
        s2Lthv(:)=0.0
        if(iZeld.gt.1) then
           theta=2.0
! theta of bpkan1 is probably much smaller, simple model had 3R_pk^2/r_cv^2/(1+(1+n_eff)(3+2*n_eff)) relating to velocity coherence scale, gamma = 3+n_eff
           fcollD=-(strain_bar(1,jp)+strain_bar(2,jp)+strain_bar(3,jp))
           ! we check here to make sure the convention of strain_bar is 
           ! the deformation tensor of the displacement field, -von
           ! i.e. tr(strain_bar)<0 for an overdense region
           ! M.A.A. 18/06/14
           
           !c P.B. 2015.06.18 Switched to local shear measurement
           !c if(fcollD<0) then
           !c   write(*,*) 'ERROR, fcolld < 0, exiting...'
           !c   stop
           !c endif
           if(iZeld.eq.2) then        
              !c facnl=f_disp*theta/14.0*3.0/2.0
              !c s2Lthv(1)=facnl*fcollD*vxon(jp)+facnl*(strain_bar(1,jp)*vxon(jp)+strain_bar(5,jp)*vzon(jp)+strain_bar(6,jp)*vyon(jp))
              !c s2Lthv(2)=facnl*fcollD*vyon(jp)+facnl*(strain_bar(2,jp)*vyon(jp)+strain_bar(4,jp)*vzon(jp)+strain_bar(6,jp)*vxon(jp))
              !c s2Lthv(3)=facnl*fcollD*vzon(jp)+facnl*(strain_bar(3,jp)*vzon(jp)+strain_bar(4,jp)*vyon(jp)+strain_bar(5,jp)*vxon(jp))
              s2Lthv(1) = strain_bar(1, jp)
              s2Lthv(2) = strain_bar(2, jp)
              s2Lthv(3) = strain_bar(3, jp)
           elseif(iZeld.eq.6) then
              s2Lthv(1)=f_disp*vx2on(jp)
              s2Lthv(2)=f_disp*vy2on(jp)
              s2Lthv(3)=f_disp*vz2on(jp)
           else
              if(iZeld.eq.3) then
                 facnl=f_disp*theta/14.0*(1.0-3.0*F_ev(kon)**2-F_pv(kon)**2)
              elseif(iZeld.eq.4) then
                 facnl=f_disp*theta/14.0*(1.0-sqrt(3.0*F_ev(kon)**2+F_pv(kon)**2))
              elseif(iZeld.eq.5) then
                 facnl=f_disp*theta/14.0
              endif              
              s2Lthv(1)=-facnl*fcollD*vxon(jp)
              s2Lthv(2)=-facnl*fcollD*vyon(jp)
              s2Lthv(3)=-facnl*fcollD*vzon(jp)
           endif
        endif
        xon2(jp2)=xon(jp)-s1Lthv(1)-s2Lthv(1)           
        yon2(jp2)=yon(jp)-s1Lthv(2)-s2Lthv(2)           
        zon2(jp2)=zon(jp)-s1Lthv(3)-s2Lthv(3)           
        vxon2(jp2)=-HDa*s1Lthv(1)-2.0*HDa*s2Lthv(1)           
        vyon2(jp2)=-HDa*s1Lthv(2)-2.0*HDa*s2Lthv(2)           
        vzon2(jp2)=-HDa*s1Lthv(3)-2.0*HDa*s2Lthv(3)           
        vrms=vrms+((xon2(jp2)-xon(jp))**2+(yon2(jp2)-yon(jp))**2+(zon2(jp2)-zon(jp))**2)
        Fcollv2(jp2)=Fcollv(jp)
        RTHLv2(jp2)=RTHLv(jp)
        vTHvir2(jp2)=vTHvir(jp)
        
        if(iLmrg.eq.1.or.iFmrg.eq.2) then
           RTHLv2(jp2)=RTHL3(jp)**onethd
           if(iredvth.eq.1) vTHvir2(jp2)=vTHvir(jp)*RTHL3(jp)**onethd/RTHLv(jp)
        endif
        if(ibad(jp).eq.-1) then
           iwas2(jp2)=-iabs(iwas(jp))
        else
           iwas2(jp2)=iabs(iwas(jp))
        endif
        if(iLagrp.eq.1.and.iFmrg.eq.0) then
           lagrange2(jp2)=jndx
        else
           lagrange2(jp2)=lagrange(jp)
        endif
        if(ifmt.eq.2.or.ifmt.eq.3) hp2(jp2)=hp(jp)
        if(ifmt.eq.3) then
           F_ev2(jp2)=F_ev(jp)
           F_pv2(jp2)=F_pv(jp)
        endif
        id=abs(iwas(jp))
        if(nic.gt.0.and.id.gt.0.and.id.le.nclmax) then
           ic=icindx(id)
           if(ibad(jp).eq.0) then
              nicv(ic)=nicv(ic)+1
           else if(ibad(jp).eq.-1) then
              nmrgv(ic)=nmrgv(ic)+1
           endif
        endif
120     continue
     enddo
     Non2=jp2
     vrms=sqrt(vrms/float(Non2))
!     write(*,*) ' Zeldovich displacement rms = ',vrms
  else
     jp2=0
     do jndx=1,nindx
        jp=indxr(jndx)
        if(ibad(jp).eq.1) goto 123
        if(iFmrg.eq.2.and.ibad(jp).eq.-1) goto 123
        fxpk = abs(2*(xon(jp)-xtile)/dcore_box)
        fypk = abs(2*(yon(jp)-ytile)/dcore_box)
        fzpk = abs(2*(zon(jp)-ztile)/dcore_box)
        if(fxpk>=1) goto 123
        if(fypk>=1) goto 123
        if(fzpk>=1) goto 123
        jp2=jp2+1
        xon2(jp2)=xon(jp)
        yon2(jp2)=yon(jp)
        zon2(jp2)=zon(jp)
        vxon2(jp2)=vxon(jp)
        vyon2(jp2)=vyon(jp)
        vzon2(jp2)=vzon(jp)
        if(ioutshear.eq.1) then
           do kk=1,6
              strain_bar2(kk,jp2)=strain_bar(kk,jp)
           enddo
        endif
        Fcollv2(jp2)=Fcollv(jp)
        RTHLv2(jp2)=RTHLv(jp)
        vTHvir2(jp2)=vTHvir(jp)
        if(iLmrg.eq.1.or.iFmrg.eq.2) then
           RTHLv2(jp2)=RTHL3(jp)**onethd
           if(iredvth.eq.1) vTHvir2(jp2)=vTHvir(jp)*RTHL3(jp)**onethd/RTHLv(jp)
        endif
        if(ibad(jp).eq.-1) then
           iwas2(jp2)=-iabs(iwas(jp))
        else
           iwas2(jp2)=iabs(iwas(jp))
        endif
        if(iLagrp.eq.1.and.iFmrg.eq.0) then
           lagrange2(jp2)=jndx
        else
           lagrange2(jp2)=lagrange(jp)
        endif
        if(ifmt.eq.2.or.ifmt.eq.3.or.ifmt.eq.5) hp2(jp2)=hp(jp)
        if(ifmt.eq.3.or.ifmt.eq.5) then
           F_ev2(jp2)=F_ev(jp)
           F_pv2(jp2)=F_pv(jp)
        endif
        id=abs(iwas(jp))
        if(nic.gt.0.and.id.gt.0.and.id.le.nclmax) then
           ic=icindx(id)
           if(ibad(jp).eq.0) then
              nicv(ic)=nicv(ic)+1
           else if(ibad(jp).eq.-1) then
              nmrgv(ic)=nmrgv(ic)+1
           endif
        endif
123     continue
     enddo
     Non2=jp2
  endif

  ! DOMAIN DECOMPOSITION FILTERING OF OUTPUT DATA
  Npk = Non2
  m = Npk_totl
  do i=1,Npk

     m            = m + 1
     xout(m)   = xon2(i)
     yout(m)   = yon2(i)
     zout(m)   = zon2(i)
     vxout(m)  = vxon2(i)
     vyout(m)  = vyon2(i)
     vzout(m)  = vzon2(i)
     rout(m)   = RTHLv2(i)
              
  enddo
  Npk_totl = m

  enddo
  ! -----------------------------------------------------------------------
  ! END LOOP OVER TILES 
  ! -----------------------------------------------------------------------

  201 continue

  call mpi_barrier(mpi_comm_world,ierr)

  call mpi_allreduce(Npk_totl,Npk_tot,1,mpi_integer,mpi_sum,&
                      mpi_comm_world,ierr)

  if(myid==0) write(*,102) Npk_tot

  allocate(Npk_eachtaskl(0:ntasks-1))
  allocate(Npk_eachtask (0:ntasks-1))
  allocate(Npk_begin    (0:ntasks-1))

  Npk_begin           = 0
  Npk_eachtaskl       = 0
  Npk_eachtaskl(myid) = Npk_totl

  call mpi_allreduce(Npk_eachtaskl,Npk_eachtask,ntasks,mpi_integer,mpi_sum,&
                      mpi_comm_world,ierr)

  Npk_begin(0) = 0
  do i=1,ntasks-1
     Npk_begin(i)=Npk_eachtask(i-1)+Npk_begin(i-1)
  enddo

  RTHLminl = min(RTHLminl,minval(rout))
  RTHLmaxl = max(RTHLmaxl,maxval(rout))

  ! WRITE PEAKS 
  pos_offset = int(4,8)*(2+7*Npk_begin(myid))+1
  write(1,pos=pos_offset) (  xout(i), yout(i), zout(i),&
            vxout(i),vyout(i),vzout(i),&
            rout(i),i=1,Npk_totl)

  deallocate(Npk_eachtaskl)
  deallocate(Npk_eachtask )
  deallocate(Npk_begin    )


  call mpi_reduce(RTHLminl,RTHLmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(RTHLmaxl,RTHLmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)

  ! WRITE HEADER IF myid == 0       
  if(myid==0) write(1,pos=1) Npk_tot,RTHLmax

  ! CLOSE OUTPUT FILE
  close(1)

  ! PRINT FINAL STATS
  if(myid==0) then
     write(*,101) ntile,ntasks
     write(*,102) Npk_tot
!     write(*,103) int(RTHLmin),RTHLmin-int(RTHLmin),&
!                  int(RTHLmax),RTHLmax-int(RTHLmax)
     write(*,103) RTHLmin,RTHLmax
  endif
     
  ! DATE AND TIME
  if(myid==0) then
     call date_and_time(dnt_date,dnt_time,dnt_zone,dnt_values)
     call system_clock(final_time,count_rate)

     elapsed_time    = (final_time - initial_time) / count_rate
     elapsed_hours   = int(elapsed_time / 3600)
     elapsed_minutes = int((elapsed_time - elapsed_hours * 3600)/60)
     elapsed_seconds = int(elapsed_time - elapsed_hours * 3600 &
                                        - elapsed_minutes * 60)

     write(*,401) dnt_values(3),dnt_values(2),dnt_values(1),&
       dnt_values(5),dnt_values(6),dnt_values(7),&
       dnt_zone,elapsed_hours,elapsed_minutes,elapsed_seconds

  endif

  return

101 format(3x,'ntile, ntasks:   ',2(i0,1x))
102 format(3x,'number of peaks: ',i0)
103 format(3x,'Rmin = ', f7.3,' Rmax = ',f7.3)
 12 format(/,3x,61('-'),/,3x,'Peak-patch MERGEPKVD running on',/,&
         3x,i0.2,'.',i0.2,'.',i4,1x,'at ',&
         i0.2,':',i0.2,':',i0.2,1x,'UTC',a,/,&
         3x,61('-'),/)           
401 format(/,3x,61('-'),/,3x,&
         'Peak-patch MERGEPKVD finished running on',/,&
         3x,i0.2,'.',i0.2,'.',i4,1x,'at ',&
         i0.2,':',i0.2,':',i0.2,1x,'UTC',a,/,&
         40x,'Total run time: ',i2.2,':',i2.2,':',i2.2,/&
         3x,61('-'),/)

end subroutine merge_pkvd

SUBROUTINE INDEXD3(N,A1,A2,A3,INDX)
  USE intreal_types

  !c Sort in descending order using A1,A2,A3,then INDX itself (ascend)
  !c Assumes INDX already initialized in order and points to elements 
  !c of ARRINn arrays.

  !c From Press routine INDEXX stm 16 Jun 1992 cita

  DIMENSION A1(*),A2(*),A3(*)
  INTEGER(I4B) INDX(N)

  L=N/2+1
  IR=N
10 CONTINUE
  IF(L.GT.1)THEN
     L=L-1
     INDXT=INDX(L)
  ELSE
     INDXT=INDX(IR)
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
        IF(A1(INDX(J)).GT.A1(INDX(J+1))) THEN
           J=J+1
        ELSEIF(A1(INDX(J)).EQ.A1(INDX(J+1))) THEN
           IF(A2(INDX(J)).GT.A2(INDX(J+1))) THEN
              J=J+1
           ELSEIF(A2(INDX(J)).EQ.A2(INDX(J+1))) THEN
              IF(A3(INDX(J)).GT.A3(INDX(J+1))) THEN
                 J=J+1
              ELSEIF(A3(INDX(J)).EQ.A3(INDX(J+1))) THEN
                 IF(INDX(J).LT.INDX(J+1)) J=J+1
              ENDIF
           ENDIF
        ENDIF
     ENDIF

     IF(A1(INDXT).GT.A1(INDX(J)))THEN
        INDX(I)=INDX(J)
        I=J
        J=J+J
     ELSEIF(A1(INDXT).EQ.A1(INDX(J))) THEN
        IF(A2(INDXT).GT.A2(INDX(J)) )THEN
           INDX(I)=INDX(J)
           I=J
           J=J+J
        ELSEIF(A2(INDXT).EQ.A2(INDX(J)) )THEN
           IF(A3(INDXT).GT.A3(INDX(J)))THEN
              INDX(I)=INDX(J)
              I=J
              J=J+J
           ELSEIF(A3(INDXT).EQ.A3(INDX(J))) THEN
              IF(INDXT.LT.INDX(J))THEN
                 INDX(I)=INDX(J)
                 I=J
                 J=J+J
              ELSE
                 J=IR+1
              ENDIF
           ELSE
              J=IR+1
           ENDIF
        ELSE
           J=IR+1
        ENDIF
     ELSE
        J=IR+1
     ENDIF
     GO TO 20
  ENDIF
  INDX(I)=INDXT
  GO TO 10
END SUBROUTINE INDEXD3

! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-                      

