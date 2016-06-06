program merge
  call  merge_pkvd
end program merge

!C need to modify for the no-eovlution case, ie fixed redshift in the box?? 

subroutine merge_pkvd
  USE intreal_types
!  use input_parameters
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
  parameter (nmax=30000000)

  real(sp) xon(nmax),yon(nmax),zon(nmax)
  real(sp) vxon(nmax),vyon(nmax),vzon(nmax)
  real(sp) vx2on(nmax),vy2on(nmax),vz2on(nmax)
  real(sp) strain_bar(6,nmax)
  real(sp) Fcollv(nmax),RTHLv(nmax),vTHvir(nmax)
  real(sp) hp(nmax),F_ev(nmax),F_pv(nmax)
  integer(i4b) iwas(nmax),lagrange(nmax)

  real(sp) xon2(nmax),yon2(nmax),zon2(nmax)
  real(sp) vxon2(nmax),vyon2(nmax),vzon2(nmax)
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

  character(LEN=80) seqfile,seqfile_out,inline

  integer(i4b) icindx(nclmax),nmrgv(nclmax),nbadv(nclmax)
  integer(i4b) idclus(nclmax),ic_cutv(nclmax),nicv(nclmax)
  real(sp) fcritv(nclmax),Rfclv(nclmax)
  real(sp) dLmaxv(nclmax),dFmaxv(nclmax),vtemp(nclmax),szeld(nclmax)

  logical selectRTH

  character*128 outcode

  myid=i4arg(1,1)-1
  ntasks=i4arg(2,1)

  amtorad=atan(1.0)/45.0/60.0
  onethd=1.0/3.0
  do ic=1,nclmax
     nicv(ic)=0
     nmrgv(ic)=0
     nbadv(ic)=0
     icindx(ic)=0
     idclus(ic)=0
  enddo

  read(998,'(a)') inline
  read(998,*) idovoid,ioutshear
  read(998,'(a)') inline
  read(998,*) Omx,OmB,Omhdm,fhdmclus,Omvac,Omcurv,h
  call read_cosmology(Omt,Omnr,Omx,OmB,h,Omcurv,Omvac,Omhdm,fhdmclus,iamcurved,dcurv)
  h0=h
  call Dlinear_cosmology(Omt,Omnr,Omx,OmB,h,Omcurv,Omvac,Omhdm,fhdmclus,iamcurved,dcurv)

  read(998,*) ifmt

  if(myid<10) then
     write(outcode,81) myid
  elseif(myid<100) then
     write(outcode,82) myid
  elseif(myid<1000) then
     write(outcode,83) myid
  else
     write(outcode,84) myid
  endif

81 format('000',i1)
82 format( '00',i2)
83 format(  '0',i3)
84 format(      i4)

  read(998,'(a)') seqfile_out
  if(ntasks>1) seqfile_out=seqfile_out(1:indlnb(seqfile_out))//'.'//outcode

  nic = 0
  !c Hardwire: No cluster file
  !c      call peak_init(nic,idclus,fcritv,Rfclv)

  if(nic.gt.nclmax) then
     write(*,*) '**Too many cluster types, max = ',nclmax
     write(*,*) '**Continuing with nic=0'
     nic=0
  else if(nic.gt.0) then
     idmax=0
     do ic=1,nic
        id=idclus(ic)
        icindx(id)=ic
        if(id.gt.idmax) idmax=id
     enddo
     write(*,*) 'Maximum IWAS in cfile = ',idmax
  endif

  idoindx=-1
  iLexc=-1
  iLmrg=-1

!  Lagrangian merging op-codes
!  #1 Index   | 0=None 1=RTHL 2=RTHf 3=Rf 
!  #2 Exclude | 0=None 1=(D=0) 2=(D<R1-R2) ','3=(D<R1) 4=(D<R1+R2) 
!  #3 Merge   | 0=None 1=BinaryReduceMass 
  read(998,'(a)',err=12) inline
  if(inline.ne.' ') then
     read(inline,*,err=12) iLindx,iLexc,iLmrg 
     if(iLindx.lt.0.or.iLindx.gt.3) iLindx=0
     if(iLexc.lt.0.or.iLexc.gt.4) iLexc=0
     if(iLmrg.lt.0.or.iLmrg.gt.1) iLmrg=0
  else
     iLindx=1
     iLexc=4
     iLmrg=0
     write(*,*) 'Lagrange op-code = 1 4 0 '
  endif

  iFindx=-1
  iFexc=-1
  iFmrg=-1

! Same as Lagrangian merging op-codes but for Final State Eulerian
  read(998,'(a)',err=13) inline
  if(inline.ne.' ') then
     read(inline,*,err=13) iFindx,iFexc,iFmrg
     if(iFindx.lt.0.or.iFindx.gt.4) iFindx=0
     if(iFexc.lt.0.or.iFexc.gt.1) iFexc=0
     if(iFmrg.lt.0.or.iFmrg.gt.2) iFmrg=0
  else
     iFindx=0
     iFexc=0
     iFmrg=0
     write(*,*) 'Final-state op-code = 0 0 0 '
  endif

  if(iFexc.eq.1.and.(iFmrg.eq.1.or.iFmrg.eq.2)) then
  ! Mult. factors for RTHf [exclude, merge] = [1, 2]
     read(998,'(a)',err=14) inline
     if(inline.ne.' ') then
        read(inline,*,err=14) fexc,fmrg
     else
        fexc=1.0
        fmrg=2.0
        write(*,*) 'RTHf exclude/merge factors = 1.0,2.0 '
     endif
  elseif(iFexc.eq.1.or.iFmrg.eq.1.or.iFmrg.eq.2) then
  ! Mult. factor for RTHf exclude/merge = [2] 
     read(998,'(a)',err=15) inline
     if(inline.ne.' ') then
        read(inline,*,err=15) fmrg
        fexc=fmrg
     else
        fexc=1.0
        fmrg=1.0
        write(*,*) 'RTHf exclude/merge factor = 1.0 '
     endif
  endif

  !c Hardwire: keep Lagrange
  iLagrp=0
  !c Hardwire: change vTHvir when reducing peaks on merge
  iredvth=1
  !c Hardwire: put link parent ID in Lagrange
  ikeeplnk=1
  !c Hardwire: do not apply Zeldovich at this stage
  read(998,*) iZeld
  f_disp=1.0
  if(iZeld.ne.0) then
     write(*,*) 'Multiplication factor for displacement : ', f_disp
     !cread(998,*) f_disp
  endif

  !c Hardwire: do not cut out clusters
  ncut_ic=0
  if(ncut_ic.ge.1) then
     do ic=1,ncut_ic
        write(*,*) ' cut out the next cluster type (id) '
        read(998,*) ic_cutv(ic)
     enddo
  endif

  Rmgmin=0.0
  selectRTH=.false.
  ! Min,Max RTHL for select [all]
  read(998,'(a)',err=21) inline
  if(inline.ne.' ') then
     read(inline,*,err=21) Rmgmin,Rmgmax
     if(Rmgmax.gt.Rmgmin) selectRTH=.true.
  else
     Rmgmin=0.
     Rmgmax=1.e20
     selectRTH=.false.
  endif

  ! Some Tree parameters
  iwrap_merge=0

  Non=0
  Nfile=0
  jp=0
  jf=0

  ! Read input filename
  read(998,'(a)') seqfile
  if(ntasks>1) seqfile=seqfile(1:indlnb(seqfile))//'.'//outcode

1 continue

  open(4,file=seqfile,status='old',form='unformatted')
  write(*,'Opened file (a)') seqfile
  ! New PEAKS format with shear output
  read(4) Non2,zzt,h0,omb,&
	  (xon2(i),yon2(i),zon2(i),vxon2(i),vyon2(i),vzon2(i),&
	   Fcollv2(i),RTHLv2(i),vTHvir2(i),iwas2(i),lagrange2(i),&
	   hp2(i),F_ev2(i),F_pv2(i),(strain_bar2(k,i),k=1,6),&
	   vx2on2(i),vy2on2(i),vz2on2(i),&
           i=1,Non2) 
  write(*,*) ' new Non  = ', Non2

  ! ------------------------------------------------
  ! Here we are reversing the sign of strain_bar2
  ! so strain_bar2_ij = -von2_i,j
  !
  strain_bar2 = - strain_bar2
  !
  ! Change by M.A.A. 18/06/14
  ! ------------------------------------------------

  Fmin=1.e24
  Fmax=-1.e24
  do jp2=1,Non2
     Fmax=max(Fcollv2(jp2),Fmax)
     Fmin=min(Fcollv2(jp2),Fmin)
  enddo
  write(*,*)  'Fmin,Fmax  ',Fmin,Fmax

  RTHLmin=1.e24
  RTHLmax=-1.e24
  do jp2=1,Non2
     RTHLmax=max(RTHLv2(jp2),RTHLmax)
     RTHLmin=min(RTHLv2(jp2),RTHLmin)
  enddo
  write(*,*)  'RTHLmin,RTHLmax  ',RTHLmin,RTHLmax,sum(RTHLv2)/Non2

  vTHvirmin=1.e24
  vTHvirmax=-1.e24
  do jp2=1,Non2
     vTHvirmax=max(vTHvir2(jp2),vTHvirmax)
     vTHvirmin=min(vTHvir2(jp2),vTHvirmin)
  enddo
  write(*,*)  'vTHvirmin,vTHvirmax  ',vTHvirmin,vTHvirmax

  jf=jf+1

  nrejc=0
  nmrg=0
  do jp2=1,Non2
     if(iwas2(jp2).le.-10000) then
        nrejc=nrejc+1
        goto 100
     endif
     if(ncut_ic.gt.0) then
        do ic=1,ncut_ic
           if(iwas2(jp2).eq.ic_cutv(ic)) goto 100
        enddo
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
     if(ifmt.eq.2.or.ifmt.eq.3.or.ifmt.eq.5) hp(jp)=hp2(jp2)
     if(ifmt.eq.3.or.ifmt.eq.5) then
        F_ev(jp)=F_ev2(jp2)
        F_pv(jp)=F_pv2(jp2)
     endif
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
  write(*,*) ' Completed file ',jf,' Non = ', Non
  if(nrejc.gt.0) write(*,*) ' Nbad = ',nrejc
  if(nmrg.gt.0) write(*,*) ' Nmrg = ',nmrg

2 write(*,*) 'Enter next seqfile name for input (none=end) : '
  read(998,'(a)') seqfile
  if(seqfile.ne.' ') go to 1

  !c  Files read into merged arrays

  Nfile=jf

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
  write(*,*) 'nrms ',nrms
  write(*,*) 'There are ', iquit, ' negative x peaks deleted '
  write(*,*) 'There are ', iwasbad,'negative iwas2 values'
  write(*,*) 'There are ', irthbad,'negative RTHLv values'
  vthrms=sqrt(vthrms/float(nrms))
  szzrms=sqrt(szzrms/float(nrms))
  write(*,*) 'For ',nrms,' peaks rms vTH,S = ',vthrms,szzrms

  dLmax=2.0*radLmax
  dFmax=2.0*radFmax
  write(*,*) 'Diameters dLmax,dFmax = ',dLmax,dFmax

  if(nic.gt.0) then
     write (*,*) ' id, nicv, nmrgv, RTHavg, vTHavg, Savg '
     do ic=1,iabs(nic)
        id=idclus(ic)
        if(nicv(ic).gt.0.or.nmrgv(ic).gt.0) then
           dLmaxv(ic)=dLmaxv(ic)/(nicv(ic)+nmrgv(ic))
           dFmaxv(ic)=dFmaxv(ic)/(nicv(ic)+nmrgv(ic))
           vtemp(ic)=vtemp(ic)/(nicv(ic)+nmrgv(ic))
           szeld(ic)=sqrt(szeld(ic)/(nicv(ic)+nmrgv(ic)))
        endif
        write(*,5001) id,nicv(ic),nmrgv(ic),dLmaxv(ic),vtemp(ic),szeld(ic)
5001    format(1x,i3,2x,i6,2x,i6,3(2x,f8.2))
     enddo
  endif

  !c Index peaks (Lagrangian)
  write(*,*) 
  if(iLindx.eq.3) then
     write(*,*) 'Lagrangian indexing by Rf rank '
  else if(iLindx.eq.2) then
     write(*,*) 'Lagrangian indexing by RTHf rank '
  else if(iLindx.eq.1) then
     write(*,*) 'Lagrangian indexing by RTHL rank '
  else
     write(*,*) 'Lagrangian indexing by particle number '
  endif

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
  write(*,*) 'Lagrangian clusters selected = ',nindx
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

     write(*,*) 'Lagrangian Index-ordering completed.'

  endif

  !c  Lagrangian exclusion
  if(iLexc.eq.4) then
     write(*,*) 'Lagrangian exclusion for D < R1+R2 '
  elseif(iLexc.eq.3) then
     write(*,*) 'Lagrangian exclusion for D < R1 '
  elseif(iLexc.eq.2) then
     write(*,*) 'Lagrangian exclusion for D < R1-R2 '
  elseif(iLexc.eq.1) then
     write(*,*) 'Lagrangian exclusion for D = 0 '
  else
     write(*,*) 'No Lagrangian exclusion '
  endif

  if(iLexc.ne.0) then
     write(*,*) "DEBUG: iLexc=", iLexc
     write(*,*) "DEBUG: nindx=", nindx
     write(*,*) "DEBUG: min/max of indxr:", minval(indxr),maxval(indxr)

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
     write(*,*) 'Lagrangian Oct-tree completed' 

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
     write(*,*) 'Lagrangian exclusion completed'
     write(*,*) ' number   trimmed = ',ntrim
     write(*,*) ' number remaining = ',ngood
  endif

  !c  Lagrangian merging
  write(*,*) 
  if(iLmrg.eq.1) then
     write(*,*) 'Binary Lagrangian merging with reducued mass.'
  else
     write(*,*) 'No Lagrangian merging '
  endif

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
        write(*,*) 'Lagrangian Oct-tree completed' 
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
     write(*,*) 'Lagrangian merging completed'
     write(*,*) ' number    merged = ',nmerg
     write(*,*) ' number destroyed = ',ndest
     write(*,*) ' number remaining = ',ngood
     write(*,*) ' number  unmerged = ',nunmg
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

     write(*,*) ' there are ', nbadt, ' total trimmed peaks '
     write(*,*) ' there are ', nmrgt, ' total merged peaks '
     write(*,*) ' there are ', ngoodt, ' total original peaks '

     if(nic.gt.0) then
        write (*,*) ' id, nicv, nmrgv, nbadv '
        do ic=1,iabs(nic)
           id=idclus(ic)
           write(*,*) id,nicv(ic),nmrgv(ic),nbadv(ic)
           nicv(ic)=0
           nmrgv(ic)=0
           nbadv(ic)=0
        enddo
     endif

  endif

  if(iFexc.ne.0.or.iFmrg.ne.0) then
     !c  Final-state exclusion and merging
     write(*,*)
     !c  Index peaks (Final state)
     if(iFindx.eq.4) then
        write(*,*) 'Final state indexing by current RTHL rank '
     elseif(iFindx.eq.3) then
        write(*,*) 'Final state indexing by Rf rank '
     elseif(iFindx.eq.2) then
        write(*,*) 'Final state indexing by RTHf rank '
     elseif(iFindx.eq.1) then
        write(*,*) 'Final state indexing by original RTHL rank '
     else
        write(*,*) 'Final state indexing by particle number '
     endif

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
     write(*,*) 'Final-state clusters selected = ',nindx
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

        write(*,*) 'Final-state Index-ordering completed.'
     endif
  else
     write(*,*) 
     write(*,*) 'No Final-state indexing '
  endif

  !c  Final-state (Eulerian) exclusion
  write(*,*)
  if(iFexc.eq.1) then
     write(*,*) 'Final-state exclusion for D < fR*(R1+R2) '
  else
     write(*,*) 'No Final-state exclusion '
  endif

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
     write(*,*) 'Final-state Oct-tree completed' 

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
     write(*,*) 'Final-state exclusion completed'
     write(*,*) ' number   trimmed = ',ntrim
     write(*,*) ' number remaining = ',ngood
  endif

  !c  Final-state (Eulerian) merging
  write(*,*)
  if(iFmrg.eq.2) then
     write(*,*) 'Final-state full merging '
  elseif(iFmrg.eq.1) then
     write(*,*) 'Final-state linking '
  else
     write(*,*) 'No Final-state merging '
  endif

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
        write(*,*) 'Final-state Oct-tree completed' 
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

     write(*,*) 'Final-state merging completed'
     write(*,*) ' nthrow,nlinks = ',nthrow,nlinks
     write(*,*) ' number clusters merged = ',nbadmg
     if(iFmrg.eq.2) write(*,*) ' Summed cluster masses on link'
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

     write(*,*) ' there are ', nbadt, ' trimmed peaks '
     write(*,*) ' there are ', nmrg, ' merged peaks '
     write(*,*) ' there are ', ngoodt, ' unmerged peaks '

     if(nic.gt.0) then
        write (*,*) ' id, nomrgv, nmrgv, nbadv '
        do ic=1,nic
           id=idclus(ic)
           write(*,*) id,nicv(ic),nmrgv(ic),nbadv(ic)
           nicv(ic)=0
           nmrgv(ic)=0
           nbadv(ic)=0
        enddo
     endif

  endif

  !c Write out in indexed order
  chin=6000.0/h0
  if(iZeld.ne.0) then
     write(*,*) 
     write(*,*) 'Applying Zeldovich displacements, iZeld 1 zeld, 2 2nd order Lagrangian, 3,4,5 variants '
     write(*,*) 'Factor = ',f_disp
  
     vrms=0.0
     jp2=0
     do jndx=1,nindx
        jp=indxr(jndx)
        if(ibad(jp).eq.1) goto 120
        if(iFmrg.eq.2.and.ibad(jp).eq.-1) goto 120
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
           if(fcollD<0) then
              write(*,*) 'ERROR, fcolld < 0, exiting...'
              stop
           endif
           if(iZeld.eq.2) then        
              facnl=f_disp*theta/14.0*3.0/2.0
              s2Lthv(1)=facnl*fcollD*vxon(jp)+facnl*(strain_bar(1,jp)*vxon(jp)+strain_bar(5,jp)*vzon(jp)+strain_bar(6,jp)*vyon(jp))
              s2Lthv(2)=facnl*fcollD*vyon(jp)+facnl*(strain_bar(2,jp)*vyon(jp)+strain_bar(4,jp)*vzon(jp)+strain_bar(6,jp)*vxon(jp))
              s2Lthv(3)=facnl*fcollD*vzon(jp)+facnl*(strain_bar(3,jp)*vzon(jp)+strain_bar(4,jp)*vyon(jp)+strain_bar(5,jp)*vxon(jp))
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
        vrms=vrms+((xon2(jp2)-xon(jp))*2+(yon2(jp2)-yon(jp))**2+(zon2(jp2)-zon(jp))**2)
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
     write(*,*) ' Zeldovich displacement rms = ',vrms
  else
     jp2=0
     do jndx=1,nindx
        jp=indxr(jndx)
        if(ibad(jp).eq.1) goto 123
        if(iFmrg.eq.2.and.ibad(jp).eq.-1) goto 123
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

  write(*,*) 
  write(*,*) 'Final results : '
  if(nic.gt.0) then
     write(*,*) ' id, nicv, nmrgv, nictot, nmrgtot'
     nictot=0
     nmrgtot=0
     do ic=1,nic
        id=idclus(ic)
        nictot=nictot+nicv(ic)
        nmrgtot=nmrgtot+nmrgv(ic)
        write(*,*) id,nicv(ic),nmrgv(ic),nictot,nmrgtot
     enddo
  endif
  write(*,*) 'Non2 = ', Non2

  !c Write output files 
201 continue
  if(seqfile_out.ne.' ') then
     open(4,file=seqfile_out,status='unknown',access='stream')
     if(ios.ne.0) then
        read(998,'(a)') seqfile_out
        goto 201
     endif

     if(idovoid.eq.1) then
        do i=1,Non2
           Fcollv2(i)=-Fcollv2(i)
           vxon2(i)=-vxon2(i)
           vyon2(i)=-vyon2(i)
           vzon2(i)=-vzon2(i)
           if(ioutshear.eq.1) then
              do kk=1,6
                 strain_bar2(kk,i)=-strain_bar2(kk,i)
              enddo
           endif
        enddo
     endif

     write(*,*) 'min, max RTHLv2 = ',minval(RTHLv2),maxval(RTHLv2)
     ! New PEAKS format with shear output
     write(4) Non2,zzt,h0,omB,(xon2(i),yon2(i),zon2(i),vxon2(i),vyon2(i),vzon2(i),Fcollv2(i),RTHLv2(i),& 
             vTHvir2(i),iwas2(i),lagrange2(i),hp2(i),F_ev2(i),F_pv2(i),(strain_bar2(k,i),k=1,6),i=1,Non2)
 
     close(4)

  endif

  return
end subroutine merge_pkvd

subroutine peak_init(nic,idclus,fcritv,Rfclv)
  USE intreal_types

  dimension fcritv(*),Rfclv(*),idclus(*)
  character(LEN=24) cfile
  character(LEN=70) irchar

  nic=0
  write(*,*) 'Enter cfile (.001) name for peaks [<none>]'
  read(998,'(a)') cfile
  if(cfile.eq.' ') return

  open(1,file = cfile,status = 'old',form = 'formatted')

  read(1,'(a)') irchar
  read(1,*) nic
  read(1,'(a)') irchar

  do ic=1,nic
     read(1,*) idclus(ic),fcritv(ic),Rfclv(ic)
     write(*,*) ic,idclus(ic),fcritv(ic),Rfclv(ic)
  enddo
  close(1)

  return
end subroutine peak_init


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
