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

  use textlib
  use mpivars
  use exclusion

  parameter (nmax=6000000)

  real x(nmax),y(nmax),z(nmax),vx(nmax),vy(nmax),vz(nmax)
  real r(nmax)

  real xo(nmax),yo(nmax),zo(nmax),vxo(nmax),vyo(nmax),vzo(nmax)
  real ro(nmax)

  ! DOMAIN DECOMPOSITION VARIABLES
  integer, allocatable :: tilelist(:)

  real boxsize,Rbuff  
  integer itile,jtile,ktile,ntile,numtiles,tile,neach,tile_index

  integer Npk, Non2
  integer, allocatable :: Npk_eachtask(:), Npk_eachtaskl(:), Npk_begin(:)
  real    Rmin,  Rmax, Rmax_in
  real    Rminl, Rmaxl
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

  ! I/O
  character *128 filein, fileout

  ! REPORT DATE AND TIME
  if(myid==0) then
     call date_and_time(dnt_date,dnt_time,dnt_zone,dnt_values)
     write(*,12) dnt_values(3),dnt_values(2),&
          dnt_values(1),dnt_values(5),dnt_values(6),&
          dnt_values(7),dnt_zone
     call system_clock(initial_time,count_rate)
  endif

  call getarg(1,filein)
  call getarg(2,fileout)
  boxsize = r4arg(3,512.)
  ntile   = i4arg(4,8)

  numtiles  = ntile**3
  dcore_box = boxsize / ntile

  ! Random tile list
  iseed=13580
  allocate(tilelist(numtiles))
  do i=1,numtiles
    tilelist(i)=i
  enddo
  do i=numtiles,2,-1
    xcur = ran(iseed)
    j = int(xcur * i) + 1
    if(j<1.or.j>i) then
      write(*,*) 'j out of bounds',x,j
      call mpi_finalize(ierr)
      stop
    endif
    m = tilelist(j)
    tilelist(j) = tilelist(i)
    tilelist(i) = m    
  enddo

  ! Open output file 
  open(1,file=fileout,status='unknown',access='stream')
  ! -----------------------------------------------------------------------
  ! LOOP OVER TILES
  ! -----------------------------------------------------------------------
  neach  = int(numtiles/ntasks)+1
  noutl = 0
  Rmaxl = -1e10
  Rminl = 1e10

  do tile_index=myid+1,ntasks*neach,ntasks
  
  if(tile_index>numtiles) goto 201
  tile = tilelist(tile_index)

  ktile = (tile - 1 )                     / ntile**2 + 1
  jtile = (tile - ntile**2*(ktile-1) - 1) / ntile    + 1
  itile =  tile - ntile**2*(ktile-1) - (jtile-1)*ntile     

  open(4,file=filein,status='unknown',access='stream')
  read(4) Npk,Rmax_in
  Rbuff  = 2 * Rmax_in
  dL_box = dcore_box + 2 * Rbuff
  xtile = (itile - 0.5) * dcore_box - boxsize/2 
  ytile = (jtile - 0.5) * dcore_box - boxsize/2 
  ztile = (ktile - 0.5) * dcore_box - boxsize/2 

  ! DOMAIN DECOMPOSITION FILTERING OF INPUT DATA
  nin = 0
  do i=1,npk
        
     read(4) xpk,ypk,zpk,vxpk,vypk,vzpk,rpk
     fxpk = abs(2*(xpk-xtile)/dL_box)
     if(fxpk<=1) then
        
        fypk = abs(2*(ypk-ytile)/dL_box)
        if(fypk<=1) then
           
           fzpk = abs(2*(zpk-ztile)/dL_box)
           if(fzpk<=1) then
              
              nin     = nin + 1
              x(nin)  = xpk
              y(nin)  = ypk
              z(nin)  = zpk
              vx(nin) = vxpk
              vy(nin) = vypk
              vz(nin) = vzpk              
              r(nin)  = rpk 
              
           endif
        endif
     endif
     
  enddo
  
  close(4)
  if(myid==0) write(*,*) nin,'halos read...'
  
  if(nin > nmax) then
     if(myid==0) write(*,*) 'Tiles too large, exiting... nin, nmax = ',nin,nmax
     call mpi_finalize(err)
     stop
  endif

  call lagrangian_exclusion(x,y,z,vx,vy,vz,r,nin,nmerge)

  noutl=0
  do i=1,nmerge
     xpk=x(i)
     ypk=y(i)
     zpk=z(i)
     vxpk=vx(i)
     vypk=vy(i)
     vzpk=vz(i)
     rpk=r(i)     
     fxpk = abs(2*(xpk-xtile)/dcore_box)
     if(fxpk<=1) then
        
        fypk = abs(2*(ypk-ytile)/dcore_box)
        if(fypk<=1) then
           
           fzpk = abs(2*(zpk-ztile)/dcore_box)
           if(fzpk<=1) then

              noutl=noutl+1
              xo(noutl)=xpk
              yo(noutl)=ypk
              zo(noutl)=zpk
              vxo(noutl)=vxpk
              vyo(noutl)=vypk
              vzo(noutl)=vzpk
              ro(noutl)=rpk

           endif
        endif
     endif
  enddo

  write(*,*) 'myid,nin = ',nin,noutl,npk

  enddo
  ! -----------------------------------------------------------------------
  ! END LOOP OVER TILES 
  ! -----------------------------------------------------------------------

  201 continue

  call mpi_barrier(mpi_comm_world,ierr)

  allocate(Npk_eachtaskl(0:ntasks-1))
  allocate(Npk_eachtask (0:ntasks-1))
  allocate(Npk_begin    (0:ntasks-1))

  Npk_begin           = 0
  Npk_eachtaskl       = 0
  Npk_eachtaskl(myid) = noutl

  call mpi_allreduce(Npk_eachtaskl,Npk_eachtask,ntasks,mpi_integer,mpi_sum,&
                      mpi_comm_world,ierr)

  Npk_begin(0) = 0
  do i=1,ntasks-1
     Npk_begin(i)=Npk_eachtask(i-1)+Npk_begin(i-1)
  enddo

  Rminl = min(Rminl,minval(ro))
  Rmaxl = max(Rmaxl,maxval(ro))

  ! WRITE PEAKS 
  pos_offset = int(4,8)*(2+7*Npk_begin(myid))+1
  write(1,pos=pos_offset) (xo(i),yo(i),zo(i),vxo(i),vyo(i),vzo(i),&
       ro(i),i=1,noutl)
  nout = Npk_begin(ntasks-1) + Npk_eachtask(ntasks-1) 


  deallocate(Npk_eachtaskl)
  deallocate(Npk_eachtask )
  deallocate(Npk_begin    )

  call mpi_reduce(Rminl,Rmin,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
  call mpi_reduce(Rmaxl,Rmax,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)

  ! WRITE HEADER IF myid == 0       
  if(myid==0) write(1,pos=1) nout,Rmax

  ! CLOSE OUTPUT FILE
  close(1)

  ! PRINT FINAL STATS
  if(myid==0) then
     write(*,101) ntile,ntasks
     write(*,102) nout
     write(*,103) Rmin,Rmax
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
