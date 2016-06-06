module SlabToCube

  !=======================================================================
  !
  ! ROUTINES FOR INTERFACING SLAB AND CUBE DOMAIN DECOMPOSITIONS
  !
  !                                                AUTHOR: MARCELO ALVAREZ
  !                                                LAST EDIT:     02.13.15
  !
  !=======================================================================

  integer, allocatable :: s2c_tidl(:,:,:),s2c_tid(:,:,:)
  integer, allocatable :: s2c_mapl(:,:,:),s2c_map(:,:,:)

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine SlabToCube_init(nlx,nly,nlz)

    use mpivars
    use tiles
    use fftw_interface

    !---------------------------------------------------------------------
    ! IMPLICIT STATEMENT
    !---------------------------------------------------------------------

    implicit none

    integer nlx,nly,nlz

    allocate(s2c_tidl(nlx,nly,nlz))
    allocate(s2c_tid (nlx,nly,nlz))

    allocate(s2c_mapl(0:ntasks-1,0:ntasks-1,3))
    allocate(s2c_map (0:ntasks-1,0:ntasks-1,3))

    allocate(requests(0:ntasks-1))

    allocate(sendbuffer(nmesh,nmesh,local_nz))
    allocate(recvbuffer(nmesh,nmesh,nmesh))

    return

  end subroutine SlabToCube_init

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine SlabToCube_map(ix,iy,iz,nlx,nly,nlz)

    use grid
    use tiles
    use fftw_interface
    use mpivars 
    use globalvars

    !---------------------------------------------------------------------
    ! IMPLICIT STATEMENT
    !---------------------------------------------------------------------

    implicit none

    integer ix,iy,iz,nlx,nly,nlz,tid,round

    logical overlap
    
    s2c_tidl = 0
    s2c_mapl = 0

    s2c_tidl(ix,iy,iz) = myid + 1
    call mpi_allreduce(s2c_tidl,s2c_tid,nlx*nly*nlz,mpi_integer,&
         mpi_sum,mpi_comm_world,ierr)
    s2c_tid = s2c_tid - 1

    kslab1=local_z_start+1
    kslab2=local_z_start+local_nz

    round = 0
    ! LOOP OVER LAYERS (STACKED IN Z-DIRECTION)
    do ktile=1,ntile

       ! DO TILES IN THIS LAYER OVERLAP WITH SLAB?
       klayer1 = (ktile-1) * nsub + 1
       klayer2 = klayer1 + nmesh - 1
       
       ! FIND BOUNDARIES OF LAYER IN SLAB
       klins1 = max(1,klayer1-kslab1+1)
       klins2 = min(int(local_nz),klayer2-kslab1+1)

       ! FIND BOUNDARIES OF SLAB IN LAYER
       ksinl1 = max(1,kslab1-klayer1+1)
       ksinl2 = min(int(nmesh),kslab2-klayer1+1)
       nksinl = ksinl2 - ksinl1 + 1

       ! LOCAL OFFSET FOR TILES
       offset_local = ksinl1
       
       overlap = .true.
       if(klayer1>kslab2.or.klayer2<kslab1) overlap=.false.

       ! LOOP OVER TILES IN LAYER (TILED IN X-Y PLANE)
       do jtile=1,ntile
          j1=(jtile-1)*nsub+1
          do itile=1,ntile
             i1=(itile-1)*nsub+1

             tid=s2c_tid(itile,jtile,ktile)
             if(tid<0) cycle

             ! CHECK FOR OVERLAP
             if(tid>=0.and.overlap) then
                s2c_mapl(myid,tid,1) = ksinl1 + 1
                s2c_mapl(myid,tid,2) = ksinl2 + 1
                s2c_mapl(myid,tid,3) = round  + 1
             endif

             round = round + 1

          enddo
       enddo
    enddo

    call mpi_allreduce(s2c_mapl,s2c_map,ntasks*ntasks*3,mpi_integer,&
         mpi_sum,mpi_comm_world,ierr)
    s2c_map = s2c_map - 1

    return

  end subroutine SlabToCube_map

  subroutine SlabToCube_exchange(tilep,slabp)

    use grid
    use tiles
    use fftw_interface
    use mpivars 
    use globalvars

    !---------------------------------------------------------------------
    ! IMPLICIT STATEMENT
    !---------------------------------------------------------------------

    implicit none

    real(C_FLOAT), pointer :: tilep(:,:,:), slabp(:,:,:)

    integer iorig,round,tid,fid,data_count,k1,k2
    logical overlap

    integer nrecv

    kslab1=local_z_start+1
    kslab2=local_z_start+local_nz

    round = 0
    ! LOOP OVER LAYERS (STACKED IN Z-DIRECTION)
    do ktile=1,ntile

       ! DO TILES IN THIS LAYER OVERLAP WITH SLAB?
       klayer1 = (ktile-1) * nsub + 1
       klayer2 = klayer1 + nmesh - 1
       
       ! FIND BOUNDARIES OF LAYER IN SLAB
       klins1 = max(1,klayer1-kslab1+1)
       klins2 = min(int(local_nz),klayer2-kslab1+1)

       ! FIND BOUNDARIES OF SLAB IN LAYER
       ksinl1 = max(1,kslab1-klayer1+1)
       ksinl2 = min(int(nmesh),kslab2-klayer1+1)
       nksinl = ksinl2 - ksinl1 + 1

       ! LOCAL OFFSET FOR TILES
       offset_local = (ksinl1-1)*nmesh*nmesh

       ! LOOP OVER TILES IN LAYER (TILED IN X-Y PLANE)
       do jtile=1,ntile
       
          j1=(jtile-1)*nsub+1
          do itile=1,ntile
             i1=(itile-1)*nsub+1

             tid = s2c_tid(itile,jtile,ktile)
             if(tid<0) cycle

             ! CHECK FOR OVERLAP
             overlap = .true.
             if(klayer1>kslab2.or.klayer2<kslab1) overlap=.false.
             if(overlap) then

                ! PACK DATA
                
                data_count = 0
                do k=klins1,klins2
                   kk=k
                   do j=1,nmesh
                      jj=j1+j-1
                      do i=1,nmesh
                         ii=i1+i-1
                         sendbuffer(i,j,k)=slabp(ii,jj,kk)
                         data_count = data_count + 1
                      enddo
                   enddo
                enddo

                ! SEND DATA
                if(myid.ne.tid) then
                   call mpi_isend(sendbuffer(1,1,klins1),data_count,mpi_real,&
                        tid,0,mpi_comm_world,request,ierr)
                else
                   overlap=.false.
                   k1 = s2c_map(myid,myid,1)
                   k2 = s2c_map(myid,myid,2)
                   recvbuffer(:,:,k1:k2)=sendbuffer(:,:,klins1:klins2)
                endif
             endif

             ! RECEIVE DATA FOR THIS ROUND
             nrecv = 0
             do fid=0,ntasks-1
                if(s2c_map(fid,myid,3)==round.and.fid.ne.myid) then
                   nrecv=nrecv+1
                   k1 = s2c_map(fid,myid,1)
                   k2 = s2c_map(fid,myid,2)
                   data_count = (k2 - k1 + 1)*nmesh**2
                   call mpi_irecv(recvbuffer(1,1,k1),data_count,mpi_real,&
                        fid,0,mpi_comm_world,requests(fid),ierr)
                endif
             enddo

             ! THIS ROUND IS DONE ONCE ALL RECEIVING IS COMPLETE
             if(nrecv>0) then
                do fid=0,ntasks-1
                   if(s2c_map(fid,myid,3)==round.and.fid.ne.myid) &
                      call mpi_wait(requests(fid),mpi_status,ierr)
                enddo
             endif

             if(overlap) call mpi_wait(request,mpi_status,ierr)
             
             round = round + 1

          enddo
       enddo
    enddo

    do k=1,nmesh
       do j=1,nmesh
          do i=1,nmesh
             tilep(i,j,k)=recvbuffer(i,j,k)
          enddo
       enddo
    enddo

    return
    
  end subroutine SlabToCube_exchange
  
end module SlabToCube


