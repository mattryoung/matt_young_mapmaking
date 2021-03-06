module grid

  use types

  IMPLICIT NONE

  real, allocatable :: delta(:),delta_sub(:),delta_sub_local(:)

  integer(i8b) i,j,k,n
  real dx,boxsize,local_z_origin
  integer(i8b) offset,length,index
  integer total_sub_size    

  integer(i8b) n12,n21,n2p1
  integer nsub12,nsub21,nsub2p1,nsub3
  integer ngrids  
  integer llx,lly,llz

  real xo,yo,zo
  integer gridnum 

end module grid

