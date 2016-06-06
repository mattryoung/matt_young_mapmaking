module grid

  use intreal_types

  IMPLICIT NONE

  real(C_FLOAT),            pointer ::  delta(:,:,:)
  real(C_FLOAT),            pointer :: deltag(:,:,:)
  complex(C_FLOAT_COMPLEX), pointer :: deltac(:,:,:)

  real, allocatable :: delta_sub(:),delta_sub_local(:)

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

