module peak_gathering

  ! CLUSTER VARIABLES

  integer  Non
  integer  iwas(nmax),lagrange(nmax)

  real     xon(nmax),yon(nmax),zon(nmax)
  real     vx2on(nmax),vy2on(nmax),vz2on(nmax)
  real     vxon(nmax),vyon(nmax),vzon(nmax)
  real     strain_bar(6,nmax)
  real     Fcollv(nmax),RTHLv(nmax),vTHvir(nmax)
  real     hp(nmax),F_ev(nmax),F_pv(nmax)

end module peak_gathering

