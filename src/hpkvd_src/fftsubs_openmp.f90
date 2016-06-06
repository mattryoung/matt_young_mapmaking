module fftsubs

end module fftsubs
!c Press and Teulolsky article on rlft3: sent Aug 26 1989
!c Modfied to go faster using HMPC four21
!c Version 12 Feb 1993 stm      DEC Alpha compatible

subroutine rlft21(data,speq,n1,n2,n3,isign)
  USE intreal_types
  use timing_diagnostics
  real(dp) wr,wi,wpr,wpi,wtemp,theta
  integer(i4b) nn(3)
  complex c1,c2,h1,h2,w
  complex data(n1/2,n2,n3),speq(n2,n3)
  nn(1)=n1
  nn(2)=n2
  nn(3)=n3
  c1=cmplx(0.5,0.0)
  c2=cmplx(0.0,-0.5*isign)
  theta=6.28318530717959d0/dble(isign*nn(1))
  wpr=-2.0d0*dsin(0.5d0*theta)**2
  wpi=dsin(theta)
  if(isign.eq.1)then
     nn(1)=nn(1)/2
     call four21(data,nn,3,isign)
     nn(1)=nn(1)*2
     do i3=1,nn(3)
        do i2=1,nn(2)
           speq(i2,i3)=data(1,i2,i3)
        enddo
     enddo
  endif
  do i3=1,nn(3)
     j3=1
     if (i3.ne.1) j3=nn(3)-i3+2
     wr=1.0d0
     wi=0.0d0
     do i1=1,nn(1)/4+1
        j1=nn(1)/2-i1+2
        do i2=1,nn(2)
           j2=1
           if (i2.ne.1) j2=nn(2)-i2+2
           if(i1.eq.1)then
              h1=c1*(data(1,i2,i3)+conjg(speq(j2,j3)))
              h2=c2*(data(1,i2,i3)-conjg(speq(j2,j3)))
              data(1,i2,i3)=h1+h2
              speq(j2,j3)=conjg(h1-h2)
           else
              h1=c1*(data(i1,i2,i3)+conjg(data(j1,j2,j3)))
              h2=c2*(data(i1,i2,i3)-conjg(data(j1,j2,j3)))
              data(i1,i2,i3)=h1+w*h2
              data(j1,j2,j3)=conjg(h1-w*h2)
           endif
        enddo
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
        w=cmplx(sngl(wr),sngl(wi))
     enddo
  enddo
  if(isign.eq.-1)then
     nn(1)=nn(1)/2
     call four21(data,nn,3,isign)
     nn(1)=nn(1)*2
  endif
  return
end subroutine rlft21


subroutine four21(data,nn,ndim,isign)
  USE intreal_types
  use timing_diagnostics
  complex data(*),t(1048)
  integer nn(3)

  if(ndim.ne.3 .or. nn(3).gt.1024)stop 'four21 error'
  n2=nn(1)*nn(2)
  n3=n2*nn(3)
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) &
!$OMP SCHEDULE(STATIC) &
!$OMP SHARED(data,n3,n2,nn,isign)
  do i=1,n3,n2
     call fourn(data(i),nn,2,isign)
  enddo
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE) &
!$OMP SCHEDULE(STATIC) &
!$OMP SHARED(data,n2,nn,isign)
  do i=0,n2-1
     nj=-i
     do j=1,nn(3)
        nj=nj+n2
        t(j)=data(nj)
     enddo
     call four1(t,nn(3),isign)
     !c     call fourn(t,nn(3),1,isign)
     nj=-i
     do j=1,nn(3)
        nj=nj+n2
        data(nj)=t(j)
     enddo
  enddo
  return
end subroutine four21


subroutine four1(data,nn,isign)
  USE intreal_types
  real(dp) wr,wi,wpr,wpi,wtemp,theta
  dimension data(*)
  n=2*nn
  j=1
  do i=1,n,2
     if(j.gt.i)then
        tempr=data(j)
        tempi=data(j+1)
        data(j)=data(i)
        data(j+1)=data(i+1)
        data(i)=tempr
        data(i+1)=tempi
     endif
     m=n/2
1    if ((m.ge.2).and.(j.gt.m)) then
        j=j-m
        m=m/2
        go to 1
     endif
     j=j+m
  enddo
  mmax=2
2 if (n.gt.mmax) then
     istep=2*mmax
     theta=6.28318530717959d0/(isign*mmax)
     wpr=-2.d0*dsin(0.5d0*theta)**2
     wpi=dsin(theta)
     wr=1.d0
     wi=0.d0
     do m=1,mmax,2
        do i=m,n,istep
           j=i+mmax
           tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
           tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
           data(j)=data(i)-tempr
           data(j+1)=data(i+1)-tempi
           data(i)=data(i)+tempr
           data(i+1)=data(i+1)+tempi
        enddo
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
     enddo
     mmax=istep
     go to 2
  endif
  return
end subroutine four1

subroutine fourn(data,nn,ndim,isign)
  USE intreal_types
  real(dp) wr,wi,wpr,wpi,wtemp,theta
  dimension nn(ndim),data(*)
  ntot=1
  do idim=1,ndim
     ntot=ntot*nn(idim)
  enddo
  nprev=1
  do idim=1,ndim
     n=nn(idim)
     nrem=ntot/(n*nprev)
     ip1=2*nprev
     ip2=ip1*n
     ip3=ip2*nrem
     i2rev=1
     do i2=1,ip2,ip1
        if(i2.lt.i2rev)then
           do i1=i2,i2+ip1-2,2
              do i3=i1,ip3,ip2
                 i3rev=i2rev+i3-i2
                 tempr=data(i3)
                 tempi=data(i3+1)
                 data(i3)=data(i3rev)
                 data(i3+1)=data(i3rev+1)
                 data(i3rev)=tempr
                 data(i3rev+1)=tempi
              enddo
           enddo
        endif
        ibit=ip2/2
1       if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
           i2rev=i2rev-ibit
           ibit=ibit/2
           go to 1
        endif
        i2rev=i2rev+ibit
     enddo
     ifp1=ip1
2    if(ifp1.lt.ip2)then
        ifp2=2*ifp1
        theta=isign*6.28318530717959d0/(ifp2/ip1)
        wpr=-2.d0*dsin(0.5d0*theta)**2
        wpi=dsin(theta)
        wr=1.d0
        wi=0.d0
        do i3=1,ifp1,ip1
           do i1=i3,i3+ip1-2,2
              do i2=i1,ip3,ifp2
                 k1=i2
                 k2=k1+ifp1
                 tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                 tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                 data(k2)=data(k1)-tempr
                 data(k2+1)=data(k1+1)-tempi
                 data(k1)=data(k1)+tempr
                 data(k1+1)=data(k1+1)+tempi
              enddo
           enddo
           wtemp=wr
           wr=wr*wpr-wi*wpi+wr
           wi=wi*wpr+wtemp*wpi+wi
        enddo
        ifp1=ifp2
        go to 2
     endif
     nprev=n*nprev
  enddo
  return
end subroutine fourn

