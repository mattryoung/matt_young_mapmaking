module dens_power2
end module dens_power2

subroutine setup_dens_pspec1(gmLB,gmUB,biasnew1,biasold1)
  USE intreal_types
  use etafiles1

  real(sp) gmLB,gmUB
  write(*,*) 'READ POWER SPECTRUM 1 FROM FILE: '
  read(999,'(a)') inputfile
  write(*,'(a)') inputfile
  write(*,*) ' biasnew,biasold for power1'
  read(999,*)  biasnew1,biasold1 
  rfmin=-1
  rfmax=-1
  ifil=0
  write(*,*) ' k range in Mpc^-1 for read gmfmin,gmfmax ' 
  read(999,*) gmfmin,gmfmax
  write(*,*) ' k range in Mpc^-1 for powfns gmLB,gmUB ' 
  read(999,*) gmLB,gmUB

  call powerinit1(ifil,gmfmin,gmfmax,rfmin,rfmax,biasnew1,biasold1)

  return
end subroutine setup_dens_pspec1


subroutine setup_dens_pspec2and1(gmLB,gmUB)
  USE intreal_types
  use etafiles1
  use etafiles2
  real(sp) gmLB,gmUB
  write(*,*) 'READ POWER SPECTRUM 1 FROM FILE: '
  read(999,'(a)') inputfile
  write(*,'(a)') inputfile
  write(*,*) ' biasnew,biasold for power1'
  read(999,*)  biasnew1,biasold1 
  write(*,*) 'READ POWER SPECTRUM 2 FROM FILE: '
  read(999,'(a)') filepow2
  write(*,'(a)') filepow2
  write(*,*) ' biasnew,biasold for power2'
  read(999,*)  biasnew2,biasold2 
  rfmin=-1
  rfmax=-1
  ifil=0
  write(*,*) ' k range in Mpc^-1 for read gmfmin,gmfmax ' 
  read(999,*) gmfmin,gmfmax
  write(*,*) ' k range in Mpc^-1 for powfns gmLB,gmUB ' 
  read(999,*) gmLB,gmUB

  call powerinit1(ifil,gmfmin,gmfmax,rfmin,rfmax,biasnew1,biasold1)

  call powerinit2(ifil,gmfmin,gmfmax,rfmin,rfmax,biasnew2,biasold2)
  return
end subroutine setup_dens_pspec2and1


subroutine powerinit1()
  USE intreal_types
  use input_parameters
  use etafiles1

  !C     SETS UP TABLE FOR function power. MUST BE CALLED FIRST
  use etafiles
  use powerc1
  maxgm=4
  lnk=0
  maxgmh=maxgm/2
  biasfac=(biasold/biasnew)**2
  open(unit=30,file=trim(pkfile),status='old',form='formatted')
  m=0
  gmzp=0.0
  do i = 1, 10000
     read(30,*,end = 11) gmz,gmln,prho,prholn
     m  = m  + 1
     gmv1(m) = gmz
     prho=biasfac*prho
     if((gmz.ge.gmfmin.and.gmz.le.gmfmax).or.(gmzp.le.gmfmin.and.gmz.gt.gmfmax)) then
        if(rfmax.lt.0.0) then
            fmax=1.0
        else
           u=gmz*rfmax
           fmax=(1.-filter(u))**2
        endif
        if(rfmin.lt.0.0) then
           fmin=1.0
        else
           u=gmz*rfmin
           fmin=filter(u)**2
        endif
        px1(m) = prho*fmin*fmax
     else
        px1(m) = 0.0
     endif
     nout=200
  enddo

11 continue
  close(30)

  nptsk1 = m

  kkmax=nptsk1-maxgmh
  sump=0.0
  sumv=0.0
  do m=1,nptsk1-1
     temp=alog(gmv1(m+1)/gmv1(m))*(px1(m)+px1(m+1))
     sump=sump+temp
     sumv=sumv+temp/gmv1(m+1)/gmv1(m)
  enddo
  sump=sqrt(sump*0.5)
  sumv=sqrt(sumv*0.5)

  return
end subroutine powerinit1

subroutine renorm_powerinit2to1(biasnew1,biasold1)
  USE intreal_types
  use powerc1
  use powerc2

  nptsk1=nptsk2
  biasfac=(biasold1/biasnew1)**2
  do i=1,nptsk2
     px1(i)=px2(i)*biasfac
     gmv1(i)=gmv2(i)
  enddo
  return
end subroutine renorm_powerinit2to1

function power1(gmLB,gmUB,gm)
  USE intreal_types
  !c     this returns d sigma_rho^2 /d lnk
  use powerc1
  dimension gxa(10),pya(10)
  dimension ptab(4)
  real pint
  power1=0.0
  if(gm.gt.gmLB.and.gm.lt.gmUB) then           
     call HUNTpow(gmv1,nptsk1,gm,kk)
     if(kk.lt.2.or.kk.gt.nptsk1-2) then
        nlint=2
        kkl=kk-1
     else 
        nlint=4
        kkl=kk-2
     endif
     do ikk=1,nlint
        gxa(ikk)=gmv1(kkl+ikk)
        ptab(ikk)=px1(kkl+ikk)
     enddo
     call lagrintpow(gxa,ptab,pint,nlint,gm)
     power1=pint
  elseif(gm.eq.gmLB) then   
     power1=px1(1)
  elseif(gm.eq.gmUB) then   
     power1=px1(nptsk1)
  endif
  return
end function power1

subroutine powerinit2(ifil,gmfmin,gmfmax,rfmin,rfmax,biasnew,biasold)
  USE intreal_types

  !C     SETS UP TABLE FOR function power. MUST BE CALLED FIRST
  use etafiles2
  use powerc2
  maxgm=4
  lnk=0
  maxgmh=maxgm/2
  biasfac=(biasold/biasnew)**2
  open(unit=30,file=trim(filepow2),status='old',form='formatted')
  m=0
  gmzp=0.0
  do i = 1, 10000
     read(30,*,end = 11) gmz,gmln,prho,prholn
     m  = m  + 1
     gmv2(m) = gmz
     prho=biasfac*prho
     if((gmz.ge.gmfmin.and.gmz.le.gmfmax).or.(gmzp.le.gmfmin.and.gmz.gt.gmfmax)) then
        if(rfmax.lt.0.0) then
           fmax=1.0
        else
           u=gmz*rfmax
           fmax=(1.-filter(u,ifil))**2
        endif
        if(rfmin.lt.0.0) then
           fmin=1.0
        else
           u=gmz*rfmin
           fmin=filter(u,ifil)**2
        endif
        px2(m) = prho*fmin*fmax
     else
        px2(m) = 0.0
     endif
     nout=200
     if(m.eq.m/nout*nout) write(*,*) ' gm,px ',gmv2(m),px2(m)
  enddo
11 continue
  close(30)

  nptsk2 = m

  kkmax=nptsk2-maxgmh
  sump=0.0
  sumv=0.0
  do m=1,nptsk2-1
     temp=alog(gmv2(m+1)/gmv2(m))*(px2(m)+px2(m+1))
     sump=sump+temp
     sumv=sumv+temp/gmv2(m+1)/gmv2(m)
  enddo
  sump=sqrt(sump*0.5)
  sumv=sqrt(sumv*0.5)
  write(*,*) 'in powerinit: sump = ', sump
  write(*,*) 'in powerinit: sumv = ', sumv

  return
end subroutine powerinit2

function power2(gmLB,gmUB,gm)
  USE intreal_types
  !c    this returns d sigma_rho^2 /d lnk
  use powerc2
  dimension gxa(10),pya(10)
  dimension ptab(4)
  real(sp) pint
  power=0.0
  if(gm.gt.gmLB.and.gm.lt.gmUB) then           
     call HUNTpow(gmv2,nptsk2,gm,kk)
     if(kk.lt.2.or.kk.gt.nptsk2-2) then
        nlint=2
        kkl=kk-1
     else 
        nlint=4
        kkl=kk-2
     endif
     do ikk=1,nlint
        gxa(ikk)=gmv2(kkl+ikk)
        ptab(ikk)=px2(kkl+ikk)
     enddo
     call lagrintpow(gxa,ptab,pint,nlint,gm)
     power2=pint
  elseif(gm.eq.gmLB) then   
     power2=px2(1)
  elseif(gm.eq.gmUB) then   
     power2=px2(nptsk2)
  endif
  return
end function power2


function filter(xx)
  USE intreal_types
  use input_parameters

  filter=0.
  if(ifil.eq.1) then
     if(xx.le.1.0) filter=1.
     return
  else if(ifil.eq.0) then
     if(xx.le.5.) then
        filter=exp(-xx**2*0.5)
     endif
  else if(ifil.eq.2) then
     if(xx.le.10.) then
        s=sin(xx)
        c=cos(xx)
        filter=3.*(s-xx*c)/xx**3
     else
        filter=-3./1.414214/xx**2
     endif
  else if(ifil.eq.3) then
     !C**   THIS IS THE FFT OF THE SPH FILTER FN W_4
     !C**   gmf IS THEN h_sph
     !C***  WHY DOESN'T THE SPH FILTER FN IN PSLANG AGREE? -OTHER CASE?
     if(xx.le.0.01) then
        filter=1.0-13.0/84.0*xx**2
     else if(xx.ge.7.0) then
        !c**   SHOULD BE NEGLECTABLE SINCE R_f(Gaussian)=.56h
        !c**   sq5_2=1.5811388; THIS ASYMPTOTIC IS <filter**2>_2\pi average
        filter=-1.5811388*12.0/xx**5
     else                 
        filter=12.0/xx**5*(-2.0*sin(xx) +sin(2.0*xx)+(6.0-8.0*cos(xx) +2.0*cos(2.0*xx))/xx)
     endif
  endif
  return
end function filter

  !c**   PRESS ETAL ROUTINE ******************
SUBROUTINE HUNTpow(XX,N,X,JLO)
  USE intreal_types
  DIMENSION XX(N)
  LOGICAL ASCND
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
END SUBROUTINE HUNTpow

subroutine lagrintpow(xa,val,est,n,x)
  USE intreal_types
  dimension xa(n),val(n)
  real(sp) xa,est,val
  !C     Lagrange 4 pt interpolation or 2 point 

  if(n.eq.4) then
     a1=(x-xa(2))/(xa(1)-xa(2))*(x-xa(3))/(xa(1)-xa(3))*(x-xa(4))/(xa(1)-xa(4))
     a2=(x-xa(1))/(xa(2)-xa(1))*(x-xa(3))/(xa(2)-xa(3))*(x-xa(4))/(xa(2)-xa(4))
     a3=(x-xa(1))/(xa(3)-xa(1))*(x-xa(2))/(xa(3)-xa(2))*(x-xa(4))/(xa(3)-xa(4))
     a4=(x-xa(1))/(xa(4)-xa(1))*(x-xa(2))/(xa(4)-xa(2))*(x-xa(3))/(xa(4)-xa(3))
     est=a1*val(1)+a2*val(2)+a3*val(3)+a4*val(4)
  else 
     a1=(x-xa(2))/(xa(1)-xa(2))
     a2 = 1-a1
     est=a1*val(1)+a2*val(2)
  endif
  return
end subroutine lagrintpow

