MODULE merge_subs
END MODULE merge_subs

subroutine Neighbors(ip,rnbr2)

  !c Find nearest neighbors within radius rnbr to particle ip and put in ilist.
  !c Non-wrapping version.

  !c Recursive version                 1 Mar 1992 stm cita
  !c Non-recursive version             5 Mar 1992
  !c Separate sibling cells           11 Mar 1992
  !c N-nearest nbrs                   12 Mar 1992
  !c Optimise Nfa!c                     2 Apr 1992
  !c No cloud search, just tree        3 Apr 1992
  !c No sort, from TreeSubs            5 May 1992
  !c Link-listed                      17 Jun 1992

  USE intreal_types
  !c modules needed for for Tree routines
  USE mparam
  use rpart
  USE ipart
  USE mnodes
  USE mbound
  USE mlists

  real(sp) rnbr2
  integer(i4b) ip
  integer(i4b) iTg(3),jv(3),ioff(3)
  integer(i4b) iL,jaddr,jL,jaddr2,jL2,m
  logical OctBit,found

  integer(i4b) jarry(nbitmax+1)
  !c Init counters
  nlist=0
  ilist(1)=0
  rlist(1)=0.0

  if (ip.lt.1.or.ip.gt.npmax) then
     write(*,*) 'Bad particle number ',ip
     return
  endif

  !c Get the node containing particle ip and store branches in jarry

  iL=0
  do L=1,3
     iTg(L)=ixv(L,ip)
  enddo

  call FindNodeN(iTg,iL,jarry,jL)
  jaddr=jarry(nbitmax-jL+1)

  found=.false.
  if(nocc(jaddr).ne.1) then
     write(*,*) 'Error: FindNode failed to find particle DEBUG',ip
     return
  elseif(iaddr(jaddr).eq.ip) then
     found=.true.
  endif

  nsrch=0
  !c Search for link-listed particles
  jp=ip
  do while (ixlink(jp).ne.0) 
     jpp=ixlink(jp)
     if(jpp.ne.ip) then
        rad2=(r(1,jpp)-r(1,ip))**2+(r(2,jpp)-r(2,ip))**2+(r(3,jpp)-r(3,ip))**2
        if(rad2.lt.rnbr2) then
           nsrch=nsrch+1
           nlist=nlist+1
           ilist(nlist)=jpp
           rlist(nlist)=rad2
        endif
     else
        found=.true.
     endif
     jp=jpp
  enddo

  !c Check for found particle
  if(.not.found) then
     write(*,*) 'Error: FindNode failed to find particle ',ip
     write(*,*) 'Continuing'
  endif

  !c Search back up tree for particles
  ioff(1)=0
  ioff(2)=0
  ioff(3)=0
  jL2=jL
  do while(jL2.lt.nbitmax)
     jc=jL2/3
     L=jL2-3*jc+1
     !c Get the address of current branch, then step over to sibling
     jaddr2=jarry(nbitmax-jL2+1)
     if(OctBit(jL2,iTg)) then
        jaddr2=jaddr2-1
        ioff(L)=-1
     else
        jaddr2=jaddr2+1
        ioff(L)=1
     endif
     if(nocc(jaddr2).gt.0) then
        nsrch=nsrch+nocc(jaddr2)
        call Offset(iTg,jc,ioff,jv)
        call GetNbrs(ip,iTg,jv,jaddr2,jL2,rnbr2)
     endif
     ioff(L)=0
     jL2=jL2+1
  enddo

  return
end subroutine Neighbors



subroutine FindNodeN(iBr,iL,jarry,jL)

  !c Find jaddr of target node iBr to BinLevel iL, or particle or stub above it.
  !c Node location is jaddr at level jL.

  !c Non-recursive version            5 Mar 1992 stm cita
  !c Keep parent address jarry()     11 Mar 1992
  !c revised                         16 Apr 1992

  USE intreal_types
  USE mparam
  use rpart
  USE ipart
  USE mnodes

  integer(i4b) iBr(3),iL,jarry(nbitmax+1),jL,jL2,jp,jaddr
  logical Octbit,go

  jL=nbitmax
  jaddr=1
  jarry(1)=1

  if(iL.lt.0.or.iL.gt.nbitmax) return

  go=.true.
  do while (go.and.jL.gt.iL)
     if(nocc(jaddr).gt.1) then
        !c Node - continue descent
        jL=jL-1
        jaddr=iaddr(jaddr)
        if(Octbit(jL,iBr)) jaddr=jaddr+1
        !c Add address to heap
        jarry(nbitmax-jL+1)=jaddr
     else
        !c Particle or stub - quit
        go=.false.
     endif
  enddo

  return
end subroutine FindNodeN

RECURSIVE subroutine GetNbrs(ip,iTg,jv,jaddr,jL,rnbr2)

  !c Search particles under target address jv, and place in
  !c ilist if distance from particle ip is less than rnbr.
  !c Start from location jaddr at level jL.

  !c Set up for recursive operation    4 Mar 1992 stm cita
  !c Non-recursive version             7 Mar 1992
  !c N-nearest neighbors              12 Mar 1992
  !c neighbors within rnbr2, no sort   5 May 1992
  !c added link-list                  13 Jun 1992

  USE intreal_types
  USE mparam
  use rpart
  USE ipart
  USE mnodes
  USE mbound
  USE mlists

  integer(i4b) ip,iTg(3),jv(3),jv2(3),jaddr,jL
  real(sp) rdx(3),rad2,rnbr2

  logical Octbit

  !c check to see if this cell is too far away
  if(jL.lt.0) return

  call irdist(iTg,jv,jL,rdx)
  rad2=(rdx(1)*dx(1))**2+(rdx(2)*dx(2))**2+(rdx(3)*dx(3))**2
  if(rad2.gt.rnbr2) return

  !c loop down tree, look for particles
  if(nocc(jaddr).gt.1) then
     !c node, follow branches
     jL2=jL-1
     jaddr2=iaddr(jaddr)
     jL3=jL2
     jaddr3=jaddr2+1
     do L=1,3
        jv2(L)=jv(L)
     enddo
     if(nocc(jaddr2).gt.0.and.nocc(jaddr3).gt.0) then
        !c Look at closest cell first
        L=jL2-3*(jL2/3)+1
        if(ixv(L,ip).le.jv2(L)) then
           call addbit0(jL2,jv2)
           call GetNbrs(ip,iTg,jv2,jaddr2,jL2,rnbr2)
           call addbit1(jL3,jv2)
           call GetNbrs(ip,iTg,jv2,jaddr3,jL3,rnbr2)
        else
           call addbit1(jL3,jv2)
           call GetNbrs(ip,iTg,jv2,jaddr3,jL3,rnbr2)
           call addbit0(jL2,jv2)
           call GetNbrs(ip,iTg,jv2,jaddr2,jL2,rnbr2)
        endif
     elseif(nocc(jaddr2).gt.0) then
        call addbit0(jL2,jv2)
        call GetNbrs(ip,iTg,jv2,jaddr2,jL2,rnbr2)
     elseif(nocc(jaddr3).gt.0) then
        call addbit1(jL3,jv2)
        call GetNbrs(ip,iTg,jv2,jaddr3,jL3,rnbr2)
     endif
  elseif(nocc(jaddr).eq.1) then
     !c found particle, add to list if close enough
     jp=iaddr(jaddr)
     rad2=(r(1,jp)-r(1,ip))**2+(r(2,jp)-r(2,ip))**2+(r(3,jp)-r(3,ip))**2
     if(rad2.lt.rnbr2) then
        nlist=nlist+1
        ilist(nlist)=jp
        rlist(nlist)=rad2
     endif
     !c include linked particles
     do while (ixlink(jp).ne.0) 
        jpp=ixlink(jp)
        rad2=(r(1,jpp)-r(1,ip))**2+(r(2,jpp)-r(2,ip))**2+(r(3,jpp)-r(3,ip))**2
        if(rad2.lt.rnbr2) then
           nlist=nlist+1
           ilist(nlist)=jpp
           rlist(nlist)=rad2
        endif
        jp=jpp
     enddo
  endif
  !c stub, end search

  return
end subroutine GetNbrs

subroutine BuildTree

  !c Construct Oct-tree from integer coordinates

  !c param.h : nbitmax     (parameter) max number of TOTAL bits
  !c           nbitwrd     (parameter) number of bits per word
  !c           npmax       (parameter) max number particles
  !c           nodemax     (parameter) max number of nodes

  !c parts.h : ixv(3,np)   integer particle Oct-tree locations
  !c           np          number of particles

  !c nodes.h : nocc(ntr)   occupation number of tree cell 
  !c           iaddr(ntr)  location of subnode (or particle if nocc=1)
  !c           ntr         number of tree nodes (.ge.np)
  !c           ndp         lowest depth of tree (bits)

  !c Written - recursive                25 Feb 1992 stm cita
  !c Non-recursive version               6 Mar 1992
  !c Deal with identical positions      13 Jun 1992

  USE intreal_types
  USE mparam
  USE ipart
  USE mnodes

  integer(i4b) jmax,ip,jaddr,jL,j

  write(*,*) 'beginning oct-tree build at level ',nbitmax
  !c build tree particle by particle
  ndp=nbitmax
  nplnk=0
  ntr=0
  nocc(1)=0
  iaddr(1)=0
  !c keep track of next available node storage
  jmax=1
  do ip=1,nobjs
     !c add particle to tree
     call AddPart(ip,jaddr,jL,jmax)
     if(jL.lt.ndp) ndp=jL

  enddo

  ntr=jmax
  !c report
  write(*,*) 'Oct-tree completed Nnode = ',ntr,' Depth = ',ndp
  write(*,*) 'Oct-tree Nocc(1) = ',nocc(1)
  write(*,*) 'Link-listed ',nplnk,' particles'

  return
end subroutine BuildTree



subroutine AddPart(ip,jaddr,jL,jmax)

  !c Add particle ip to the tree. Return location jaddr, level jL
  !c The maximum used node is jmax.  Non-recursive version.
  !c Checks for coincident particles if Nocc=1 (13 Jun 1992 stm)

  USE intreal_types
  USE mparam
  USE ipart
  USE mnodes

  integer(i4b) ip,jaddr,jL,jmax,jp,jaddr2
  logical Octbit,go,check

  jL=nbitmax
  jaddr=0
  if(ip.lt.1.or.ip.gt.nobjs) return

  jaddr=1
  go=.true.
  check=.true.
  do while (go.and.jL.gt.0)
     if(nocc(jaddr).gt.1) then
        !c Node - continue descent
        nocc(jaddr)=nocc(jaddr)+1
        jL=jL-1
        jaddr=iaddr(jaddr)
        if(Octbit(jL,ixv(1,ip))) jaddr=jaddr+1
     elseif(nocc(jaddr).eq.1) then
        !c Particle - check if coincident
        jp=iaddr(jaddr)
        if(check) then
           if(ixv(1,ip).eq.ixv(1,jp).and.ixv(2,ip).eq.ixv(2,jp).and.ixv(3,ip).eq.ixv(3,jp)) then
              !c Coincident - link-list new particle to old one
              jpp=jp
              do while(ixlink(jpp).ne.0) 
                 jpp=ixlink(jpp)
              enddo
              ixlink(jpp)=ip
              nplnk=nplnk+1
              go=.false.
           endif
           check=.false.
        endif
        !c Not coincident - create new nodes, put occupant ahead, continue
        if(go) then
           jL=jL-1
           jaddr2=jmax+1
           nocc(jaddr)=nocc(jaddr)+1
           iaddr(jaddr)=jaddr2
           jaddr=jaddr2
           nocc(jaddr2)=0
           iaddr(jaddr2)=0
           nocc(jaddr2+1)=0
           iaddr(jaddr2+1)=0
           jmax=jmax+2

           if(Octbit(jL,ixv(1,jp))) jaddr2=jaddr2+1
           nocc(jaddr2)=1
           iaddr(jaddr2)=jp

           if(Octbit(jL,ixv(1,ip))) jaddr=jaddr+1
        endif
     else
        !c Empty cell - put particle in and quit
        nocc(jaddr)=nocc(jaddr)+1
        iaddr(jaddr)=ip
        go=.false.
     endif
  enddo

  if(go.and.jL.eq.0) then
     !c Reached zero-bit level - link list
     if(nocc(jaddr).eq.1) then
        jpp=iaddr(jaddr)
        do while(ixlink(jpp).ne.0) 
           jpp=ixlink(jpp)
        enddo
        ixlink(jpp)=ip
        nplnk=nplnk+1
     elseif(nocc(jaddr).eq.0) then
        nocc(jaddr)=nocc(jaddr)+1
        iaddr(jaddr)=ip
     elseif(nocc(jaddr).gt.1) then
        write(*,*) 'Reach zero bit ip,nocc = ',ip,nocc(jaddr)
     endif
  endif

  return
end subroutine AddPart



subroutine Bounds(nobj,r)

  !c Find the min,max of array r
  USE intreal_types
  USE mbound

  integer(i4b) nobj,i,j
  real(sp) r(3,nobj)

  do i=1,3
     rmin(i)=r(i,1)
     rmax(i)=r(i,1)
  enddo

  do j=1,nobj
     do i=1,3
        rmin(i)=min(rmin(i),r(i,j))
        rmax(i)=max(rmax(i),r(i,j))
     enddo
  enddo

  write(*,*) 'Minimum bounds ',rmin(1),rmin(2),rmin(3)
  write(*,*) 'Maximum bounds ',rmax(1),rmax(2),rmax(3)

  return
end subroutine Bounds




subroutine XToIx

  !c Convert positions r to integer positions ixv using bounding box 
  !c rmin,rmax

  !c Scale x,y,z directions separately

  USE intreal_types
  USE mparam
  use rpart
  USE ipart
  USE mbound

  integer(i4b) i,j,jbmax,ixmin(3),ixmax(3)
  real(sp) center(3),boxlen(3),bmax

  ibmax=2**nbitwrd-1
  bmax=float(ibmax)
  jbmax=int(bmax)
  write(*,*) 'Allowed integer,real max = ',ibmax,jbmax

  if(iwrap_merge.ne.1) then
     do i=1,3
        boxlen(i)=rmax(i)-rmin(i)
        center(i)=0.5*(rmin(i)+rmax(i))
        !c Add 2% padding to box size for nonperiodic
        boxlen(i)=1.02*boxlen(i)
     enddo
  else
     write(*,*) 'Enter periodic box xc,yc,zc,len :'
     read(999,*) center(1),center(2),center(3),boxlen(1)
     boxlen(2)=boxlen(1)
     boxlen(3)=boxlen(1)
  endif

  do i=1,3
     dx(i)=boxlen(i)/bmax
     dx_1(i)=bmax/boxlen(i)
     xmin(i)=center(i)-0.5*boxlen(i)
     xmax(i)=center(i)+0.5*boxlen(i)
  enddo

  !c store integer positions and initialize link-list
  do j=1,nobj
     do i=1,3
        ixv(i,j)=int((r(i,j)-xmin(i))*dx_1(i))
     enddo
     ixlink(j)=0
  enddo

  do i=1,3
     ixmin(i)=int((rmin(i)-xmin(i))*dx_1(i))
     ixmax(i)=int((rmax(i)-xmin(i))*dx_1(i))
  enddo
  write(*,*) 'Integer min bounds ',ixmin(1),ixmin(2),ixmin(3)
  write(*,*) 'Integer max bounds ',ixmax(1),ixmax(2),ixmax(3)

  nobjs=nobj

  return
end subroutine XToIx


SUBROUTINE SORT2RI(N,RA,IB)
  USE intreal_types
  DIMENSION RA(N),IB(N)
  L=N/2+1
  IR=N
10 CONTINUE
  IF(L.GT.1)THEN
     L=L-1
     RRA=RA(L)
     IIB=IB(L)
  ELSE
     RRA=RA(IR)
     IIB=IB(IR)
     RA(IR)=RA(1)
     IB(IR)=IB(1)
     IR=IR-1
     IF(IR.EQ.1)THEN
        RA(1)=RRA
        IB(1)=IIB
        RETURN
     ENDIF
  ENDIF
  I=L
  J=L+L
20 IF(J.LE.IR)THEN
     IF(J.LT.IR)THEN
        IF(RA(J).LT.RA(J+1))J=J+1
     ENDIF
     IF(RRA.LT.RA(J))THEN
        RA(I)=RA(J)
        IB(I)=IB(J)
        I=J
        J=J+J
     ELSE
        J=IR+1
     ENDIF
     GO TO 20
  ENDIF
  RA(I)=RRA
  IB(I)=IIB
  GO TO 10
END SUBROUTINE SORT2RI
