module cloud

  use array_bounds

  integer(kind=4) ixs(np),iys(np),izs(np)
  integer(kind=4) irs2(np),ixsvec(np,3)

  equivalence (ixs,ixsvec(1,1))
  equivalence (iys,ixsvec(1,2))
  equivalence (izs,ixsvec(1,3))

end module cloud

