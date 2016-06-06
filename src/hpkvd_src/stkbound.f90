PROGRAM run_stkbound 
  call stkbound
END PROGRAM run_stkbound

SUBROUTINE stkbound 

  !     DETERMINES THE GRAPHICAL SOLUTION TO THE                          
  !     3 k BINS IF WE DO ln, hybrid and FFT sampling                     

  EXTERNAL F0 
  COMMON / clouds / lovera, icloud, rc 

  COMMON / fvals / f3, f1 
  COMMON / param / enkrad1, enkang1, enkrad2, enkang2, ensamp2, alpha2, xmax                                                      
  COMMON / gmbndry / gm0, gm_Nyq, gm_fund 
  COMMON / gmf / gmfmin, gmfmax 
  COMMON / cetaf / netaf,nkradv(3),nkangv(3),lsamlogv (3),ensampv (3),alphv(3), gmLBv(3)                                        

  !     METHOD: 1. read in number of k's ; nk+=nkrad*ifix(sqrt(nkang))**2 
  !     2. form ln(k_1/k_0)^3  using the k_1 (F_1) and k_2 bndaries (F_3) 
  !     3. form difference function F_3-F_1=F                             
  !     4. choose dx. start at x=1+eps and increment by dx until          
  !     F(x+dx)*F(x) changes sign.                                        
  !        If x exceeds sqrt3*Nyquist/k_0 no solution                     
  !     5. call RTBIS to bisect until accuracy is achieved                
  !     6. evaluate F_3 to get k_2; F_1 to get k_1                        

  read(999,*) boxsize,rc,Nfftbox,netaf
  DO ietaf = 1, netaf 
     read(999,*) nkradv(ietaf),nkangv(ietaf),lsamlogv(ietaf),ensampv(ietaf),&
               alphv(ietaf),gmLBv(ietaf)
  END DO

  lovera=Nfftbox
  dx=1e-2
  xtol=1e-5
  pi = 4.0 * atan (1.0) 
  gm0 = gmLBv(netaf) 
  alatt = 2.0 * rc / float (lovera) 
  gm_Nyq = pi / alatt 
  gm_fund = gm_Nyq * 2.0 / float (Nfftbox) 
  gm_fund_rc = gm_Nyq * 2.0 / float (lovera) 
  xmax = sqrt (3.0) * gm_Nyq / gm0 
  IF (dx.gt.xmax) dx = xmax / 3.0 
  enkang1 = ifix (sqrt (float (nkangv (2) ) ) ) **2 
  enkrad1 = nkradv (2) 
  enkang2 = ifix (sqrt (float (nkangv (1) ) ) ) **2 
  enkrad2 = nkradv (1) 
  alpha2 = alphv (1) 
  ensamp2 = ensampv (1) 

  xp = 1.00001 
  fp = F0(xp) 
  WRITE(0,*) xp,f3,f1,fp 

122 xn = xp + dx 
  fn = F0 (xn) 
  WRITE(0,*) xn,f3,f1,fn 
  IF (fp*fn.le.0.0) goto 123 
  xp = xn 
  fp = fn 
  IF (xn.gt.xmax) then 
     WRITE(0,*) ' sign of F0 never changes ' 
     STOP 
  ENDIF
  GOTO 122 
123 CONTINUE 
  xbis = RTBIS_F0 (xp, xn, xtol) 
  f31 = F0 (xbis) 
  WRITE(0,*) ' final tolerance, xbis :  ',f31,xbis 
  xm = f3 
  xmin = f1 
  gm1 = gm0 * exp (f1 / 3.0) 
  gm2 = gm1 * xbis 

  WRITE(*,6) gm1,gm2
6 FORMAT(2(1pe11.3)) 
  STOP 

END SUBROUTINE stkbound
                                                                        
                                                                        
FUNCTION F0 (x) 
  COMMON / fvals / f3, f1 
  COMMON / param / enkrad1, enkang1, enkrad2, enkang2, ensamp2, alpha2, xmax                                                      
  COMMON / gmbndry / gm0, gm_Nyq, gm_fund 
  DOUBLE PRECISION xns, xns2, alpha2_1alpha2, one_alpha2_inv, dLam1, &
       dLam1_inv, dLam2, cnst, cnstln, ensamp2invx3, ensamp2inv, pi,     &
       enA1_enA2, beta, vol2_1fn, vol2_2fn, y3ln2, y3ln1                 
  DATA lone / 1 / 
  save
  IF (lone.eq.1) then 
     one_alpha2_inv = 1.0d0 / (1.0d0 - dble (alpha2) ) 
     alpha2_1alpha2 = dble (alpha2) * one_alpha2_inv 
     dLam1 = 1.0d0 / dble (enkrad1) 
     dLam1_inv = dble (enkrad1) 
     dLam2 = 1.0d0 / dble (enkrad2) 
     ensamp2inv = 1.0d0 / dble (ensamp2) 
     ensamp2invx3 = 3.0d0 * ensamp2inv 
     pi = 4.0d0 * atan (1.0d0) 
     lone = 0 
     enA1_enA2 = dble (enkang1) / dble (enkang2) 
     cnst = (dble (gm0) / dble (gm_fund) ) **3 / dble (enkang2)  * (2.0d0 * pi) / (3.0d0)                                        
     cnstln = - log (cnst) 
     dpensamp2 = dble (ensamp2) 
     dpalpha2 = dble (alpha2) 
  ENDIF
  dpx = dble (x) 
  xns = dpx**dpensamp2 
  xns2 = xns * xns 
  beta = 0.5d0 * (1.0d0 + xns) * alpha2_1alpha2 
  vol2_1fn = (sqrt (beta**2 + one_alpha2_inv * (1.0d0 + dLam2 * (xns2 - 1.0d0) + dpalpha2 * xns) ) - beta) **ensamp2invx3 - 1.0d0 
  vol2_2fn = dpx**3 - (sqrt (beta**2 + one_alpha2_inv * (xns2 - dLam2 * (xns2 - 1.0d0) + dpalpha2 * xns) ) - beta) **ensamp2invx3 
  y3ln2 = cnstln - log (vol2_2fn) 
  y3ln1 = - dLam1_inv * log (1.0d0 - enA1_enA2 * vol2_1fn) 

  F0 = sngl (y3ln2 - y3ln1) 
  f3 = sngl (y3ln2) 
  f1 = sngl (y3ln1) 
  RETURN 
END FUNCTION F0
                                                                        
                                                                                                                                               
FUNCTION RTBIS_F0 (X1, X2, XACC) 
  PARAMETER (JMAX = 40) 
  FMID = F0 (X2) 
  F = F0 (X1) 
  IF (F * FMID.GE.0.) THEN
     WRITE(*,*) 'Root must be bracketed for bisection.' 
  ENDIF
  IF (F.LT.0.) THEN 
     RTBIS_F0 = X1 
     DX = X2 - X1 
  ELSE 
     RTBIS_F0 = X2 
     DX = X1 - X2 
  ENDIF
  DO 11 J = 1, JMAX 
     DX = DX * .5 
     XMID = RTBIS_F0 + DX 
     FMID = F0 (XMID) 
     IF (FMID.LT.0.) RTBIS_F0 = XMID 
     IF (ABS (DX) .LT.XACC.OR.FMID.EQ.0.) RETURN 
11 END DO
  WRITE(*,*) 'too many bisections' 
END FUNCTION RTBIS_F0
