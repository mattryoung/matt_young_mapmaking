!********************************************************************** 
!     BINIT.F                                                           
!     THIS PROGRAM EVALUATES Psi_{j,L},  sigma_j^2, and power spectra   
!     FOR A WIDE RANGE OF FLUCTUATION SPECTRA                           
!     ITS OUTPUT IS INPUT FOR INIT.F IN THE SPH CODE                    
!********************************************************************** 
!     SPECTRUM OUT ON 30, SPEC PARAM ON 33, XI(R) ON 31                 

PROGRAM binit 
  CALL readpar 
  CALL worker 
  CLOSE (1) 
  STOP 
END PROGRAM binit
                                                                        
!***********************************************************************
SUBROUTINE readpar 
!***********************************************************************
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)

  COMMON / warmdm / gxdec 
  COMMON / quant / gmn, fl, ns 
  COMMON / lad / lad 
  COMMON / xtrapI / lxtra, lxtrag, mount 
  COMMON / xtrap / peps1, peps2, qa1, qa2, fr 
  COMMON / omg / omg, omx, omnr, omnu, Gamma_ebw 
  
  COMMON / adisoc / brat, brat2 
  COMMON / h0 / h 
  COMMON / j3norm / rnorm, ampn, amp 
  COMMON / misc / ngm, ltrans, lpsi, lonly 
  COMMON / ifil / ifil 
  COMMON / tableI / itable, normon
  COMMON / table / amp2 
  COMMON / hbha / hbha 
  COMMON / krange / dgmlog, dgmln, gmf, gmspecl, gmspecu, gmb, gmu 
  COMMON / rcurv / rcurv, rls 
  COMMON / finI / ifluc
  COMMON / fin / fin, finp 
  COMMON / lspecI / lspec
  COMMON / lspec / pphi 
  COMMON / lspec2 / akspecmin, akspecmax, dlogakspec 
  COMMON / normI / inorm 
  COMMON / norm /  sig8in
  COMMON / holtzmann_paramsI / initialize
  COMMON / holtzmann_params / t2, t3, t4, t5, aktfnmax 
  COMMON / be92_params / bet1, bet2, bet3, anu1, anu, anu_1,akfnmax                                                           
  COMMON / BE84 / initialize_BE84 
  COMMON / ttable / initialize_table 

  REAL(dp) ns 
!  real(dp) log
  CHARACTER(LEN=100) dummy 
  CHARACTER(LEN=80) for01, fsfile, filen, filet 
                                                                        
  initialize_table = 1 
  initialize = 1 
  initialize_BE84 = 1 
  
!  WRITE ( * , * ) ' name of binit parameter file ' 
!  READ ( * , '(a)') for01 
  OPEN (unit = 1, file = 'binit.dat', status = 'old',form = 'formatted')                                               
  READ (1, 100) dummy 
  READ (1, * ) lad, itable, ns 
  WRITE ( *, 100) dummy 
  WRITE ( * , * ) 'lad,itable,ns', lad, itable, ns 
  READ (1, 100) dummy 
  WRITE ( *, 100) dummy 
  READ (1, 100) dummy 
  WRITE ( *, 100) dummy 
  READ (1, 100) dummy 
  WRITE ( *, 100) dummy 
  READ (1, * ) lxtra, lxtrag 
  WRITE ( *, * ) lxtra, lxtrag 
  READ (1, 100) dummy 
  WRITE ( *, 100) dummy 
  IF (lxtra.eq.1) then 
     READ (1, * ) peps1, peps2, ra2, ra1, fr, mount 
     qa1 = 1. / ra1 
     qa2 = 1. / ra2 
     WRITE ( *, * ) peps1, peps2, ra2, ra1, fr, mount 
  ELSE 
     READ (1, 100) dummy 
     WRITE ( *, 100) dummy 
  ENDIF
                                                                        
  READ (1, 100) dummy 
  WRITE ( *, 100) dummy 
  READ (1, * ) gxdec, brat 
  WRITE ( *, * ) gxdec, brat 
  brat2 = brat * brat 
  ibar = 1 
  READ (1, 100) dummy 
  WRITE ( *, 100) dummy 
  omnu = 0.0 
  IF (lad.eq.10.or.lad.eq.11) then 
     READ (1, * ) t2, t3, t4, t5, aktfnmax 
  ELSEIF (lad.eq.13) then 
     READ (1, * ) bet1, bet2, bet3, anu, anu1 
     akfnmax = 5.0 
     anu_1 = 1.0 / anu 
  ELSE 
     READ (1, 100) dummy 
  ENDIF
  READ (1, 100) dummy 
  WRITE ( *, 100) dummy 
  READ (1, * ) omg, omnr, omx, omnu, h 
  gmn = (omx + omnu) * h * h 
  IF (lad.eq.12.or.lad.eq.14) then 
     omb = omnr - (omx + omnu) 
     gmn = omnr * h * h * exp (- (omb * (1.0 + sqrt (2.0 * h)/ omnr) - 0.06))
     Gamma_ebw = omnr * h * exp (-(omb * (1.0 + sqrt (2.0 * h)/ omnr) - 0.06)) 
  ENDIF
  WRITE ( * , * ) 'omg,omnr,omx,omnu,h,gmn ', omg, omnr, omx, omnu, h, gmn                                                            
  READ (1, 100) dummy 
  WRITE ( *, 100) dummy 
  READ (1, * ) dgmlog 
  dgmln = dgmlog * dlog (10.d0) 
  WRITE ( * , * ) 'dgmlog,dgmln', dgmlog, dgmln 
  READ (1, 100) dummy 
  WRITE ( *, 100) dummy 
  READ (1, * ) Rmspecl, Rmspecu 
  WRITE ( * , * ) 'Rmspecl,Rmspecu', Rmspecl, Rmspecu 
  READ (1, 100) dummy 
  WRITE ( *, 100) dummy 
  READ (1, * ) Rmf, Rmb, Rmu 
  gmf = 1. / Rmf 
  gmspecl = 1. / Rmspecl 
  gmspecu = 1. / Rmspecu 
  gmb = 1. / Rmb 
  gmu = 1. / Rmu 
  WRITE ( * , * ) 'Rmf,Rmb,Rmu ', Rmf, Rmb, Rmu 
  WRITE ( * , * ) 'gmf,gmspecl,gmspecu,gmb,gmu', gmf, gmspecl, gmspecu, gmb, gmu                                                 
  ifluc = 0 
  rhor0 = 6000. / h / sqrt (omnr) 
  rhor = rhor0 
  READ (1, 100) dummy 
  WRITE ( *, 100) dummy 
  READ (1, * ) sig8in, inorm 
  WRITE (*, * ) ' sig8in,inorm= ', sig8in, inorm 
  rnorm = 10. / h 
  ampn = 0.81*sig8in**2 
  IF (inorm.eq.1) then 
     rnorm = 25. / h 
     ampn = 0.1488*sig8in**2  
  ENDIF
  WRITE ( * , * ) 'rnorm,ampn= ', rnorm, ampn 
                                                                        
  READ (1, 100) dummy 
  WRITE ( *, 100) dummy 
  READ (1, * ) ifil 
  WRITE ( *, * ) ifil 
  READ (1, 100) dummy 
  WRITE ( *, 100) dummy 
      
  READ (1, * ) hbha 
  WRITE ( *, * ) hbha	 
  READ (1, 100) dummy 
  WRITE ( *, 100) dummy 
  IF (itable.eq.1.or.itable.eq.2.or.itable.eq.3.or.itable.eq.4) then                                                              
     READ (1, '(a)') fsfile 
     WRITE ( * , '(a)') fsfile 
     OPEN (14, file = fsfile, status = 'old',form = 'formatted')                                            
     gm = 1. 
     t = tr (gm, gmn) 
     CLOSE (14) 
  ELSE 
     READ (1, 100) dummy 
     WRITE ( *, 100) dummy 
  ENDIF
  READ (1, 100) dummy 
  WRITE ( *, 100) dummy 
  READ (1, * ) lspec, akspecmin, akspecmax, dlogakspec 
  WRITE ( *, * ) lspec, akspecmin, akspecmax, dlogakspec 
  READ (1, 100) dummy 
  WRITE ( *, 100) dummy 
  IF (lspec.ne.0) then 
     READ (1, '(a)') filen 
     WRITE ( * , '(a)') filen 
     OPEN (30, file = filen, status = 'unknown', form = 'formatted')                                 
  ELSE 
     READ (1, 100) dummy 
     WRITE ( *, 100) dummy 
  ENDIF
  READ (1, 100) dummy 
  WRITE ( *, 100) dummy 
  READ (1, * ) ltrans 
  WRITE ( *, * ) ltrans 
  READ (1, 100) dummy 
  WRITE ( *, 100) dummy 
  IF (ltrans.eq.1) then 
     READ (1, '(a)') filet 
     WRITE ( * , '(a)') filet 
     OPEN (51, file = filet, status = 'unknown', form = 'formatted')                                 
  ELSE 
     READ (1, 100) dummy 
     WRITE ( *, 100) dummy 
  ENDIF
  READ (1, 100) dummy 
  WRITE ( *, 100) dummy 
  READ (1, * ) lpsi 
  WRITE ( *, * ) lpsi 
100 FORMAT(a) 
                                                                        
  RETURN 
END SUBROUTINE readpar
                                                                        
!***********************************************************************
SUBROUTINE worker 
!***********************************************************************
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)
  save
  REAL(dp) ns 
!  real(dp) log
  COMMON / ncorr / nc 
  COMMON / lspecI / lspec
  COMMON / lspec / pphi 
  COMMON / lspec2 / akspecmin, akspecmax, dlogakspec 
  COMMON / h0 / h0 
  COMMON / omg / omg, omx, omnr, omnu, Gamma_ebw 
  COMMON / j3norm / rnorm, ampn, amp 
  COMMON / quant / gmn, fl, ns 
  COMMON / ifil / ifil 
  COMMON / hbha / hbha 
  COMMON / ibar / ibar 
  COMMON / lad / lad 
  COMMON / rcurv / rcurv, rls 
  COMMON / finI / ifluc
  COMMON / fin / fin, finp 
  COMMON / krange / dgmlog, dgmln, gmf, gmspecl, gmspecu, gmb, gmu 
  DIMENSION rrv (100), flv (300) 
  DIMENSION corr (100), fint (100) 
  DIMENSION umax (0:6) 
  COMMON / misc / ngm, ltrans, lpsi, lonly 
  COMMON / normI / inorm 
  COMMON / norm / sig8in 
  DIMENSION XK ( - 4:4, 0:4) 
  CHARACTER(LEN=80) filen1 
  CHARACTER(LEN=100) dummy 
  ilog = 1 
!     ILOG=1 IS DLNK INTEGRATION; ILOG=0 IS DK INTEGRATION                
                                                                        
!     this is for normalization purposes  - we need the spectrum from   
!     0.5 Mpc to at least 1000 Mpc                                      
  umax (0) = 5.0 
  umax (1) = 1.0 
  umax (2) = 50.0 
  umax (3) = 2.0 
  umax (4) = 2.0 
  umax (5) = 25.0 
  umax (6) = 25.0 
  gmin = 1.e-04 
  gmax = 1000.0 
  ngm = idint (dlog (gmax / gmin) / dgmln) 
  h = dgmln 
  WRITE ( * , * ) ' ngm,h = ', ngm, h 
  nc = 6 
  ncor = nc 
!     don't need to calculate all of the radial ones for norm           
  amp = 1. 
  lfirst = 1 
                                                                        
!     LIGHT TRACES MASS NORMALIZATION TO J3 AT 10h^-1 MPC inorm=0       
!     LIGHT TRACES MASS NORMALIZATION TO J3 AT 25h^-1 MPC inorm=1       
!     LIGHT TRACES MASS NORMALIZATION TO DM/D AT 8h^-1 MPC inorm=2      
                                                                        
  ifil_store = 0 
  IF (ifil.eq.6) then 
     ifil_store = ifil 
     ifil = 2 
  ENDIF
  fl = 1.e-02 
  IF (ifil.eq.1) fl = 1. / gmax 
  IF (lad.eq.2) fl = 1. 
                                                                        
!     THIS IS DONE FOR ACCURACY. WE KNOW THAT THIS DOESN'T AFFECT NORM  
  r = rnorm 
  CALL integ (r, ilog, ngm, h, gmin, gmax, corr, ncor) 
  amp = ampn / corr (2) 
!     **  cannot have J_3 going negative for this normalization         
  r = 8. / h0 
  CALL integ (r, ilog, ngm, h, gmin, gmax, corr, ncor) 
  deltam = sqrt (corr (3) ) 
  WRITE ( * , * ) 'amp , deltam(wrt J_3(10)) ', amp, deltam 
  IF (inorm.eq.2) then 
     amp = amp*(sig8in /deltam)**2 
  ENDIF
37 CONTINUE 
                                                                        
  gmin = gmb 
  gmax = gmf 
  ngm = idint (dlog (gmax / gmin) / dgmln) 
  h = dgmln 
  WRITE ( * , * ) ' ngm = ', ngm 
  r = - 1. 
  fl = 1. / gmf 
  CALL integ (r, ilog, ngm, h, gmin, gmax, corr, ncor) 
                                                                        
  sig02 = corr (1) 
  sig0 = sqrt (sig02) 
  sig12 = ( - 3.) * corr (6) 
  sig1 = sqrt (sig12) 
  sig22 = 5. * corr (5) 
  sig2 = sqrt (sig22) 
  rstar = 1.73205 * sig1 / sig2 / fl 
  gamma = sig12 / sig2 / sig0 
  sigv2 = corr (4) 
  sigv = sqrt (sigv2) 
  sigm1 = sigv 
  gammav = sig02 / sigv / sig1 * h0 
  sigv = sigv * hbha 
  rcv = 1.73205 * sigm1 / sig0 
  rscoh = rstar / gamma * 1.290994 
  rk2sq = 1.732051 * rstar / gamma 
  WRITE (6, * ) 'fluc: fl,sig0,gamma,rstar,sigv,gammav ' 
  WRITE (6, 11) fl, sig0, gamma, rstar, sigv, gammav 
  WRITE (6, * ) 'fluc: gmb,gmf,rscoh,rcv,<k^2>^{-1/2} ' 
  WRITE (6, 11) gmb, gmf, rscoh, rcv, rk2sq 
                                                                        
  gmin = gmu 
  gmax = gmb 
  ngm = idint (dlog (gmax / gmin) / dgmln) 
  h = dgmln 
  WRITE ( * , * ) ' ngm = ', ngm 
                                                                        
  IF (ifil_store.eq.6) then 
     ifil = ifil_store 
  ENDIF
  IF (lspec.eq.1) then 
     gm_want = akspecmin 
     dakspec = 10.0**dlogakspec 
     gmln = dlog (gmspecu) - dgmln 
731  CONTINUE 
     gmln = gmln + dgmln 
     gm = exp (gmln) 
     CALL igrnd (r, gm, ilog, fint) 
     f4 = fint (4) 
     f4 = f4 * hbha * hbha 
     gmh = gm / h0 
     powln = dlog (fint (1) + 1.0e-33) 
!     write(6,15) gm,gmln,fint(1),powln                                 
     IF ( (gm.ge.akspecmin.and.gm.le.akspecmax) .and.gm.ge.gm_want) then                                                           
        gm_want = gm_want * dakspec 
        WRITE (30, 15) gm, gmln, fint (1), powln 
     ENDIF
     IF (gm.le.gmspecl) goto 731 
     lspec = 0 
     CLOSE (30) 
  ENDIF
                                                                        
  IF (ltrans.eq.1) then 
     gm = 1.e-05 
     tn = trans (gm) 
     gmlog = - 4. 
     dgml = 0.2 
732  CONTINUE 
     gm = 10.**gmlog 
     t = trans (gm) / tn 
     WRITE (6, * ) gm, t 
     WRITE (51, * ) gm, t 
     gmlog = gmlog + dgml 
     IF (gmlog.le.3.2) goto 732 
     CLOSE (51) 
     ltrans = 0 
  ENDIF
                                                                        
  IF (lpsi.eq.1) then 
233  CONTINUE 
     r = - 1. 
     READ (1, 100) dummy 
     WRITE ( *, 100) dummy 
     READ (1, 100) dummy 
     WRITE ( *, 100) dummy 
     READ (1, * ) lonly 
     WRITE ( *, * ) lonly 
     nc = 6 
     IF (lonly.eq.3) nc = 8 
     ncor = nc 
     READ (1, 100) dummy 
     WRITE ( *, 100) dummy 
     READ (1, * ) nfil 
     WRITE ( *, * ) nfil 
     READ (1, 100) dummy 
     WRITE ( *, 100) dummy 
     READ (1, * ) (flv (j), j = 1, nfil) 
     WRITE ( *, * ) (flv (j), j = 1, nfil) 
     READ (1, 100) dummy 
     WRITE ( *, 100) dummy 
     READ (1, '(a)') filen1 
     WRITE ( * , '(a)') filen1 
     
     IF (lonly.eq.0) then 
                                                                        
        OPEN (35, file = filen1, status = 'unknown', form = 'unformatted')                   
!     open(36,file=filen1(1:lnblnk(filen1))//'.036',                    
!     x    status='unknown',form='formatted')                           
                                                                        
     ELSEIF (lonly.eq.2.or.lonly.eq.3) then 
                                                                        
        OPEN (19, file = filen1, status = 'unknown', form = 'formatted')                     
     ENDIF
                                                                        
     DO jfil = 1, nfil 
        fl = flv (jfil) 
        IF (ifil.eq.1.and.fl.eq.0.0) fl = 1. / gmax 
        gmin = gmu 
        gmax = umax (ifil) / fl 
        ngm = idint (dlog (gmax / gmin) / dgmln) 
        h = dgmln 
        WRITE ( * , * ) ' ngm,fl,gmin,gmax = ', ngm, fl, gmin, gmax 
        CALL integ (r, ilog, ngm, h, gmin, gmax, corr, ncor) 
                                                                        
        sig02 = corr (1) 
        sig0 = sqrt (sig02) 
        sig12 = ( - 3.) * corr (6) 
        sig1 = sqrt (sig12) 
        sig22 = 5. * corr (5) 
        sig2 = sqrt (sig22) 
        rstar = 1.732050808 * sig1 / sig2 / fl 
        gamma = sig12 / sig2 / sig0 
        sigv2 = corr (4) 
        sigv = sqrt (sigv2) 
        sigm1 = sigv 
        gammav = sig02 / sigv / sig1 * h0 
        sigv = sigv * hbha 
        rcv = 1.73205 * sigm1 / sig0 
        rscoh = rstar / gamma * 1.290994 
        rk2sq = 1.732051 * rstar / gamma 
        rfilh = fl * h0 
        WRITE (6, * ) 'bgnd: rfilh,sig0,gamma,rstar,sigv,gammav ' 
        WRITE (6, 11) rfilh, sig0, gamma, rstar, sigv, gammav 
        WRITE (6, * ) 'bgnd: gmb,gmf,rscoh,rcv,<k^2>^{-1/2} ' 
        WRITE (6, 11) gmb, gmf, rscoh, rcv, rk2sq 
        IF (lonly.eq.2) then 
                                                                        
           WRITE (19, 11) rfilh, sig0, gamma, rstar, sigv, gammav 
        ELSEIF (lonly.eq.3) then 
           rfilh = fl * h0 
           amu = corr (7) / corr (1) 
           dsig_dlnr = sqrt (corr (8) - amu * amu * corr (1) ) 
           WRITE (19, 11) rfilh, sig0, gamma, rstar, amu, dsig_dlnr 
        ENDIF
     END DO
     nc = 13 
     ncor = nc 
                                                                        
     IF (lonly.eq.0) then 
        READ (1, 100) dummy 
        WRITE ( *, 100) dummy 
        READ (1, * ) nrrv 
        WRITE ( *, * ) nrrv 
        READ (1, 100) dummy 
        WRITE ( *, 100) dummy 
        READ (1, * ) (rrv (kk), kk = 1, nrrv) 
        WRITE ( *, * ) (rrv (kk), kk = 1, nrrv) 
        DO kk = 1, nrrv 
           r = rrv (kk) 
           CALL integ (r, ilog, ngm, h, gmin, gmax, corr, ncor) 
           DO jej = - 4, 4 
              DO lel = 0, 4 
                 XK (jej, lel) = 0.0 
              END DO
           END DO
           
           XK (0, 0) = corr (1) 
           XK ( - 1, 1) = r * corr (2) / 3. 
           XK (2, 0) = corr (4) 
           XK (4, 0) = corr (5) 
           XK (4, 1) = corr (6) 
           XK (2, 2) = corr (7) 
           XK (4, 2) = corr (8) 
           XK (4, 3) = corr (9) 
           XK (4, 4) = corr (10) 
           XK ( - 1, 3) = corr (11) 
           XK ( - 2, 2) = corr (12) 
           XK ( - 2, 0) = corr (13) 
                                                                        
           xi = corr (1) 
           psi = xi / sig02 
           psi00 = psi 
!     psi20=corr(4)/sig12                                               
!     psi40=corr(5)/sig22                                               
!     psi41=corr(6)/sig22                                               
!     psi22=corr(7)/sig12                                               
!     psi42=corr(8)/sig22                                               
!     psi43=corr(9)/sig22                                               
!     psi44=corr(10)/sig22                                              
           rh = r * h0 
           IF (r.eq.0.) then 
              rlog = - 25 
           ELSE 
              rlog = dlog10 (r) 
           ENDIF
!     NOTATION IS AS IN BLM NOTES i.e. Psi_{j,L}                        
!     write(6,11) r,rlog,XK(-2,0),XK(4,1),XK(0,0),XK(2,0)               
           WRITE (35) r, rlog, ( (XK (jej, lel), jej = - 4, 4), lel = 0, 4)                                              
!     psi00,psi20,psi40,psi22                                           
!     write(6,11) r,rlog,psi41,psi42,psi43,psi44                        
!     write(36,11) r,rlog,psi41,psi42,psi43,psi44                       
                                                                        
        END DO
                                                                        
        CLOSE (35) 
        CLOSE (36) 
     ELSE 
        READ (1, 100) dummy 
        WRITE ( *, 100) dummy 
        READ (1, 100) dummy 
        WRITE ( *, 100) dummy 
        READ (1, 100) dummy 
        WRITE ( *, 100) dummy 
        READ (1, 100) dummy 
        WRITE ( *, 100) dummy 
     ENDIF
     CONTINUE 
     READ (1, 100) dummy 
     WRITE ( *, 100) dummy 
     READ (1, * ) lyes 
     WRITE ( *, * ) lyes 
     READ (1, 100) dummy 
     WRITE ( *, 100) dummy 
     IF (lyes.eq.1) goto 233 
  ENDIF
11 FORMAT(1x,7(1pe11.3)) 
15 FORMAT(6(1pe12.4)) 
100 FORMAT(a) 
                                                                        
  RETURN 
END SUBROUTINE worker
                                                                        
!***********************************************************            
SUBROUTINE integ (r, ilog, ngm, h, gmin, gmax, corr, ncor) 
!***********************************************************            
!     SIMPSON RULE INTEGRATION ROUTINE                                  
!     THIS IS SET UP ONLY IN THE ILOG=1 MODE                             
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)
  save
  DIMENSION corr (100), fint (100) 
  COMMON / nc / nc 
                                                                        
  nc = ncor 
  e = exp (h) 
  eh = sqrt (e) 
  h3 = h / 3. 
  gm = gmin 
  CALL igrnd (r, gm, ilog, fint) 
     write(*,*) gm,fint(2)                                             
  DO i = 1, ncor 
     corr (i) = fint (i) / 2. 
  END DO
  gm = gm * eh 
  CALL igrnd (r, gm, ilog, fint) 
  DO i = 1, ncor 
     corr (i) = corr (i) + fint (i) * 2. 
  END DO
  ngm1 = ngm - 1 
  DO j = 1, ngm1 
     gm = gm * eh 
     CALL igrnd (r, gm, ilog, fint) 
     DO i = 1, ncor 
        corr (i) = corr (i) + fint (i) 
     END DO
     gm = gm * eh 
     CALL igrnd (r, gm, ilog, fint) 
     DO i = 1, ncor 
        corr (i) = corr (i) + fint (i) * 2. 
     END DO
  END DO
  gm = gm * eh 
  CALL igrnd (r, gm, ilog, fint) 
  DO i = 1, ncor 
     corr (i) = h3 * (corr (i) + fint (i) / 2.) 
  END DO
  RETURN 
END SUBROUTINE integ
                                                                        
!*******************************************************************    
SUBROUTINE igrnd (r, gm, ilog, fint) 
!*******************************************************************    
                                                                        
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)
  PARAMETER (pi2 = 19.73920) 
  REAL(dp) ns 
  save
  DIMENSION win (100), fint (100) 
  COMMON / nc / ncor 
  COMMON / xtrapI / lxtra, lxtrag, mount 
  COMMON / xtrap / peps1, peps2, qa1, qa2, fr 
  COMMON / omg / omg, omx, omnr, omnu, Gamma_ebw 
  COMMON / h0 / h 
  COMMON / quant / gmn, fl, ns 
  COMMON / lspecI / lspec
  COMMON / lspec / pphi 
  COMMON / j3norm / rnorm, ampn, amp 
  COMMON / misc / ngm, ltrans, lpsi, lonly 
  eps = 0. 
  qc = 0. 
  tt = trans (gm) 
  del2 = gm** (2 + ilog + ns) * tt / pi2 
  IF (lxtra.eq.1.and.lxtrag.eq.0) then 
     IF (gm.lt.qa2) then 
        eps = 1. 
        IF (gm.le.qa1) then 
           qc = (qa2 / qa1) **peps2 - 1. 
           IF (mount.eq.1) eps = (gm / qa1) **peps1 
        ELSE 
           qc = (qa2 / gm) **peps2 - 1. 
        ENDIF
     ENDIF
     del2 = fr * eps * del2 * qc + del2 
  ENDIF
  IF (lxtra.eq.1.and.lxtrag.eq.1) then 
     IF (gm.lt.4. * qa2) then 
        qc2 = (gm / qa2) **2 
        e2 = exp ( - qc2) 
        eps = e2 
        qc = qc2 
        IF (gm.le.4. * qa1) then 
           qc1 = (gm / qa1) **2 
           e1 = exp ( - qc1) 
           IF (mount.eq.1) then 
              eps = (1. - e1) * e2 
           ELSE 
              eps = (1. - e1) **peps1 * e2 
           ENDIF
        ENDIF
     ENDIF
     del2 = fr * eps * del2 / qc + del2 
  ENDIF
  delout2 = amp * del2 
  fil = filter (gm) 
  del2fil = delout2 * fil * fil 
  IF (r.lt.0) then 
     CALL wind (gm, win) 
  ELSE 
     CALL windr (gm, r, win) 
  ENDIF
  IF (lpsi.eq.1.and.lonly.eq.3) ncor = ncor - 2 
  DO i = 1, ncor 
     fint (i) = del2fil * win (i) 
  END DO
                                                                        
  IF (lpsi.eq.1.and.lonly.eq.3) then 
     ncor = ncor + 2 
     dfil = fil_dlnr (gm) 
     fint (ncor - 1) = delout2 * dfil * fil * win (1) 
     fint (ncor) = delout2 * dfil * dfil * win (1) 
  ENDIF
  IF (lspec.ne.0) then 
     pphi = fint (1) * (1.5 * omnr * h * h / (3000. * gm) **2) **2 
!     if(lspec.eq.1.and.tt.ne.0.) then                                  
!     pphi=pphi/tt                                                      
!     endif                                                             
  ENDIF
  RETURN 
END SUBROUTINE igrnd
                                                                        
!********************************************************************   
FUNCTION trans (gm) 
!********************************************************************   
                                                                        
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)
  REAL(dp) ns 
  save
  COMMON / adisoc / brat, brat2 
  COMMON / quant / gmn, fl, ns 
  COMMON / ibar / ibar 
  COMMON / tableI / itable, normon
  COMMON / table / amp2 
  COMMON / lad / lad 
  
  IF (lad.eq. - 1) then 
     trans = 1.0 
  ELSEIF (lad.eq.1) then 
     trans = tad (gm, gmn) **2 
  ELSEIF (lad.eq.12) then 
     trans = tad_ebw92 (gm, gmn) **2 
  ELSEIF (lad.eq.14) then 
     Dnu_choice = Dnu (gm, Dnu_cdm, Dcdm_cdm) 
     trans = tad_ebw92 (gm, gmn) **2 * Dnu_choice**2 
  ELSEIF (lad.eq.8) then 
     trans = tadBE (gm, gmn) **2 
  ELSEIF (lad.eq.9) then 
     trans = tad_be83 (gm, gmn) **2 
  ELSEIF (lad.eq.10) then 
     trans = tfn_holtzmann (gm) **2 
  ELSEIF (lad.eq.11) then 
     trans = tfn_holtzmann_steep (gm) **2 
  ELSEIF (lad.eq.13) then 
     trans = tad_be92 (gm) **2 
  ELSEIF (lad.eq.14) then 
     trans = tadBE (gm, gmn) **2 * Dnu (gm, Dnu_cdm, Dcdm_cdm) **2 
  ELSEIF (lad.eq.0.or.lad.eq.2) then 
     trans = tisoc (gm, gmn) **2 
  ELSEIF (lad.eq.4) then 
     trans = (tad (gm, gmn) * fnu (gm, gmn) ) **2 
  ELSEIF (lad.eq.5) then 
     trans = (tad (gm, gmn) * fwarm (gm, gmn) ) **2 
  ELSEIF (lad.eq.6) then 
     trans = ( (tisoc (gm, gmn) * brat) **2 + tad (gm, gmn) **2) 
  ELSEIF (lad.eq.7) then 
     trans = 1. 
  ENDIF
  IF (itable.eq.1.or.itable.eq.2.or.itable.eq.3.or.itable.eq.4) then                                                              
     trans = tr (gm, gmn) **2 
  ENDIF
  IF (ibar.eq.1) then 
     tb = tbar (gm, gmn) 
     trans = trans / tb / tb 
  ENDIF
  RETURN 
END FUNCTION trans
                                                                        
!******************************************************************     
FUNCTION tbar (gm, gmn) 
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)

  tbar = 1. 
  IF (gm.ge.10.) tbar = 1. + (gm * 1.2e-03) **2 / gmn 
  RETURN 
END FUNCTION tbar
                                                                        
!*******************************************************************    
FUNCTION tad (gm, gmn) 
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)

  q = gm / gmn 
  y = 1. + 3.89 * q + (16.1 * q) **2 + (5.46 * q) **3 + (6.71 * q)**4                                                               
  g = 1. 
  IF (q.gt.0.001) g = dlog (1. + 2.34 * q) / 2.34 / q 
  y = g / y**0.25 
  tad = y 
  RETURN 
END FUNCTION tad
                                                                        
!*******************************************************************    
FUNCTION tadBE (gm, gmn) 
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)
  save
  PARAMETER (NBE = 8) 
  COMMON / omg / omg, omx, omnr, omnu, Gamma_ebw 
  COMMON / h0 / h 
                
  REAL(dp), dimension(NBE), PARAMETER :: av = (/ 11.3, 23.1, 23.5, 37.1, 48.2, 76.2, 37.9, 207.0 /)
  REAL(dp), dimension(NBE), PARAMETER :: bv = (/ 5.29, 11.4, 17.3, 21.1, 31.6, 71.1, 18.4, 47.9 /) 
  REAL(dp), dimension(NBE), PARAMETER :: cv= (/ 3.10, 6.48, 7.15, 10.8, 15.6, 27.7, 10.1, 38.5 /) 
  REAL(dp), dimension(NBE), PARAMETER :: pnuv = (/ 1.13, 1.25, 1.07, 1.12, 1.30, 2.15, 1.21, 1.73 /) 
  REAL(dp), dimension(NBE), PARAMETER :: omtv = (/ 1.0, 1.0, 0.4, 0.3, 0.2, 0.2, 0.2, 0.2 /) 
  REAL(dp), dimension(NBE), PARAMETER :: ombv = (/ 0.03, 0.03, 0.03, 0.03, 0.03, 0.10, 0.03, 0.03 /) 
  REAL(dp), dimension(NBE), PARAMETER :: h0v = (/ 0.75, 0.50, 0.75, 0.75, 0.75, 0.75, 1.00, 0.50 /) 
  DATA lone / 1 / 
                                      
  IF (lone.eq.1) then 
!     write(*,*) ' which BE84 case ? menu: omt, omb, h0 '               
     omb = omg - omx 
     menu = 0 
     DO jj = 1, NBE 
        IF ((abs(h-h0v(jj)).lt.1.e-03).and.(abs(omg-omtv(jj)).lt.1.e-03) &
             .and.(abs(omb-ombv(jj)).lt.1.e-03) ) then
           menu = jj 
        ENDIF
!     write(*,*) jj, ': ',omtv(jj),ombv(jj),h0v(jj)                     
                                                                        
     END DO
     IF (menu.eq.0) then 
        WRITE ( * , * ) ' PARAMETERS WRONG FOR BE84 CASE ' 
        STOP 
     ENDIF
!     read(999,*) menu                                                    
!     write(*,*) ' BE84 menu choice = ', menu                           
     a = av (menu) 
     b = bv (menu) 
     c = cv (menu) 
     pnu = pnuv (menu) 
     pnu_1 = 1.0 / pnu 
     lone = 0 
  ENDIF
!      y=(1.0+(a*gm+(b*gm)**1.5+(c*gm)**2)**pnu)**pnu_1                 
!      tadBE=1.0/y                                                      
  tadBE = a * gm + (b * gm) **1.5 + (c * gm) **2 
  tadBE = tadBE**pnu 
  tadBE = 1.0 / (1.0 + tadBE) **pnu_1 
  RETURN 
END FUNCTION tadBE
                                                                        
!*******************************************************************    
FUNCTION tad_ebw92 (gm, gmn) 
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)
  save
  PARAMETER (a = 25.6, b = 12.0, c = 6.8, pnu = 1.13, pnu_1 = 1.0 / pnu)                                                              
  COMMON / omg / omg, omx, omnr, omnu, Gamma_ebw 
  q = 0.25 * gm / gmn 
!      y=(1.0+(a*q+(b*q)**1.5+(c*q)**2)**pnu)**pnu_1                    
!      tad_ebw92=1.0/y                                                  
  tad_ebw92 = a * q + (b * q) **1.5 + (c * q) **2 
  tad_ebw92 = tad_ebw92**pnu 
  tad_ebw92 = 1.0 / (1.0 + tad_ebw92) **pnu_1 
  RETURN 
END FUNCTION tad_ebw92
                                                                        
!*******************************************************************    
FUNCTION Dnu (gm, Dnu_cdm, Dcdm_cdm) 
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)
  save
  PARAMETER (c1 = 0.0015, c2 = -0.1207, c3 = 0.1015, c4 = -0.01618, c5 = 0.001711)                                           
  PARAMETER (b1 = 0.01647, b2 = 2.803e-05, b3 = 10.90, b4 = 3.259) 
  COMMON / omg / omg, omx, omnr, omnu, Gamma_ebw 
  COMMON / h0 / h 
  COMMON / redshift / ZZ 
  omgamma = (2.46e-5) / .25 
  beta = (5 - sqrt (25 - 24 * omnu) ) / 4 
  tgamma = 0.9967 
  aratio = 1 / (1.681 * omgamma) 
  Rnr = 10.8 * tgamma * tgamma 
  B = Rnr * (1 + .1435) / (omnu + .1435) 
  Rstar = 69.52 * tgamma * tgamma 
  tmp = sqrt (omnu * (1 - .9465 * omnu) ) / (1 + (9.259 * omnu) **2) 
  A = Rstar * tmp * (1 + 10.912 * omnu) 
  ak = gm 
  tm1 = 1 + (A * ak) **2 + ( (1 / aratio) * (B * ak) **4) * (1 - omnu) ** (1 / beta) 	                                             
  tm2 = 1 + (B * ak) **2 - (B * ak) **3 + (B * ak) **4 
  Dnu = (tm1 / tm2) **beta 
  Gamma_nu = omnu * h * h 
  ZZ = 1.0 
  x = gm / Gamma_nu / h / sqrt (ZZ) 
  tm1 = exp ( - c1 * x) 
  tm2 = 1 + c2 * sqrt (x) + c3 * x + c4 * x * sqrt (x) + c5 * x * x 
  Dnu_cdm = sqrt (tm1 / tm2) 
  x0 = x 
  tm1 = 1.0 + b1 * x** (b4 / 2) + b2 * x**b4 
  tm2 = 1 + b3 * x0**b4 
  Dcdm_cdm = sqrt ( (tm1 / tm2) ** (omnu**1.05) ) 
  RETURN 
END FUNCTION Dnu
                                                                        
!*******************************************************************    
FUNCTION tad_be83 (gm, gmn) 
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)
  save
  q = gm / gmn 
  y = 1. + 1.7 * q + 9.0 * q**1.5 + q**2 
  tad_be83 = 1.0 / y 
  RETURN 
END FUNCTION tad_be83
                                                                        
!*******************************************************************    
FUNCTION fnu (gm, gmn) 
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)
  save
 
  y = 0. 
  z = gm * 2.6 / gmn 
  IF (z.le.10.) then 
     y = exp ( - 0.16 * z - 0.5 * z**2) 
  ENDIF
  fnu = y 
  RETURN 
END FUNCTION fnu
                                                                        
!*******************************************************************    
FUNCTION fwarm (gm, gmn) 
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)
  save
   
  COMMON / warmdm / gxdec 
  DATA lone / 1 / 
  IF (lone.eq.1) then 
                                                                        
     rfw = 0.2 / gmn * (100. / gxdec) **1.33333 
     lone = 0 
  ENDIF
  y = 0. 
  z = gm * rfw 
  IF (z.le.10.) then 
     y = exp ( - 0.5 * z - 0.5 * z**2) 
  ENDIF
  fwarm = y 
  RETURN 
END FUNCTION fwarm
                                                                        
                                                                        
!*******************************************************************    
FUNCTION tisoc (gm, gmn) 
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)
  save
   
  q = gm / gmn 
  y = (40. * q) **2 / (1. + 215. * q + (16. * q) **2 / (1. + 0.5 * q) )                                                              
  y = 1. + y + (5.6 * q) **1.6 
  y = 1. / y**1.25 
!     EB formula                                                        
!d    y=15.*q+(0.86*q)**1.5 +(5.6*q)**2                                 
!d    y=1.+y**1.24                                                      
!d    y=1./y**0.807                                                     
                                                                        
  tisoc = y 
  RETURN 
END FUNCTION tisoc
                                                                        
                                                                        
!*************************************************************          
                                                                        
FUNCTION tr (gm, gmn) 
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)
  save
!  real(dp) log
                                                                        
  COMMON / lad / lad 
  COMMON / out / tx (10000), gmv (10000), gmlnv (10000) 
  COMMON / tableI / itable, normon
  COMMON / table / amp2 
  
  COMMON / omg / omg, omx, omnr, omnu, Gamma_ebw 
  COMMON / h0 / h 
  COMMON / ttable / initialize_table 
  DIMENSION txis (5000) 
                                                                        
  DATA lnint / 1 / 
                                                                        
  IF (initialize_table.eq.1) then 
                                                                        
     m = 0 
     IF (itable.eq.2) then 
        READ (14, * ) tauf1 
        chsi = (1. - omg) * h**2 / .312 / 57300. 
        IF (omg.lt.1.) tauf1 = sinh (sqrt (chsi) * tauf1) / sqrt (chsi)                                                       
        arg = 1.0 / 6.581 
     ELSEIF (itable.eq.3) then 
        READ (14, * ) ixon, omg, omb, h, tauf1 
        chsi = (1. - omg) * h**2 / .312 / 57300. 
        IF (omg.lt.1.) tauf1 = sinh (sqrt (chsi) * tauf1) / sqrt (chsi)                                                       
! TO BE FIXED FOR MORE CORRECT NUMBERS                                  
        arg = 1.0 / 6.581 
     ENDIF
     DO i = 1, 100000 
        IF (itable.eq.1) then 
           READ (14, *, end = 10) gmz, deltax 
        ELSEIF (itable.eq.2) then 
           READ (14, *, end = 10) gmz, D2LZ, drmsz, drmspz, deltz, deltax                                                   
           gmz = gmz * arg 
        ELSEIF (itable.eq.3) then 
           READ (14, *, end = 10) gmz, D2LZ, drmsz, drmspz, deltab, deltax, zz3, zz2, zz5                                    
           gmz = gmz * arg 
        ELSEIF (itable.eq.4) then 
           READ (14, *, end = 10) gmz, deltax, deltab, deltag, deltanu, deltamnu                                        
           gmz = gmz * h 
        ELSEIF (itable.eq.5) then 
           READ (14, *, end = 10) gmz, deltax, deltab, deltag, deltanu, deltamnu, deltatot                                        
           gmz = gmz * h 
        ENDIF
        m = m + 1 
        IF (lad.eq.0.or.lad.eq.2) then 
           txis (m) = deltax 
           deltax = deltax / gmz**2 
        ENDIF
        gmv (m) = gmz 
        gmlnv (m) = dlog (gmz) 
                                                                        
        tx (m) = deltax 
     END DO
!                                                                       
10   CONTINUE 
     npts = m 
     igp = 1 
     gfirst = gmv (1) 
     glast = gmv (npts) 
     tfirst = tx (1) 
     tfiso = txis (1) 
                                                                        
     gmnf = omnr * h * h 
     IF (lad.eq.1.or.lad.eq.3) then 
        tf = tad (gfirst, gmnf) 
        tl = tad (glast, gmn) 
        amp = tf / tfirst 
     ELSEIF (lad.eq.12) then 
        tf = tad_ebw92 (gfirst, gmnf) 
        tl = tad_ebw92 (glast, gmn) 
        amp = tf / tfirst 
     ELSEIF (lad.eq.14) then 
        Dnu_choice = Dnu (gfirst, Dnu_cdm, Dcdm_cdm) 
        tf = tad_ebw92 (gfirst, gmnf) * Dnu_choice 
        Dnu_choice = Dnu (glast, Dnu_cdm, Dcdm_cdm) 
        tl = tad_ebw92 (glast, gmn) * Dnu_choice 
        amp = tf / tfirst 
     ELSEIF (lad.eq.13) then 
        tf = tad_be92 (gfirst) 
        tl = tad_be92 (glast) 
        amp = tf / tfirst 
     ELSE 
        tf = tisoc (gfirst, gmnf) 
        tl = tisoc (glast, gmn) 
        amp = tf / tfiso 
     ENDIF
                                                                        
     IF(lad.eq.1.or.lad.eq.3.or.lad.eq.12.or.lad.eq.13.or.lad.eq.14) then                                                           
        DO jp = 1, npts 
           tx (jp) = tx (jp) * amp 
        enddo
     ELSE 
        DO jp = 1, npts 
           txis (jp) = txis (jp) * amp 
        enddo
     ENDIF
     tfirst = tx (1) 
     tfiso = txis (1) 
     tlast = tx (npts) 
     tliso = txis (npts) 
     tr = 1.0 
                                                                        
     initialize_table = 0 
     RETURN 
  ENDIF
  IF (lad.eq. - 1) then 
     tr = 1.0 
     RETURN 
  ENDIF
  IF (gm.le.gfirst) then 
     IF (lad.eq.1.or.lad.eq.3) then 
        t = tad (gm, gmnf) 
     ELSEIF (lad.eq.12) then 
        t = tad_ebw92 (gm, gmnf) 
     ELSEIF (lad.eq.14) then 
        Dnu_choice = Dnu (gm, Dnu_cdm, Dcdm_cdm) 
        t = tad_ebw92 (gm, gmnf) * Dnu_choice 
     ELSEIF (lad.eq.13) then 
        t = tad_be92 (gm) 
     ELSE 
        t = tisoc (gm, gmnf) 
     ENDIF
     tr = tfirst * t / tf 
     ig = 1 
     RETURN 
  ENDIF
                                                                        
  IF (gm.ge.glast) then 
     IF (lad.eq.1.or.lad.eq.3) then 
        t = tad (gm, gmn) 
     ELSEIF (lad.eq.12) then 
        t = tad_ebw92 (gm, gmn) 
     ELSEIF (lad.eq.14) then 
        Dnu_choice = Dnu (gm, Dnu_cdm, Dcdm_cdm) 
        t = tad_ebw92 (gm, gmn) * Dnu_choice 
     ELSEIF (lad.eq.13) then 
        t = tad_be92 (gm) 
     ELSE 
        t = tisoc (gm, gmn) 
     ENDIF
     tr = tlast * t / tl 
     IF (lad.eq.2) tr = tliso / gm**2 
     RETURN 
  ENDIF
  ig = igp 
11 CONTINUE 
  ig1 = ig + 1 
  IF (gmv (ig1) .gt.gm.and.gmv (ig) .le.gm) then 
     IF (lnint.eq.1) then 
        gmln = dlog (gm) 
        gg = (gmln - gmlnv (ig) ) / (gmlnv (ig1) - gmlnv (ig) ) 
     ELSE 
        gg = (gm - gmv (ig) ) / (gmv (ig1) - gmv (ig) ) 
     ENDIF
     IF (lad.ne.2) then 
        tr = tx (ig) + ( (tx (ig1) - tx (ig) ) * gg) 
     ELSE 
        tr = txis (ig) + ( (txis (ig1) - txis (ig) ) * gg) 
        tr = tr / gm**2 
     ENDIF
     GOTO 12 
  ENDIF
  ig = ig + 1 
  IF (ig.eq.npts) ig = 1 
  GOTO 11 
12 CONTINUE 
  igp = ig 
                                                                        
  RETURN 
END FUNCTION tr
                                                                        
FUNCTION tad_be92 (ak) 
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N) 
  save

  COMMON / be92_params / bet1, bet2, bet3, anu1, anu, anu_1, akfnmax                                                           
  COMMON / omg / omg, omx, omnr, omnu, Gamma_ebw 
  COMMON / h0 / h 
!      tad_be92=bet1*ak+(bet2*ak)**anu1+(bet3*ak)**2                    
  tad_be92 = bet1 * ak + (bet2 * ak) * anu1 + (bet3 * ak) **2 
  tad_be92 = tad_be92**anu 
  tad_be92 = 1.0 / (1.0 + tad_be92) **anu_1 
  RETURN 
END FUNCTION tad_be92
                                                                        
FUNCTION tfn_holtzmann (ak) 
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)
  save
  
  COMMON / holtzmann_paramsI / initialize
  COMMON / holtzmann_params / t2, t3, t4, t5, aktfnmax 
                                                                        
  COMMON / omg / omg, omx, omnr, omnu, Gamma_ebw 
  COMMON / h0 / h 
!      if(initialize.eq.1) then                                         
!         write(*,*) ' tfn parameters: t2,t3,t4,t5,aktfnmax '           
!         read(999,*) t2,t3,t4,t5,aktfnmax                                
!         initialize=0                                                  
!      endif                                                            
  sqak = sqrt (ak) 
  tfn_holtzmann = 1.0 / (1.0 + sqak * (t2 + sqak * (t3 + sqak * (t4 + sqak * t5) ) ) )                                            
  RETURN 
END FUNCTION tfn_holtzmann
                                                                        
                                                                        
FUNCTION tfn_holtzmann_steep (ak) 
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)
  save
  COMMON / holtzmann_paramsI / initialize
  COMMON / holtzmann_params / t2, t3, t4, t5, aktfnmax 
                                                                        
  COMMON / omg / omg, omx, omnr, omnu, Gamma_ebw 
  COMMON / h0 / h 
!      if(initialize.eq.1) then                                         
!         write(*,*) ' tfn parameters: t2,t3,t4,t5,aktfnmax '           
!         read(999,*) t2,t3,t4,t5,aktfnmax                                
!         initialize=0                                                  
!      endif                                                            
  tfn_holtzmann_steep = 1.0 / (1.0 + ak * (t2 + ak * (t3 + ak * (t4 + ak**3 * t5) ) ) )                                           
  RETURN 
END FUNCTION tfn_holtzmann_steep
                                                                        
                                                                        
!*******************************************************************    
FUNCTION filter (gm) 
!*******************************************************************    
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)                                                                        
  save
  PARAMETER (pi2 = 6.28318531, sq2 = 1.41421356) 
  PARAMETER (uscaleg_1 = 0.1, bertfil = 12.0) 
  REAL(dp) ns 
  COMMON / ifil / ifil 
  COMMON / quant / gmn, fl, ns 
  COMMON / h0 / h 
  u = gm * fl 
  filter = 0. 
!     fl IS IN REAL COMOVING DISTANCE UNITS                             
  IF (ifil.eq.1) then 
     IF (u.le.1.0) filter = 1. 
     RETURN 
  ELSEIF (ifil.eq.0) then 
     IF (fl.eq.0.) return 
     IF (u.le.5.0) then 
        filter = exp ( - u * u * 0.5) 
     ENDIF
     RETURN 
  ELSEIF (ifil.eq.2) then 
     IF (u.le.50.0) then 
        IF (u.lt.2.e-02) then 
           u2=u*u 
           filter=(1.-u2/10.*(1.-u2/28.))
!     this is 3*j_1/u
        ELSE
           u2=u*u 
           su=sin(u)
           cu=cos(u)
           filter=3.*(su-u*cu)/u2/u
        endif
     else
        filter=-3./sq2/u/u
     endif
     return
  else if(ifil.eq.4) then
     if(u.ge.2.0) return
     
     if(u.eq.0.) then 
        filter=3./pi2
        return
     endif
         
     u2=u*u
     qa=sqrt(4.-u2)
     qb=dlog((2.+qa)/u)
     filter=qa*(2.+3.25*u2)-u2*(6.+0.375*u2)*qb
     filter=filter/pi2
     if(u.ge.1.) return
     qa=sqrt(1.-u2)
     qb=dlog((1.+qa)/u)
     filter=filter-(qa*(1.+6.5*u2)-u2*(6.+1.5*u2)*qb)/pi2
     return
  endif
  if(ifil.eq.5) then
     if(u.le.25.0) then
        if(u.lt.2.e-02) then
           u2=u*u 
           filter=(1.-u2*0.1*(1.-u2/28.))
!     this is 3*j_1/u

        else
           u2=u*u 
           u_3=1.0/u2/u
           su=sin(u)
           cu=cos(u)
           filter=3.0*(su-u*cu)*u_3
        endif
        v=(u*uscaleg_1)**2
        amask=exp(-0.5*v)
        filter=filter*amask
     endif
     return
  endif
  if(ifil.eq.6) then
     if(u.le.25.0) then
        if(u.lt.2.e-02) then
           u2=u*u 
           filter=(1.-u2*0.1*(1.-u2/28.))
!     this is 3*j_1/u
        else
           u2=u*u 
           u_3=1.0/u2/u
           su=sin(u)
           cu=cos(u)
           filter=3.0*(su-u*cu)*u_3
        endif
        v=(gm*bertfil/h)**2
        amask=exp(-0.5*v)
        filter=filter*amask
     endif
     return
  endif
      
end FUNCTION filter
      
!*******************************************************************
function fil_dlnr(gm)
!*******************************************************************
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)
  save
!  real(dp) log
      
  parameter (pi2=6.28318531,sq2=1.41421356)
  parameter (uscaleg_1=0.1,bertfil=12.0)
!     dwin_dlnr = -udW/du
  real(dp) ns
  common/ifil/ifil
  common/quant/gmn,fl,ns
  common/h0/h
      
  u=gm*fl
  fil_dlnr=0.
!    fl IS IN REAL COMOVING DISTANCE UNITS
  if(ifil.eq.1) then
     return
  else if(ifil.eq.0) then
     if(fl.eq.0.) return
     if(u.le.5.0) then
        u2=u*u
        fil_dlnr=u2*exp(-u2*0.5)
     endif
     return
  else if(ifil.eq.2) then
     if(u.le.50.0) then
        if(u.lt.2.e-02) then
           u2=u*u 
           fil_dlnr=u2/5.*(1.-u2/14.)
!     this is 3*j_1/u
        else
           u2=u*u 
           su=sin(u)
           cu=cos(u)
           fil_dlnr=3.0*((3.0-u2)*su-3.0*u*cu)/u2/u
        endif
     else
        fil_dlnr=-6./sq2/u/u
     endif
     return
  else if(ifil.eq.4) then
     if(u.ge.2.0) return
     
     if(u.eq.0.) then 
        return
     endif
         
     u2=u*u
     qa=sqrt(4.-u2)
     qb=dlog((2.+qa)/u)
     fil_dlnr=-9.*qa*u2+u2*(12.+1.5*u2)*qb
     fil_dlnr=fil_dlnr/pi2
     if(u.ge.1.) return
     qa=sqrt(1.-u2)
     qb=dlog((1.+qa)/u)
     fil_dlnr=fil_dlnr+(18.*qa*u2-u2*(12.+6.*u2)*qb)/pi2
     return
  endif
  if(ifil.eq.5) then
     if(u.le.25.0) then
        if(u.lt.2.e-02) then
           u2=u*u 
           fil=(1.-u2*0.1*(1.-u2/28.))
           fil_dlnr=u2*0.2*(1.-u2/14.)
        else
           u2=u*u 
           u_3=1.0/u2/u
           su=sin(u)
           cu=cos(u)
           fil=3.0*(su-u*cu)*u_3
           fil_dlnr=3.0*((3.0-u2)*su-3.0*u*cu)*u_3
        endif
        v=(u*uscaleg_1)**2
        amask=exp(-0.5*v)
        fil_dlnr=(fil_dlnr+fil*v)*amask
     endif
     return
  endif
      
  if(ifil.eq.6) then
     if(u.le.25.0) then
        if(u.lt.2.e-02) then
           u2=u*u 
           fil=(1.-u2*0.1*(1.-u2/28.))
           fil_dlnr=u2*0.2*(1.-u2/14.)
        else
           u2=u*u 
           u_3=1.0/u2/u
           su=sin(u)
           cu=cos(u)
           fil=3.0*(su-u*cu)*u_3
           fil_dlnr=3.0*((3.0-u2)*su-3.0*u*cu)*u_3
        endif
        v=(gm*bertfil/h)**2
        amask=exp(-0.5*v)
        fil_dlnr=(fil_dlnr+fil*v)*amask
     endif
     return
  endif
      
end function fil_dlnr
      
!*******************************************************************
subroutine wind(gm,win)
!*******************************************************************
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)                                                                        
  save      
  common/h0/h
  common/lad/lad
  COMMON / finI / ifluc
  COMMON / fin / fin, finp 
  dimension win(100)
  win(1)=1.
  win(2)=1.
  win(3)=1.
  win(4)=h*h/gm/gm
  win(5)=gm**4/5.
  win(6)=gm*gm/(-3.)
  return
end subroutine wind
      
!*******************************************************************
subroutine windr(gm,r,win)
!*******************************************************************
  USE intreal_types
  IMPLICIT REAL(DP) (A-H,O-Z)
  IMPLICIT INTEGER(i4b) (I-N)                                                                        
  save
  parameter (fisoc=36.15)        
  common/ncorr/nc
  common/h0/h
  common/lad/lad
  COMMON / finI / ifluc
  COMMON / fin / fin, finp 
  dimension win(100)
  
  x=gm*r
  x2=x*x
  if(x.ge.0.049) then
     s=sin(x)
     c=cos(x)
  endif
  if(x.lt.0.05) then
     bj0=1.-x2/6.*(1.-x2/20.)
  else
     bj0=s/x
  endif
  if(x.lt.0.05) then
     bj1=x/3.*(1.-x2/10.*(1.-x2/28.))
  else
     bj1=(s-x*c)/x/x
  endif
  if(x.lt.1.0) then
     bj2=x2/15.*(1.-x2/14.*(1.-x2/36.))
  else
     bj2=(s*(3.-x2)-c*3.*x)/x/x2
  endif
  if(x.lt.1.08) then
     bj3=x2*x/105.*(1.-x2/18.*(1.-x2/44.))
  else
     bj3=5.*bj2/x-bj1
  endif
  if(x.lt.1.08) then
     bj4=x2*x2/945.*(1.-x2/22.*(1.-x2/52.))
  else
     bj4=7.*bj3/x-bj2
  endif
      
  win(1)=bj0
!     THIS GIVES xi
  if(x.eq.0.) then
     win(2)=1.
  else
     win(2)=3.*bj1/x
  endif
!     THIS GIVES 3*J_3/r^3
!     THIS ALSO GIVES (3/r) K_{-1,1}
  win(3)=win(2)**2
!     THIS GIVES DELTA(r)^2
  gm2=gm*gm
  gm4=gm2*gm2
  win(4)=bj0*gm2
!     THIS GIVES K_20
  win(5)=bj0*gm4
!     THIS GIVES K_40
  win(6)=bj1*gm4
!    THIS GIVES K_41
  win(7)=bj2*gm2
!    THIS GIVES K_22
  win(8)=bj2*gm4
!    THIS GIVES K_42
  win(9)=bj3*gm4 
!     THIS GIVES K_43
  win(10)=bj4*gm4 
!     THIS GIVES K_44
  win(11)=bj3/gm 
!    THIS GIVES K_{-1,3}
  win(12)=bj2/gm2 
!    THIS GIVES K_{-2,2}
  win(13)=bj0/gm2 
!    THIS GIVES K_{-2,0}
  return
end subroutine windr
     
!*************************************************************
     
      
