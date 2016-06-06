module input_parameters

  character*128 &
       filein, fileout, pkfile, filterfile, etabfile

  integer &
       inputdens,ibatch,idovoid,ioutshear,&
       ifil,num_redshifts,&
       lkgaus,SFTon,netaf,FFTon,&
       iseedFFT,nkradv,nkangv,lsamlogv,iseedfv,&
       iwantdkFFT,iseed_dkfft,iwrap,&
       igetETA,&
       ivol,nlx,nly,nlz,&
       nbuff,next,ievol,igetD2F,&
       maskprev,igethom_strat,iwant_evpvtabhpkvd,Nfsc,ifclinv,&
       ichange_pspec,debug,wsmooth,ioutfield,NonGauss
  real &
       global_redshift,maximum_redshift,&
       Omx,OmB,Omhdm,fhdmclus,Omvac,Omcurv,h,gmfmin,gmfmax,&
       rfmin,rfmax,&
       gmLBFFT,gmUBFFT,&
       ensampv,alphv,gmLBv,gmUBv,&
       filter_RF,&
       bZfin,bZinit,&
       dcore_box,dL_box,cenx,ceny,cenz,&
       f_camp,&
       Fsmin,Fsmax,zinit_fac_ell,devtab,dpv_evtab,evtabmax,pv_evtabmax,&
       rmax2rs,fNL

  dimension nkradv(2), nkangv(2), lsamlogv(2), ensampv(2), alphv(2), &
       gmLBv(2), gmUBv(2), iseedfv(2)

  ! Other parameters

  real    dcurv, fbuffer_box
  real    biasold,biasnew

  integer nboxes, idocore, idoboxf, ifmt, iamcurved
  integer min_core_box, max_core_box
  logical verbose
  integer*8 ncxm,ncxp,ncym,ncyp,nczm,nczp

end module input_parameters
