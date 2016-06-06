      SUBROUTINE collapse(z, omegam, omegal, w, delc, dcoll)
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     %     collapse.f                                      %
c     %                                                     %
c     %     This subroutines calculates:                    % 
c     %     (1) linear density contrast threshold for a     %
c     %        collapsed halo at z, delc(z), and            % 
c     %     (2) spherical overdensity of a collapsed halo   %
c     %        at z, dcoll(z).                              %
c     %         Note that dcoll(z) is the overdensity with  %
c     %        respect to the "mean mass density" of the    %
c     %        universe.                                    %
c     %                                                     % 
c     %     References:                                     %
c     %        Lacey & Cole (1993); Nakamura & Suto (1997)  %
c     %                                                     %
c     %                           1999. 6.28     E.Komatsu  %
c     %		(subroutine)	  2000.11.19     E.Komatsu  %
c     %         (quintessense)    2002. 1.23     E.Komatsu  %
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      implicit none
c
c output
c
      real*8 delc  ! linear density contrast threshold at z
      real*8 dcoll ! spherical overdensity of a collapsed halo at z
c
c inputs
c 
      real*8 z                 ! redshift
      real*8 omegam, omegal, w ! cosmological parameters 
c
      real*8 x, x2, x3, x3w
      real*8 omega             ! Omega_m(z)
      real*8 eta
      real*8 pi
      parameter (pi= 3.14159265358979d0)
c
      x  = 1d0+z
      x2 = x**2d0
      x3 = x**3d0
      x3w= x3**w
c
      omega = omegam*x3/( omegam*x3 + (1d0-omegam-omegal)*x2 + omegal )
      eta= dlog( (2d0/omega-1d0) + dsqrt((2d0/omega-1d0)**2d0-1d0) )
c
      if ( (omegam.lt.1d0).and.(omegal.eq.0d0) ) then
         delc= 1.5d0
     &        *(3d0*dsinh(eta)*(dsinh(eta)-eta)
     &         /(dcosh(eta)-1d0)**2d0-2d0)
     &        *(1d0+(2d0*pi/(dsinh(eta)-eta))**(2d0/3d0))
         dcoll= 4d0*pi*pi
     &         *(dcosh(eta)-1d0)**3d0/(dsinh(eta)-eta)**2d0
      else
         delc= 3d0*(12d0*pi)**(2d0/3d0)/20d0
     &          *(1d0+0.0123d0*dlog10(omega))
         dcoll= 18d0*pi*pi
     &         *(1d0+0.4093d0*(1d0/omega-1d0)**0.9052d0)
      endif

      return
      END
