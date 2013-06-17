      MODULE xc_bdhc_mod
!----------------------------------------------------------------------
!
!   DOUBLE PRECISION FUNCTION exc_bdhc(rs,th,xi,exc0,dexc0)
!
!   INPUT: rs = Wigner-Seitz radius, scaled by the Bohr radius
!          th = Ratio of temperature to Fermi temperature
!          xi = Spin-polarization (1 = polarized, 0 = unpolarized)
!          exc0 = Exchange-correlation energy for the same system at 0T
!
!   RETURNS: exc_bdhc = Exchange-correlation energy (Rydbergs)
!
!
!   DOUBLE PRECISION FUNCTION dexc_bdhc(rs,th,xi,exc0,dexc0)
!
!   INPUT: rs = Wigner-Seitz radius, scaled by the Bohr radius
!          th = Ratio of temperature to Fermi temperature
!          xi = Spin-polarization (1 = polarized, 0 = unpolarized)
!          exc0 = Exchange-correlation energy for the same system at 0T
!          dexc0 = rs derivative of exc0
!
!   RETURNS: dexc_bdhc = rs derivative of exc_bdhc (Rydbergs/rs)
!
!   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!   EXAMPLE USAGE:
!
!      PROGRAM main
!
!         USE xc_bdhc_mod
!         IMPLICIT NONE
!         DOUBLE PRECISION, PARAMETER :: rs = 4.0, th = 1.0, xi = 1
!         DOUBLE PRECISION, PARAMETER :: exc0 = -0.323254
!         DOUBLE PRECISION, PARAMETER :: dexc0 = 0.07655
!         DOUBLE PRECISION :: exc, dexc
!
!         exc = exc_bdhc(rs,t,xi,exc0)
!         dexc = dexc_bdhc(rs,t,xi,exc0,dexc0)
!         write(*,*) exc, dexc
!
!      END PROGRAM main
!
!   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!   If you use this module in your calculations, please cite both the
!   fit and original Monte Carlo simulation:
!
!      E. W. Brown, J. L. DuBois, M. Holzmann and D. M. Ceperley
!      Exchange-correlation energy for the 3D homogeneous electron gas
!      at arbitrary temperature
!      Submitted to Phys. Rev. B. Rapid Communications (2013)
!
!      E. W. Brown, B. K. Clark, J. L. DuBois, and D. M. Ceperley
!      Path Integral Monte Carlo simulation of the warm-dense
!      homogeneous electron gas
!      Phys. Rev. Lett. 110, 146405 (2013)
!
!----------------------------------------------------------------------
      IMPLICIT NONE
      CONTAINS

        SUBROUTINE params(rs,xi,a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3)
          DOUBLE PRECISION :: rs, xi
          DOUBLE PRECISION :: a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3
          IF (xi == 1) THEN
             IF (rs > 10.0) THEN
                a1 = -7.23836
                a2 = -6.65715
                a3 = -5.89226
                b1 = 19.8258
                b2 = 19.9802
                b3 = 17.3632
                c1 = 0.254584
                c2 = 0.263629
                c3 = 0.238536
                d1 = 0.0521708
                d2 = 0.0540244
                d3 = 0.0488823
             ELSE
                a1 = -1.57839
                a2 = -1.46754
                a3 = -0.784554
                b1 = -9.99823
                b2 = -11.3387
                b3 = -11.5341
                c1 = 7.10336
                c2 = 7.85547
                c3 = 7.07407
                d1 = -2.19297
                d2 = -2.40187
                d3 = -2.17553
             END IF
          ELSE
             IF (rs > 10.0) THEN
                a1 = 4.38637
                a2 = 5.96304
                a3 = 5.43786
                b1 = 1.22928
                b2 = 0.249599
                b3 = -1.10198
                c1 = -0.789404
                c2 = -0.991637
                c3 = -0.716191
                d1 = 0.178368
                d2 = 0.220769
                d3 = 0.157061
             ELSE
                a1 = 3.94068
                a2 = 5.59666
                a3 = 8.19611
                b1 = -0.330048
                b2 = -1.39311
                b3 = -2.43483
                c1 = -0.0381205
                c2 = -0.254872
                c3 = -1.7384
                d1 = -0.0356196
                d2 = 0.00877504
                d3 = 0.383061
             END IF
          END IF
        END SUBROUTINE params

        DOUBLE PRECISION FUNCTION exc_bdhc(rs,th,xi,exc0)
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(IN) :: rs,th,xi,exc0
          DOUBLE PRECISION, PARAMETER :: PI = 3.14159265359
          DOUBLE PRECISION :: a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3
          DOUBLE PRECISION :: u1,u2,M1,M2,M3,P1,P2,TF,T

          call params(rs,xi,a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3)

          u1 = 3./(rs**3.)
          u2 = sqrt(6./(rs**3.))
          M1 = exp(a1*log(rs) + b1 + c1*rs + d1*rs*log(rs))
          M2 = exp(a2*log(rs) + b2 + c2*rs + d2*rs*log(rs))
          M3 = exp(a3*log(rs) + b3 + c3*rs + d3*rs*log(rs))
          TF = (9.*PI/(2.*(1.-xi)+2.))**(2./3.) / (rs**2.)
          T = th * TF
          P1 = (M2*u1 + M3*u2)*(T**2.) + M2*u2*(T**(5./2.))
          P2 = 1. + M1*(T**2.) + M3*(T**(5./2.)) + M2*(T**3.)

          exc_bdhc = (exc0 - P1)/P2
        END FUNCTION exc_bdhc

        DOUBLE PRECISION FUNCTION dexc_bdhc(rs,th,xi,exc0,dexc0)
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(IN) :: rs,th,xi,exc0,dexc0
          DOUBLE PRECISION, PARAMETER :: PI = 3.14159265359
          DOUBLE PRECISION :: a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3
          DOUBLE PRECISION :: u1,u2,M1,M2,M3,P1,P2,TF,T
          DOUBLE PRECISION :: du1,du2,dM1,dM2,dM3,dP1,dP2

          call params(rs,xi,a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3)

          u1 = 3./(rs**3.)
          u2 = sqrt(6./(rs**3.))
          M1 = exp(a1*log(rs) + b1 + c1*rs + d1*rs*log(rs))
          M2 = exp(a2*log(rs) + b2 + c2*rs + d2*rs*log(rs))
          M3 = exp(a3*log(rs) + b3 + c3*rs + d3*rs*log(rs))
          TF = (9.*PI/(2.*(1.-xi)+2.))**(2./3.) / (rs**2.)
          T = th * TF
          P1 = (M2*u1 + M3*u2)*(T**2.) + M2*u2*(T**(5./2.))
          P2 = 1. + M1*(T**2.) + M3*(T**(5./2.)) + M2*(T**3.)

          du1 = -9./(rs**4.)
          du2 = -3.*sqrt(3./(2.*(rs**5.)))
          dM1 = M1*((a1/rs) + c1 + d1*(1.+log(rs)))
          dM2 = M2*((a2/rs) + c2 + d2*(1.+log(rs)))
          dM3 = M3*((a3/rs) + c3 + d3*(1.+log(rs)))
          dP1 = (dM2*u1 + M2*du1 + dM3*u2 + M3*du2)*(T**2.)          &
                + (dM2*u2 + M2*du2)*(T**(5./2.))
          dP2 = dM1*(T**2.) + dM3*(T**(5./2.)) + dM2*(T**3.)

          dexc_bdhc = ((dexc0 - dP1)/P2) - ((exc0 - P1)/P2)*(dP2/P2)
        END FUNCTION dexc_bdhc

      END MODULE xc_bdhc_mod

