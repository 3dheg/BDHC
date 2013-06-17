!         USE xc_pz_mod
!         IMPLICIT NONE
!         DOUBLE PRECISION, PARAMETER :: rs = 4.0, xi = 1
!         DOUBLE PRECISION :: exc0, dexc0
!
!         exc0 = exc0_pz(rs,xi)
!         dexc0 = dexc0_pz(rs,xi)
!         write(*,*) exc0, dexc0


      PROGRAM main

         USE xc_bdhc_mod
         USE xc0_pw_mod
         IMPLICIT NONE
         DOUBLE PRECISION, Parameter :: rs = 4.0, th = 1.0, xi = 1
         DOUBLE PRECISION :: exc0, dexc0, exc, dexc, drs

         exc0 = exc0_pw(rs,xi)
         drs = 0.0001
         dexc0 = (exc0_pw(rs+drs,xi)-exc0)/drs
         exc = exc_bdhc(rs,th,xi,exc0)
         dexc = dexc_bdhc(rs,th,xi,exc0,dexc0)
         write(*,*) 'Exc0',exc0
         write(*,*) 'Exc',exc
         write(*,*) 'dExc0/drs',dexc0
         write(*,*) 'dExc/drs',dexc

      END PROGRAM main
