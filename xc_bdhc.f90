!         USE xc_pz_mod
!         IMPLICIT NONE
!         DOUBLE PRECISION, PARAMETER :: rs = 4.0, xi = 1
!         DOUBLE PRECISION :: exc0, dexc0
!
!         exc0 = exc0_pz(rs,xi)
!         dexc0 = dexc0_pz(rs,xi)
!         write(*,*) exc0, dexc0


      MODULE constants
          ! Global Constants
          integer, parameter :: MAX_INTEGER_DIGITS = 10
      END MODULE constants

      PROGRAM main

         USE xc_bdhc_mod
         USE xc0_pw_mod
         USE constants           ! Constants defined above

         ! Disable implicit declarations (i-n rule)
         IMPLICIT NONE
         DOUBLE PRECISION :: rs, t, xi
         DOUBLE PRECISION :: exc0, dexc0, exc, dexc, drs

         ! Variable defintions
         CHARACTER(MAX_INTEGER_DIGITS) :: rs_str, t_str, xi_str

         ! First command line argument is the base, second is the exponent
         call getarg(1, rs_str)
         call getarg(2, t_str)
         call getarg(3, xi_str)

         ! Make sure user provided both base and exponent
         if (rs_str=='' .OR. t_str=='' .OR. xi_str=='') then
             stop 'Usage: rs t xi'
         endif

         ! Convert strings to integers
         read (rs_str, *) rs
         read (t_str, *) t
         read (xi_str, *) xi

         ! Compute
         exc0 = exc0_pw(rs,xi)
         drs = 1.e-6
         dexc0 = (exc0_pw(rs+drs,xi)-exc0)/drs
         exc = exc_bdhc(rs,t,xi,exc0)
         dexc = dexc_bdhc(rs,t,xi,exc0,dexc0)
         write(*,*) exc, dexc

      END PROGRAM main
