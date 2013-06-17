
      MODULE xc0_pw_mod
!   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!   Perdew-Zunger
!
!      citation
!
!----------------------------------------------------------------------
      IMPLICIT NONE
      CONTAINS

!    Paolo Gio Giorgi:
!
!       citations
!
! EXAMPLE: use the subroutine spinres to compute
! the spin-resolved correlation energies at z=0.4
! as a function of rs, and write them on a file
!
!      program test
!      implicit none
!      double precision rs,z,ecud,ecuu,ecdd
!      integer i
!      open(8,file='ecss_z0.4',status='unknown')
!      z=0.4d0
!      do i=1,100
!         rs=0.5d0*i
!         call spinres(rs,z,ecud,ecuu,ecdd)
!         write(8,*) rs,ecud,ecuu,ecdd
!      enddo
!      stop
!      end
!
! END EXAMPLE
!!!!!!!!!!!!!!
      subroutine spinres(rs,z,ecud,ecuu,ecdd)
! spin-resolved correlation energy of 3D electron gas of 
! density parameter rs and spin polarization z
! Hartree atomic units used
! ec = ecud + ecuu + ecdd
! ud = up-down+down-up
! uu = up-up
! dd = down-down
! Gori-Giorgi & Perdew, PRB 69, 041103 (2004)
      implicit none
      double precision rs,z,ecud,ecuu,ecdd
      double precision Fuu,Fud,Fdd,ec
      if(abs(z).gt.1.d0) then 
         write(*,*) 'bad z input'
      else if(z.eq.1.d0) then
         Fuu=1.d0
         Fdd=0.d0
         Fud=0.d0
      else if(z.eq.-1.d0) then
         Fuu=0.d0
         Fdd=1.d0
         Fud=0.d0
      else
         call frac(rs,z,Fuu)
         call frac(rs,-z,Fdd)
         Fud=1.d0-Fuu-Fdd
      endif
      call ecPW(rs,z,ec)
      ecud=Fud*ec
      ecuu=Fuu*ec
      ecdd=Fdd*ec
      return
      end

      subroutine frac(rs,z,Fuu)
      implicit none
      double precision rs,z,Fuu
      double precision cPW0,ff,IPW,FuuHD,FuuLD,FuuSKTP,dLD,den,cf,pi,  &
                       ec1,ecz,Auu,B,C
      pi=acos(-1.d0)
      cf=(9.d0*pi/4.d0)**(1.d0/3.d0)
      B = 0.178488d0
      C = 2.856d0   
! high-density limit
      ff=((1.d0+z)**(4.d0/3.d0)+(1.-z)**(4.d0/3.d0)-2.d0)/(2.d0        &
         **(4.d0/3.d0)-2.d0)
      cPW0=0.03109d0-0.00988d0*ff*(1.d0-z**4)-0.01555d0*ff*z**4
      IPW=cPW0/0.03109d0
      FuuHD=(1.d0+z)/4.d0/IPW
! low-density limit
      dLD=0.8986d0
      den=3.d0/8.d0/pi*cf*((1.d0+z)**(4.d0/3.d0)+(1.d0-z)**(4.d0/3.d0))&
          -dLD
      FuuLD=(3.d0/8.d0/pi*cf*(1.d0+z)**(4.d0/3.d0)-dLD*((1.d0+z)/2.d0) &
           **2)/den
! SKTP
      call ecPW(3.28d0,1.d0,ec1)
      call ecPW(3.28d0,z,ecz)
      FuuSKTP=((1.d0+z)/2.d0)**(11.d0/6.d0)*ec1/ecz
! interpolation formulas
      Auu=(FuuSKTP-FuuHD)/sqrt(3.28d0)+C*FuuSKTP+B*sqrt(3.28d0)*       &
          (FuuSKTP-FuuLD)
      Fuu=(FuuHD+Auu*sqrt(rs)+FuuLD*rs*B)/(1.d0+C*sqrt(rs)+B*rs)
      return
      end


      subroutine ecPW(x,y,ec)
! correlation energy of 3D electron gas of density parameter rs
! and spin polarization z
! in Hartree atomic units: ec=ec(rs,zeta)
! x -> rs; y -> zeta
! Perdew & Wang, PRB 45, 13244 (1992)
      implicit none
      double precision x,y,ec
      double precision pi,f02,ff,aaa,G,ec0,ec1,alfac
      pi=dacos(-1.d0)
      
      f02=4.d0/(9.d0*(2.d0**(1.d0/3.d0)-1.d0))

      ff=((1.d0+y)**(4.d0/3.d0)+(1.d0-y)**(4.d0/3.d0)-                 &
         2.d0)/(2.d0**(4.d0/3.d0)-2.d0)

      aaa=(1.d0-dlog(2.d0))/pi**2
      call  GPW(x,aaa,0.21370d0,7.5957d0,3.5876d0,                     &
                1.6382d0,0.49294d0,G)
      ec0=G

      aaa=aaa/2.d0
      call GPW(x,aaa,0.20548d0,14.1189d0,6.1977d0,                     &
               3.3662d0,0.62517d0,G)
      ec1=G
      call GPW(x,0.016887d0,0.11125d0,10.357d0,3.6231d0,               &
               0.88026d0,0.49671d0,G)
      alfac=-G

      ec=ec0+alfac*ff/f02*(1.d0-y**4)+(ec1-ec0)*ff*y**4
!f2py intent(in, out) ec
      return
      end

      subroutine GPW(x,Ac,alfa1,beta1,beta2,beta3,beta4,G)
      implicit none
      double precision G,Ac,alfa1,beta1,beta2,beta3,beta4,x
      G=-2.d0*Ac*(1.d0+alfa1*x)*dlog(1.d0+1.d0/(2.d0*                  &
         Ac*(beta1*x**0.5d0+                                           &
         beta2*x+beta3*x**1.5d0+beta4*x**2)))
      return
      end

!----------------------------------------------------------------------
!
        SUBROUTINE exHF(rs,xi,ex0)
          IMPLICIT NONE
          DOUBLE PRECISION :: rs, xi, ex0

          IF (xi.eq.1.d0) THEN
            ex0 = -1.1545/rs
          ELSE
            ex0 = -0.916331/rs
          ENDIF
          CONTINUE
        END


        REAL FUNCTION exc0_pw(rs,xi)
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(IN) :: rs, xi
          DOUBLE PRECISION :: ex0, ec0

          call exHF(rs,xi,ex0)
          call ecPW(rs,xi,ec0)

          exc0_pw = ex0 + 2.d0*ec0
        END FUNCTION exc0_pw

      END MODULE xc0_pw_mod
