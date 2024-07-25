program main
implicit none
real(8) :: x3(3),bxyz(3)
integer(4) :: it, itmax
real(8) :: dt=1.d0, t, tmax=70*60.d0
itmax=tmax/dt
write(*,*) "tmax=",tmax,"itmax=",itmax
x3(1:3)=(/ 0.d0,0.d0,-8.d0 /)
x3=x3*1.d3

open(1,file="ana_bxyz.dat")
do it=1,itmax
 t=it*dt
 call calbxyzana(x3,t,bxyz)
 write(*,*)"it=",it, "bxyz=",bxyz(1:3)
 write(1,*) t,bxyz(1:3)
end do
close(1)
end program

! uses the analytical solution of Ichihara et al. (2013)
subroutine calbxyzana(x3,t,bxyz)
implicit none
real(8),intent(in) :: x3(3) ! [m]
real(8),intent(in) :: t ! [sec]
real(8),intent(out) :: bxyz(3)
real(8) :: v3(3),Fz=-30000.d0 ! [nT]
real(8),parameter :: h=8000.d0, g=9.8d0
real(8),parameter :: sigma=3.3d0, pi=4.d0*atan(1.d0)
real(8),parameter :: sigma1=0.01d0 ! 100 Ohm medium beneath the layer
real(8),parameter :: mu=4.*pi*1.d-7
complex(8) :: iunit=(0.d0,1.d0)
complex(8) :: eta,a,a1,ah,ah2,bz,bh,sh,ch
real(8) :: k,c,c_h,omega

c=sqrt(g*h) ! [m/sec]
c_h=2./mu/sigma/h
k=2.*pi/(80.d0*1.d3) ! [rad/m]
omega=k*c ![rad/sec]

!# eta
eta=cdexp(iunit*(k*x3(1)-omega*t)) ! [m]

!# Tyler (2013)
bz=c/(c+iunit*c_h)*eta/h*Fz
bxyz(3)=real(bz)
bxyz(1)=real(-iunit*bz)
bxyz(2)=0.d0

!# Ichihara et al. (2013)
      a =cdsqrt(k**2. - iunit*omega*mu*sigma)
      a1=cdsqrt(k**2. - iunit*omega*mu*sigma1)
      ah=a*h
      ah2=a*h/2.d0
      bz=-iunit*omega*mu*sigma/a*              &
        2.*sh(ah2)*(a*ch(ah2)+k*sh(ah2))/(       &
        (a1*k+a**2.)*sh(ah)+(a1*a+a*k)*ch(ah) &
        )*Fz*eta/h      
      bh=iunit*(a1/k)*bz
	bxyz(1)=real(bh)
	bxyz(3)=real(bz)
return
end

!########################################################### sinh
      function sh(a)
      implicit real(selected_real_kind(8))(a-h,o-z)
      complex(8) :: sh,a
      sh=(cdexp(a)-cdexp(-a))/2.d0
!      write(*,*) "a=",a
!      write(*,*) "cdexp(a)=",cdexp(a)
!      write(*,*) "sh=",sh
      return
      end
!########################################################### cosh
      function ch(a)
      implicit real(selected_real_kind(8))(a-h,o-z)
      complex(8) :: ch,a
      ch=(cdexp(a)+cdexp(-a))/2.d0
      return
      end 
