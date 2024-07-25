!# Coded on September 20, 2016

module solution_ana
use param ! include "precision.inc"
use constants
implicit none

contains
!###########################################################
subroutine outana(g_param)
implicit none
type(param_forward),intent(in) :: g_param
real(8),parameter :: dt = 1.0
real(8)       :: x3obs(3,g_param%nobs),x3(3), bxyz(3),sigma_ocean, sigma_crust
real(8)       :: depth,lambda,eta,Fz,tstart,tmax,t
integer(4)    :: itmax, it, ncond, i, nobs, iflag_geomag
character(70) :: outf, outfile, obsname(g_param%nobs),name

!#[0]## set input
iflag_geomag = g_param%iflag_geomag
depth        = g_param%h_anapara%depth * 1.d3 ! [km] -> [m] ocean depth
lambda       = g_param%h_anapara%lambda
eta          = g_param%h_anapara%eta
Fz           = g_param%geomagvector(3)   ! [nT] (upward positive)
tstart       = g_param%tstart
tmax         = g_param%tmax
outf         = g_param%outbxyzfolder ! 2018.11.14
nobs         = g_param%nobs
obsname      = g_param%obsname
x3obs        = g_param%xyzobs     ! [km]
sigma_ocean  = g_param%cond(2) ! [S/m]
sigma_crust  = g_param%cond(3) ! [S/m]
ncond        = g_param%ncond
itmax        = tmax/dt

if ( iflag_geomag .ne. 2 ) goto 100
if ( ncond .ne.        3 ) goto 101

!#[1]##
write(*,*) "tmax=",tmax,"itmax=",itmax

!#[2]## observatory loop
do i=1,nobs

 !#[2-1]##
 name = obsname(i)
 outfile=outf(1:len_trim(outf))//name(1:len_trim(name))//"_bxyz_ts_ana.dat"
 x3(1:3)=x3obs(1:3,i)*1.d3 ! [km] -> [m]

 open(1,file=outfile)
 do it=1,itmax
  t=it*dt
  call calbxyzana(x3,t,Fz,depth,lambda,eta,sigma_ocean,sigma_crust,bxyz)
!  write(*,*)"it=",it, "bxyz=",bxyz(1:3)
  write(1,'(4g15.7)') t,bxyz(1:3)
 end do
 close(1)
end do

return

write(*,*) "### OUTANA END!! ###"

100 continue
write(*,*) "GEGEGE! For analycical solution, static geomag field should be provided."
write(*,*) "iflag_geomag=",iflag_geomag,"should be 2!"
stop

101 continue
write(*,*) "GEGEGE! ncond",ncond,"is not supported in analytical solution in calbxyzana!"
write(*,*) "ncond should be 2!"
stop

end
!########################################################### calbxyzana
! Modified on 2016.09.20
! uses the analytical solution of Ichihara et al. (2013)
subroutine calbxyzana(x3,t,Fz,h,lambda,amp,sigma,sigma1,bxyz)
implicit none
real(8),intent(in) :: x3(3)  ! [m]
real(8),intent(in) :: t      ! [sec] time
real(8),intent(in) :: amp    ! [m] tsunami height amplitude
real(8),intent(in) :: Fz     ! [nT]
real(8),intent(in) :: h      ! [m] ocean depth
real(8),intent(in) :: lambda ! [km] wave length
real(8),intent(in) :: sigma, sigma1 ! ocean and crust conductivity [S/m], respectively
real(8),intent(out) :: bxyz(3)
real(8),   parameter :: g=9.8d0
complex(8),parameter :: iunit=(0.d0,1.d0)
complex(8) :: eta,a,a1,ah,ah2,bz,bh
real(8)    :: k,c,c_h,omega
!complex(8),external :: sh,ch

c=sqrt(g*h) ! [m/sec]
c_h=2./dmu/sigma/h
k=2.*pi/(lambda*1.d3) ! [rad/m]
omega=k*c ![rad/sec]

!# eta
eta=amp*cdexp(iunit*(k*x3(1)-omega*t)) ! [m]

!# Tyler (2013)
bz=c/(c+iunit*c_h)*eta/h*Fz
bxyz(3)=real(bz)
bxyz(1)=real(-iunit*bz)
bxyz(2)=0.d0

!# Ichihara et al. (2013)
      a =cdsqrt(k**2. - iunit*omega*dmu*sigma)
      a1=cdsqrt(k**2. - iunit*omega*dmu*sigma1)
      ah=a*h
      ah2=a*h/2.d0
      bz=-iunit*omega*dmu*sigma/a*              &
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

end module solution_ana
