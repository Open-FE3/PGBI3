!# This subroutine is used along with
!# mk3dgeo.f90 to generate gridded .pos file for background data
!# under coding on May 30, 2016
!
!########################################### OUTPOSFILE2_OBS
! This subroutine is added on May 27, 2016
! copied from volcano/3D_ana_comp/mesh_for_FEM/3dgeo.f90
!
subroutine OUTPOSFILE(bgmeshfile)
use param
implicit none
include "meshpara.h"
character(50),intent(in) :: bgmeshfile
type(param_forward),intent(in) :: g_param
real(8),intent(in) :: height1, height2, height3, height4
real(8),intent(in) :: lc1,lc2,lc3, dlen,dlen2,xout,yout,shift
character(50) :: outpos
integer(4) ::  ifile=2
real(8) :: x0,y0,z0,x,y,z,dd,dx(-3:3)
real(8),dimension(4) :: s
!integer(4) :: nx,ny,nz
integer(4) :: nx,ny,nz
integer(4),parameter :: nxmax=50,nymax=400,nzmax=50
real(8) :: xgrd(nxmax+1),ygrd(nymax+1),zgrd(nzmax+1)
!real(8),dimension(:),allocatable :: xgrd,ygrd,zgrd
integer(4) :: i,j,k,kk,jj
real(8) :: s1,s2,s3,s4,s5,s6,s7,s8,value_OBS,v(8)
!##
y2=earthrad*d2r*(nlat-latorigin)
y1=-y2
x2=earthrad*dcos(latorigin*d2r)*d2r*(elon-lonorigin)
x1=-x2
x01=x1-lenout-1
x02=x2+lenout+1
y01=y1-lenout-1
y02=y2+lenout+1
!######
nx=12
ny=12
nz=27
!####### Initial grid ######
xgrd(1:12)=(/x01,x1,4/5*x1,3/5*x1,2/5*x1,1/5*x1,
	    &  1/5*x2,2/5*x2,3/5*x2,4/5*x2,x2,x02/)
	    !#
ygrd(1:12)=(/y01,y1,4/5*y1,3/5*y1,2/5*y1,1/5*y1,
	    &  1/5*y1,2/5*y1,3/5*y1,4/5*x1,y2,y02/)
	    !#
zgrd(1:27)=(/hdown,-30.d0,-20.d0,-15.d0,-5.d0,-1.d0,-0.5d0,-0.25d0,-0.2d0,-0.15d0,-0.1d0,&
      &       -0.01d0,-0.003d0,-0.002d0,-0.001d0, 0.d0,0.001d0,0.002d0,0.003d0,0.01d0,&
	&         0.1d0,0.5d0,1.d0,5.d0,15.d0,20.d0,30.d0,hup/)!down to up

!#########
do i=1,nx
write(*,*) "xgrd (",i,")=",xgrd(i)
end do
do i=1,ny
write(*,*) "ygrd (",i,")=",ygrd(i)
end do
do i=1,nz
write(*,*) "zgrd (",i,")=",zgrd(i)
end do

!#[2]## calculate values and output .pos file
open(ifile,file=outpos)
write(ifile,*) 'View "backgraound mesh"{'
do k=1,nz-1
 do j=1,ny-1
  do i=1,nx-1                                                 !     z   y       x
  s1=value(xgrd(i)  ,ygrd(j+1),zgrd(k+1)) ! 1: up   top     left
  s2=value(xgrd(i)  ,ygrd(j)  ,zgrd(k+1)) ! 2: up   bottom  left
  s3=value(xgrd(i+1),ygrd(j)  ,zgrd(k+1)) ! 3: up   bottom  right
  s4=value(xgrd(i+1),ygrd(j+1),zgrd(k+1)) ! 4: up   top     right
  s5=value(xgrd(i)  ,ygrd(j+1),zgrd(k)  ) ! 5: down top     left
  s6=value(xgrd(i)  ,ygrd(j)  ,zgrd(k)  ) ! 6: down bottom  left
  s7=value(xgrd(i+1),ygrd(j)  ,zgrd(k)  ) ! 7: down bottom  right
  s8=value(xgrd(i+1),ygrd(j+1),zgrd(k)  ) ! 8: down top     right
  v(1:8)=(/s1,s2,s3,s4,s5,s6,s7,s8/)
  call SHWRITE(ifile,xgrd(i),xgrd(i+1),ygrd(j),ygrd(j+1),zgrd(k+1), zgrd(k), v)
  end do
 end do
end do
write(2,*) "};"
close(2)
write(*,*) "### OUTPOSFILE END!! ###"

return
end
!#######################################################    value
function value(x,y,z)
implicit none
include "meshpara.h"
real(8),intent(in) :: x,y,z,shift ! [km]
real(8),intent(in) :: lc1,lc3,dlen,dlen2
real(8) :: x1,x2,y1,y2,xyz(3,8),h0,hup,hdown
real(8) :: x01,x02,y01,y02
real(8) :: usreso,ureso,dsreso,dreso,breso
real(8) :: value_OBS, value_OBS_def
real(8) :: r,k,A,B,robs,robsmin,power
integer(4) :: i
r=dsqrt(x**2. + y**2. + z**2.)
breso=500.d0
h0=0.d0
!##
hup=upheight
hdown=downheight
ureso=upresolution
usreso=upsurfaceresolution
dsreso=downsurfaceresolution
dreso=downresolution
y2=earthrad*d2r*(nlat-latorigin)
y1=-y2
x2=earthrad*dcos(latorigin*d2r)*d2r*(elon-lonorigin)
x1=-x2
x01=x1-lenout-1
x02=x2+lenout+1
y01=y1-lenout-1
y02=y2+lenout+1

reso_surface=13.d0
!# [1] resolution = A*(z**k) + B, where k should be given as parameter
! A= (zmax -z_surf)/(breso - sreso)
! B= sreso - A*(z_sur**k)
! value=A*(z**k) + B
  A= (zmax -z_surf)/(breso - sreso)
  B= reso_surf - A*(z_sur**k)

!# horizontal position
if (  x1 .lt. x .and. x .lt. x2 .and. y1 .lt. y .and. y .lt. y2) then

end if
!# [2] A*r^3+B
!  shift=10.d0
!  power=3.d0  ! power should be around 2.5 ~ 2.7 when solved in local PC on May 15, 2016
  power=2.5d0
  A=(lc3-lc1)/((dlen2+shift)**power - (dlen+shift)**power)
  B=lc1-A*((dlen+shift)**power)
  value_OBS=A*((r+shift)**power)+B
!  if ( dabs(x) .lt. 1. .and. dabs(y) .lt. 1. ) then
!   if (dabs(z-(-0.2)) .lt. 0.051 ) value_OBS=0.5*value_OBS
!   if (dabs(z-(-0.2)) .lt. 0.02 ) value_OBS=0.02d0
!  end if
end if
!## for each observatory
value_OBS_def=value_OBS
return
end
