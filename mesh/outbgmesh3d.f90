!# Coded on 2016.09.22
subroutine outbgmesh3d(posfile,g_meshpara)
use horizontalresolution
use param_mesh
implicit none

type(meshpara),intent(in) :: g_meshpara
character(50),intent(in)  :: posfile
integer(4) :: ifile
real(8) :: si,sb
!###
integer(4) :: nx,ny,nz,i,j,k,jj,nobs
integer(4),parameter :: nxmax=300,nymax=300,nzmax=100
real(8) :: xgrd(nxmax+1),ygrd(nymax+1),zgrd(nzmax+1),dd,dx(-8:8)
real(8) :: v1,v2,v3,v4,v5,v6,v7,v8,x0,y0,z0,xx,yy,zz,r2,sigma,v(8)

!#[0]## set input
xgrd(1:4) = g_meshpara%xbound(1:4)
ygrd(1:4) = g_meshpara%ybound(1:4)
xgrd(1)   = xgrd(1) -20.0
xgrd(4)   = xgrd(4) +20.0
ygrd(1)   = ygrd(1) -20.0
ygrd(4)   = ygrd(4) +20.0
zgrd(4)   = g_meshpara%zmax + 20.0
zgrd(3)   = g_meshpara%upzin
zgrd(2)   = g_meshpara%downzin
zgrd(1)   = g_meshpara%zmin - 20.0
sigma = g_meshpara%sigma ! [km]
nobs = g_meshpara%nobs
r2=dsqrt(2.d0)

!####### Initial grid ######
ifile=1
nx=4
ny=4
nz=4
!#
dd= 0.1d0 ! 1m for threshold for distance between two grid lines
dx(-7:7)= (/-3.0,-2.0,-1.5,-1.0,-0.75,-0.5,-0.25,0.0,&
     &       0.25,0.50,0.75,1.0,1.5,2.0,3.0/) * sigma
jj=7

!#[1]## modify the grid
!goto 100
do i=1,nobs
 x0=g_meshpara%xyz(1,i) ! [km]
 y0=g_meshpara%xyz(2,i) ! [km]
 z0=g_meshpara%xyz(3,i) ! [km]
 write(*,*) "x0=",x0,"y0=",y0,"z0=",z0
 do j=-jj,jj
 !#[x]
  xx = x0 + dx(j)
  call updategrd(xgrd,nxmax,nx,xx,dd)
 !#[y]
  yy = y0 + dx(j)
  call updategrd(ygrd,nymax,ny,yy,dd)
 !#[z]
  zz = z0 + dx(j)
  call updategrd(zgrd,nzmax,nz,zz,dd)
 end do
!
end do
do i=1,nx
 write(*,*) "xgrd (",i,")=",xgrd(i)
end do
do i=1,ny
 write(*,*) "ygrd (",i,")=",ygrd(i)
end do
do i=1,nz
 write(*,*) "zgrd (",i,")=",zgrd(i)
end do


open(ifile,file=posfile)
write(ifile,*)'View "backgraound mesh"{'
do k=1,nz-1
 do j=1,ny-1
  do i=1,nx-1
   v1=value_3d(xgrd(i  ),ygrd(j+1),zgrd(k  ),g_meshpara)
   v2=value_3d(xgrd(i  ),ygrd(j  ),zgrd(k  ),g_meshpara)
   v3=value_3d(xgrd(i+1),ygrd(j  ),zgrd(k  ),g_meshpara)
   v4=value_3d(xgrd(i+1),ygrd(j+1),zgrd(k  ),g_meshpara)
  !
   v5=value_3d(xgrd(i  ),ygrd(j+1),zgrd(k+1),g_meshpara)
   v6=value_3d(xgrd(i  ),ygrd(j  ),zgrd(k+1),g_meshpara)
   v7=value_3d(xgrd(i+1),ygrd(j  ),zgrd(k+1),g_meshpara)
   v8=value_3d(xgrd(i+1),ygrd(j+1),zgrd(k+1),g_meshpara)
  v(1:8)=(/v1,v2,v3,v4,v5,v6,v7,v8/)
  call SHWRITE(ifile,xgrd(i),xgrd(i+1),ygrd(j),ygrd(j+1),zgrd(k),zgrd(k+1), v)
  end do
 end do
end do
write(ifile,*)"};"
close(ifile)
write(*,*) "### outbgmesh3d end! ###"
return
end subroutine outbgmesh3d
!######################################################### SHWRITE
subroutine SHWRITE(ifile,x1,x2,y1,y2,z1,z2,v)
! x is the horizontal axis, y is the vertical axis
implicit none
integer(4) :: ifile,i,j
real(8) :: x1,x2,y1,y2,z1,z2
real(8),dimension(8) :: x,y,z,v
x(1:8)=(/x1,x1,x2,x2,x1,x1,x2,x2/)
y(1:8)=(/y2,y1,y1,y2, y2,y1,y1,y2/)
z(1:8)=(/z1,z1,z1,z1,z2,z2,z2,z2/)
write(ifile,*) "SH(",x(1),",",y(1),",",z(1)
write(ifile,*)  (",",x(j),",",y(j),",",z(j),j=2,8)
write(ifile,*) "){",v(1)
write(ifile,*)  (",",v(i),i=2,8)
write(ifile,*)"};"
return
end subroutine SHWRITE

