!########################################################### outbgmesh14
subroutine outbgmesh14_obs(g_meshpara,posfile,pos)
use param_mesh
use horizontalresolution ! for value
use bgelem               ! 2018.08.28
implicit none
type(bgele),           intent(out)  :: pos       ! 2018.08.28
type(meshpara),        intent(in)   :: g_meshpara
character(50),         intent(in)   :: posfile
integer(4)                          :: ifile
real(8)                             :: si,sb
integer(4)                          :: nx,ny,i,j,ii,jj,nobs
integer(4),parameter                :: nxmax=300,nymax=300
real(8),   allocatable,dimension(:) :: xgrd,ygrd,dx
real(8)                             :: dd
real(8)                             :: v1,v2,v3,v4,x0,y0,z0,xx,yy,r2,sigma

!#[0]## set input
allocate( dx(-8:8)      ) ! 2018.08.28
allocate( xgrd(nxmax+1) ) ! 2018.08.28
allocate( ygrd(nymax+1) ) ! 2018.08.28
xgrd(1:4) = g_meshpara%xbound(1:4)
ygrd(1:4) = g_meshpara%ybound(1:4)
xgrd(1)   = xgrd(1) -20.0
xgrd(4)   = xgrd(4) +20.0
ygrd(1)   = ygrd(1) -20.0
ygrd(4)   = ygrd(4) +20.0
sigma     = g_meshpara%sigma ! [km]
nobs      = g_meshpara%nobs
r2        = dsqrt(2.d0)

!####### Initial grid ######
ifile=1
nx=4
ny=4
!#
dd= 0.1d0 ! 1m for threshold for distance between two grid lines
dx(-7:7)= (/-3.0,-2.0,-1.5,-1.0,-0.75,-0.5,-0.25,0.0,&
     &       0.25,0.50,0.75,1.0,1.5,2.0,3.0/) * sigma
jj=7

!#[1]## modify the grid
!goto 100
do i=1,nobs
 x0 = g_meshpara%xyz(1,i) ! [km]
 y0 = g_meshpara%xyz(2,i) ! [km]
 z0 = g_meshpara%xyz(3,i) ! [km]
write(*,*) "x0=",x0,"y0=",y0,"z0=",z0
!#[x]
do j=-jj,jj
 xx = x0 + dx(j)
 call updategrd(xgrd,nxmax,nx,xx,dd)
end do
!#[y]
do j=-jj,jj
 yy = y0 + dx(j)
 call updategrd(ygrd,nymax,ny,yy,dd)
end do
!
end do
do i=1,nx
 write(*,*) "xgrd (",i,")=",xgrd(i)
end do
do i=1,ny
 write(*,*) "ygrd (",i,")=",ygrd(i)
end do

call initbgele((nx-1)*(ny-1),pos) ! see bgele.f90 2018.08.28
open(ifile,file=posfile)
write(ifile,*)'View "backgraound mesh"{'
do j=1,ny-1
 do i=1,nx-1
  v1=value(xgrd(i  ),ygrd(j+1),g_meshpara)
  v2=value(xgrd(i  ),ygrd(j  ),g_meshpara)
  v3=value(xgrd(i+1),ygrd(j  ),g_meshpara)
  v4=value(xgrd(i+1),ygrd(j+1),g_meshpara)
  !
  call SQWRITE(ifile,xgrd(i),xgrd(i+1),ygrd(j),ygrd(j+1),(/v1,v2,v3,v4/))
  ii= (j-1)*(nx-1) + i ! 2018.08.28
  call setbgele(pos,ii,xgrd(i),xgrd(i+1),ygrd(j),ygrd(j+1),(/v1,v2,v3,v4/)) ! 2018.08.28
 end do
end do
write(ifile,*)"};"
close(ifile)

write(*,*) "### outbgmesh14_obs end! ###"
return
end

