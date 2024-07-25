! Under coding on May 19, 2016
! for use in mk3dgeo.f90

!########################################
subroutine genbgmesh9_grd(bgmeshfile)
implicit none
include "meshpara.h"
real(8) :: x1,x2,y1,y2,xyz(3,8),h0,hup,hdown
real(8) :: x01,x02,y01,y02
real(8) :: usreso,ureso,dsreso,dreso,breso
integer(4) :: i,j
character(50) :: bgmeshfile
integer(4),parameter :: nx=4,ny=4,nz=20
real(8) :: xgrd(nx+1),ygrd(ny+1),zgrd(nz+1)
breso=500.d0
h0=0.d0
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
xyz(1:2,1)=(/x1,y1/)
xyz(1:2,2)=(/x2,y1/)
xyz(1:2,3)=(/x2,y2/)
xyz(1:2,4)=(/x1,y2/)
xyz(1:2,5:8)=xyz(1:2,1:4)
xyz(3,5:8)=hup ! top
xyz(3,1:4)=h0 !z=0 surface
i=1
!# uplarge hexahedron
x01=x1-lenout-1
x02=x2+lenout+1
y01=y1-lenout-1
y02=y2+lenout+1
!#################################################################
xgrd(1:nx+1)=(/x01,x1,0.d0,x2,x02/)
ygrd(1:ny+1)=(/y01,y1,0.d0,y2,y02/)
zgrd(1:nz+1)=(/zmax,upheight,0.d0,downheight,zmin/)!down to up
!#################################################################
open(ifile,file=outpos)
write(ifile,*) 'View "backgraound mesh"{'
do k=1,nz
 do j=1,ny
  do i=1,nx                                                  !     z   y       x
  s1=value(xgrd(i)  ,ygrd(j+1),zgrd(k+1),lc1,lc3,dlen,dlen2,shift) ! 1: up   top     left
  s2=value(xgrd(i)  ,ygrd(j)  ,zgrd(k+1),lc1,lc3,dlen,dlen2,shift) ! 2: up   bottom  left
  s3=value(xgrd(i+1),ygrd(j)  ,zgrd(k+1),lc1,lc3,dlen,dlen2,shift) ! 3: up   bottom  right
  s4=value(xgrd(i+1),ygrd(j+1),zgrd(k+1),lc1,lc3,dlen,dlen2,shift) ! 4: up   top     right
  s5=value(xgrd(i)  ,ygrd(j+1),zgrd(k)  ,lc1,lc3,dlen,dlen2,shift) ! 5: down top     left
  s6=value(xgrd(i)  ,ygrd(j)  ,zgrd(k)  ,lc1,lc3,dlen,dlen2,shift) ! 6: down bottom  left
  s7=value(xgrd(i+1),ygrd(j)  ,zgrd(k)  ,lc1,lc3,dlen,dlen2,shift) ! 7: down bottom  right
  s8=value(xgrd(i+1),ygrd(j+1),zgrd(k)  ,lc1,lc3,dlen,dlen2,shift) ! 8: down top     right
  v(1:8)=(/s1,s2,s3,s4,s5,s6,s7,s8/)
  call SHWRITE(ifile,xgrd(i),xgrd(i+1),ygrd(j),ygrd(j+1),zgrd(k+1), zgrd(k), v)
  end do
 end do
end do
write(2,*) "};"
close(2)write(*,*) "### END GENBGMESH9_GRD END ###"
return
end subroutine genbgmesh9_grd

