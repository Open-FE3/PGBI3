!# adopt "precision.inc" for len_file on 2016.11.15
!# 
!# grid_data_2D is added on May 16, 2016
!#
!# Coded on Oct. 25, 2015
!# to calculate (vx,vy,vz) from COMCOT results
!# Contents :
!#   subroutine  COUNTCOMCOT
!#   subroutine  READCOMCOTGRD
!#   subroutine  READGRD
!#   subroutine  MKVXYZCOEF
!#   subroutine  READMN
!#   subroutine CALOCEANVXYZ
!#   subroutine CALOCEANFXYZ
!#  subroutine GENLVEC3COM(Fxyz, nodes, g_node, nodtot, Fxyzip, ip)
!#  subroutine OUTVXYZ
!########################################   COUNTCOMCOT
module FROMCOMCOT
use param
implicit none

type grid_xy
 integer(4) :: nx
 integer(4) :: ny
 real(8),allocatable,dimension(:)   :: xgrd,ygrd
end type grid_xy

type grid_data_2D
 type(grid_xy) :: xygrid
 real(8),allocatable,dimension(:,:,:) :: dat
 ! dat(i,j,k) is
 ! i-th component at coordinate (j,k) in nx*ny grid
end type grid_data_2D

!===
contains
!########################################   COUNTCOMCOT
subroutine COUNTCOMCOT(xygrid, xgrdfile, ygrdfile)
implicit none
type(grid_xy),intent(out) :: xygrid
character(70),intent(in) :: xgrdfile, ygrdfile
integer(4) :: nx,ny
real(8) :: tmp
nx=0 ; ny=0
open(1,file=xgrdfile)
do
  read(1,*,end=10) tmp
  nx=nx+1
end do
10 write(*,*) "nx=", nx, "file=",xgrdfile
close(1)
open(2,file=ygrdfile)
do
  read(2,*,end=20) tmp
  ny=ny+1
end do
20 write(*,*) "ny=", ny, "file=",ygrdfile
close(2)
xygrid%nx=nx
xygrid%ny=ny
write(*,*) "### COUNTCOMCOT END!! ###"
return
end
!########################################   READCOMCOTGRD
subroutine READCOMCOTGRD(xygrid, xgrdfile, ygrdfile)
implicit none
type(grid_xy),intent(inout) :: xygrid
character(70),intent(in) :: xgrdfile, ygrdfile
integer(4) :: i
! xgrd(1:nx) read
allocate(xygrid%xgrd(xygrid%nx))
allocate(xygrid%ygrd(xygrid%ny))
open(1,file=xgrdfile)
do i=1,xygrid%nx
  read(1,*) xygrid%xgrd(i)
end do
close(1)
! ygrd(1:ny) read
open(2,file=ygrdfile)
do i=1,xygrid%ny
  read(2,*) xygrid%ygrd(i)
end do
close(2)
write(*,*) "### READCOMCOTGRD END!! ###"
return
end

!######################################## ALLOCATEGRDDATA
subroutine ALLOCATEGRDDATA(h_grd,xygrid,idim)
implicit none
integer(4),intent(in) :: idim
type(grid_xy),intent(in) :: xygrid
type(grid_data_2D),intent(inout) :: h_grd

h_grd%xygrid%nx=xygrid%nx
h_grd%xygrid%ny=xygrid%ny
allocate(h_grd%xygrid%xgrd(xygrid%nx))
allocate(h_grd%xygrid%ygrd(xygrid%ny))
h_grd%xygrid%xgrd=xygrid%xgrd
h_grd%xygrid%ygrd=xygrid%ygrd
allocate(h_grd%dat(idim,xygrid%nx,xygrid%ny))

return
end
!######################################## ALLOCATEGRDDATA
subroutine DEALLOCATEGRDDATA(h_grd)
implicit none
type(grid_data_2D),intent(inout) :: h_grd

deallocate(h_grd%xygrid%xgrd)
deallocate(h_grd%xygrid%ygrd)
deallocate(h_grd%dat)

return
end
!
!########################################   READFXYZ
subroutine READGRDDATA(infile,h_grd,ncomp)
implicit none
type(grid_data_2D),intent(inout) :: h_grd
integer(4),        intent(in)    :: ncomp
character(70),     intent(in)    :: infile
integer(4) :: i,j,k
real(8) :: f
!#[1]## Fxyz read
write(*,*) "infile=",infile
open(1,file=infile)
do j=1,h_grd%xygrid%ny
  do i=1,h_grd%xygrid%nx
     read(1,*) ( h_grd%dat(k,i,j),k=1,ncomp ) !,f
!    read(1,*) Fxyz_grd%dat(2,i,j), Fxyz_grd%dat(1,i,j), Fxyz_grd%dat(3,i,j),f
!    write(*,*) Fxyz_grd%dat(2,i,j), Fxyz_grd%dat(1,i,j), Fxyz_grd%dat(3,i,j),f
  end do
end do
close(1)
!Fxyz_grd%dat(3,:,:)=-Fxyz_grd%dat(3,:,:)

write(*,*) "### READGRDDATA END!! ###"
return
end
!##########################################################
subroutine READMN(it, vxyh_grd, g_param)  ! inserted on Oct. 25, 2015
implicit none
integer(4),         intent(in)    :: it
type(param_forward),intent(in)    :: g_param
type(grid_data_2D), intent(inout) :: vxyh_grd
character(70)                     :: mpre, npre ! e.g., "m_01_"
character(70)                     :: mfile, nfile
integer(4)                        :: itime_comcot, i, j,nx,ny
character(6)                      :: num
real(8)                           :: dt_comcot ! 2018.11.15

!#[0]## set input
  mpre      = g_param%mpre
  npre      = g_param%npre
  dt_comcot = g_param%dt_comcot ! 2018.11.15

!#[1]## generate file name
  itime_comcot = it*dt_comcot  ! 2018.11.15
  nx           = vxyh_grd%xygrid%nx
  ny           = vxyh_grd%xygrid%ny
  write(num,'(i6.6)') itime_comcot
  mfile = mpre(1:len_trim(mpre))//num(1:6)//".dat"
  nfile = npre(1:len_trim(mpre))//num(1:6)//".dat"
  write(*,*) "mfile=",mfile
  write(*,*) "nfile=",nfile

!#[2]## read vxh, vyh
open(1,file=mfile)
open(2,file=nfile)
do j=1,ny
!  write(*,*) "j=",j,"ny=",ny
  read(1,'(15f9.3)') (vxyh_grd%dat(1,i,j),i=1,nx)
  read(2,'(15f9.3)') (vxyh_grd%dat(2,i,j),i=1,nx)
end do
close(1)
close(2)
write(*,*) "### READMN END!! ###"
return
end
!
!
!#################################################
subroutine GENLVEC3COM(Fxyz, nodes, g_node, nodtot, Fxyzip, ip)
implicit none
integer(4),intent(in) :: nodes, ip
real(8),intent(in) :: Fxyz(3, nodes)
integer(4),intent(in) :: nodtot, g_node(nodtot)
real(8), intent(out) :: Fxyzip(3, nodtot)
integer(4) :: i
do i=1,nodtot
  if ( g_node(i) .le. nodes) then
    Fxyzip(1:3,i) = Fxyz(1:3, g_node(i))
  else
    Fxyzip(1:3,i)=0.d0
  end if
end do
write(*,*) "### GENLVEC3COM END!! ### ip=",ip
return
end

end module FROMCOMCOT
