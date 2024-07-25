!# modified the format of input ctl file on 2017.06.27
!# Coded on September 7, 2016
module param_mesh
use constants
implicit none

type meshpara
!# topography file
 character(50) :: topofile
!# coordinate info
 real(8) :: wlon,elon,slat,nlat
 real(8) :: lonorigin,latorigin
 real(8) :: zmin,zmax
 real(8) :: xbound(4)     !
 real(8) :: ybound(4) !
 !  ----------- y02
 !  |  -----  | y2
 !  |  |   |  |
 !  |  -----  | y1
 !  ----------- y01
 !x01 x1  x2 x02
!# resolution info
 real(8) :: sizein,sizebo,sizecoastratio,lenout
 real(8) :: upzin       ! [km]
 real(8) :: downzin     ! [km]
 real(8) :: sizein3d    ! [km]
!# info of observatory
 integer(4) :: nobs
 character(5),allocatable,dimension(:) :: obsname
 real(8),allocatable,dimension(:,:) :: lonlatalt ! lon[deg],lat[deg], alt [m]
 real(8),allocatable,dimension(:,:) :: xyz ! east, north, upward [km]
!# obs mesh info
 real(8) :: A     ! [km] resolution at observatories
 real(8) :: sigma ! [km] sigma for gaussian resolution near observatories
 integer(4) :: iflag_extrude ! 0 for nothing, 1: vertical resolution is less than hori
 integer(4) :: nlayer_max    ! 2017.06.27
 real(8),allocatable,dimension(:) :: depth_nlayer   ! depth_nlayer(nlayer_max-1) 2017.06.27
!# header
 character(50) :: head    ! header like "polygon"
end type meshpara

contains
!######################################### calobsxyz
subroutine calobsxyz(g_meshpara)
use constants
implicit none
type(meshpara),intent(inout) :: g_meshpara
integer(4) :: i,nobs
real(8) :: lonlatalt(3,g_meshpara%nobs),lonorigin,latorigin

!#[1]## set input
lonorigin = g_meshpara%lonorigin
latorigin = g_meshpara%latorigin
nobs      = g_meshpara%nobs
lonlatalt = g_meshpara%lonlatalt ! lon [deg], lat [deg], alt [m]

!#[2]## calculate xyz for obs
do i=1,nobs
 g_meshpara%xyz(1,i) = (lonlatalt(1,i) - lonorigin )*d2r*earthrad*cos(latorigin * d2r)![km]
 g_meshpara%xyz(2,i) = (lonlatalt(2,i) - latorigin )*d2r*earthrad ! [km]
 g_meshpara%xyz(3,i) =  lonlatalt(3,i)*1.d-3  ! [m] -> [km]
end do

return
end
!######################################### calcartesianbound
subroutine calcartesianbound(g_meshpara)
type(meshpara),intent(inout) :: g_meshpara
real(8) :: lonorigin,latorigin,nlat,elon,lenout
real(8) :: x1,x2,y1,y2     !
real(8) :: x01,x02,y01,y02
!#[1]## set input
lonorigin = g_meshpara%lonorigin
latorigin = g_meshpara%latorigin
elon = g_meshpara%elon
nlat = g_meshpara%nlat
lenout = g_meshpara%lenout

y2=earthrad*d2r*(nlat-latorigin)
y1=-y2
x2=earthrad*dcos(latorigin*d2r)*d2r*(elon-lonorigin)
x1=-x2

x01=x1-lenout
x02=x2+lenout
y01=y1-lenout
y02=y2+lenout

g_meshpara%xbound(1:4) = (/x01,x1,x2,x02/)
g_meshpara%ybound(1:4) = (/y01,y1,y2,y02/)

return
end subroutine calcartesianbound

!######################################### readmeshpara
subroutine readmeshpara(g_meshpara,parafile)
type(meshpara),intent(out) :: g_meshpara
character(50),optional     :: parafile    ! added on 2016.09.15
integer(4) :: i,n,input = 5

if ( present(parafile) ) input = 1
if ( present(parafile) ) open(input,file=parafile)

!#[0]## skip header of controll file
write(*,*) "input # of header lines"
read(input,20) n
write(*,*) "# of header is",n
do i=1,n
 read(input,*)
end do

!#[1]## topography file
write(*,*) "input filename of topgraphyfile (lon[deg],lat[deg],alt[m])"
read(input,21) g_meshpara%topofile
write(*,*) g_meshpara%topofile

!#[2]## read wlon,elon,slat,nlat / zmin, zmax
write(*,*) "input wlon,elon,slat,nlat [deg] line by line"
read(input,22) g_meshpara%wlon
read(input,22) g_meshpara%elon
read(input,22) g_meshpara%slat
read(input,22) g_meshpara%nlat
write(*,*) "input zmin, zmax [km] (upwardpositive)"
read(input,22) g_meshpara%zmin
read(input,22) g_meshpara%zmax


!#[3]## calculate lonorigin and latorigin
g_meshpara%lonorigin = (g_meshpara%wlon + g_meshpara%elon)/2.
g_meshpara%latorigin = (g_meshpara%slat + g_meshpara%nlat)/2.
write(*,*) "lonorigin=",g_meshpara%lonorigin
write(*,*) "latorigin=",g_meshpara%latorigin

!#[4]## resolution info
write(*,*) "input sizein [resolution [km] in focus area]"
read(input,22) g_meshpara%sizein
write(*,*) "input sizebo [resolution [km] at calculation boundary]"
read(input,22) g_meshpara%sizebo
write(*,*) "input sizecoastratio [resolution ratio to sizein at coastline]"
read(input,22) g_meshpara%sizecoastratio
write(*,*) "input lenout [ditance [km] from the edge of focus area to boundary]"
read(input,22) g_meshpara%lenout

!#[5]## observatory info
write(*,*) "input nobs"
read(input,20) g_meshpara%nobs
allocate(g_meshpara%obsname(g_meshpara%nobs))
allocate(g_meshpara%lonlatalt(3,g_meshpara%nobs))
allocate(g_meshpara%xyz(3,g_meshpara%nobs))
write(*,*) "input observatory name in less than 6 characters"
write(*,*) "input observatory lon [deg], lat[deg], altitude [m]"
do i=1,g_meshpara%nobs
 read(input,'(a)') g_meshpara%obsname(i)
 read(input,23) g_meshpara%lonlatalt(1:3,i)
end do
write(*,*) "input observatory resolution A [km] and sigma [km]"
read(input,22) g_meshpara%A
read(input,22) g_meshpara%sigma
write(*,*) "input iflag_extrude: 0 for nothing, 1 for vertical resolution less than horizontal [km]"
read(input,20) g_meshpara%iflag_extrude
write(*,*) "input nlayer_max and depth_nlayer (1:nlayer_max-1)"
read(input,20) g_meshpara%nlayer_max    !
allocate(g_meshpara%depth_nlayer(g_meshpara%nlayer_max - 1))
do i=1,g_meshpara%nlayer_max - 1
 read(input,22) g_meshpara%depth_nlayer(i) ![km]
end do
call calobsxyz(g_meshpara)

!#[6]## read header
write(*,*) "input header"
read(input,21) g_meshpara%head

!#[7]## bgmesh resolution for
write(*,*) "input upzin, downzin, sizein3d"

 read(input,22)  g_meshpara%upzin    ! [km]
 read(input,22)  g_meshpara%downzin  ! [km]
 read(input,22)  g_meshpara%sizein3d ! [km]

20 format(20x,i10  )
21 format(20x,a50  )
22 format(20x,g15.7)
23 format(20x,3g15.7)

if ( present(parafile) ) close(input)

end subroutine readmeshpara

end module param_mesh
