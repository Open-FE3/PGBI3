! Coded by Takuto Minami on 2016.09.02
!
!(1) xyz : cartesian coordinate where x eastward, y northward, z upward
!       origin is at groudnsurface[km],[km],[km]
!
!(2) lon, lat, alt : londitude [deg], latitude [deg], altitude [km]
! 
!(3) xyzspehrical : cartesian coordinate where z southpole to northpole,
!                   x is set so that x,z plane is Greenwich meridian
!                   y axis points londitude = 90 deg
!                   origin is set at the center of the Earth
module spherical
use constants ! earthrad, pi, d2r, r2d
use outerinnerproduct ! see m_outerinnerproduct.f90
use mesh_type ! see m_mesh_type.f90
use param_mesh ! see m_param_mesh.f90
implicit none

real(8) :: colatorigin,lon,lat,alt,x,y,z,r,colat  ! [km]

contains

!########################################### tensorsphe2xyz
subroutine tensorsphe2xyz(lonlatalt,sphe2xyz)
implicit none
real(8),intent(in)  :: lonlatalt(3)
real(8),intent(out) :: sphe2xyz(3,3)

lon   = lonlatalt(1)*d2r
colat = (90.d0 - lonlatalt(2))*d2r

!# tensor for spherical to xyz
sphe2xyz(1,1) = -sin(lon)
sphe2xyz(1,2) =  cos(lon)
sphe2xyz(1,3) =  0.d0
sphe2xyz(2,1) = -cos(lon)*cos(colat)
sphe2xyz(2,2) = -sin(lon)*cos(colat)
sphe2xyz(2,3) =           sin(colat)
sphe2xyz(3,1) =  cos(lon)*sin(colat)
sphe2xyz(3,2) =  sin(lon)*sin(colat)
sphe2xyz(3,3) =           cos(colat)

return
end
!########################################### subroutine
! coded on 2016.11.19
subroutine meshxyz2xyzspherical(h_mesh,g_meshpara)
implicit none
type(mesh),intent(inout) :: h_mesh
type(meshpara),intent(in) :: g_meshpara
real(8) :: lonorigin,latorigin
integer(4) :: i

!#[1]##
 lonorigin = g_meshpara%lonorigin
 latorigin = g_meshpara%latorigin
 allocate(h_mesh%lonlatalt(3,h_mesh%node))
 allocate(h_mesh%xyzspherical(3,h_mesh%node))

!#[2]##
do i=1,h_mesh%node
 call xyz2lonlatalt(h_mesh%xyz(1:3,i),lonorigin,latorigin,h_mesh%lonlatalt(1:3,i))
 call lonlatalt2xyzspherical(h_mesh%lonlatalt(1:3,i),h_mesh%xyzspherical(1:3,i))
! write(*,*) "h_mesh%xyz(1:3,i)=",h_mesh%xyz(1:3,i)
! write(*,*) "->",h_mesh%xyzspherical(1:3,i)
end do
write(*,*) "### meshxyz2xyzspherical end!! ###"

return
end
!########################################### subroutine sphrcl2decar
! coded on 2016.11.19
subroutine vectorspherical2xyz(lonlatalt,vector,node)
implicit none
real(8),intent(in)    :: lonlatalt(3,node)    ! position of vector
real(8),intent(inout) :: vector(3,node) ! vector in xyz -> vector in spherical xyz
integer(4),intent(in) :: node
real(8) :: vec(3),v(3,node)
integer(4) :: i

!#[1]## calculate v
do i=1, node
 !# rotate vector
 vec(1:3) = vector(1:3,i)
 lon    = lonlatalt(1,i)*d2r
 colat  = (90.d0 - lonlatalt(2,i))*d2r
 v(1,i)=-vec(1)*sin(lon)            + vec(2)*cos(lon)
 v(2,i)=-vec(1)*cos(lon)*cos(colat) - vec(2)*sin(lon)*cos(colat) + vec(3)*sin(colat)
 v(3,i)= vec(1)*cos(lon)*sin(colat) + vec(2)*sin(lon)*sin(colat) + vec(3)*cos(colat)
end do

!#[2]## set output
vector = v

return
end

!########################################### subroutine sphrcl2decar
! coded on 2016.11.19
subroutine vectorxyz2spherical(lonlatalt,vector,node)
implicit none
real(8),            intent(in)     :: lonlatalt(3,node)    ! position of vector
real(8),            intent(inout)  :: vector(3,node)!vector in xyz-> vector in spherical xyz
integer(4),         intent(in)     :: node
real(8),allocatable,dimension(:,:) :: v ! 2018.08.28
real(8)    :: vec(3)
integer(4) :: i

allocate(v(3,node)) ! 2018.08.28

!#[1]## calculate v
do i=1, node
 !# rotate vector
 vec(1:3) = vector(1:3,i)
 lon    = lonlatalt(1,i)*d2r
 colat  = (90.d0 - lonlatalt(2,i))*d2r
 v(1,i) = -vec(1)*sin(lon) - vec(2)*cos(colat)*cos(lon) + vec(3)*sin(colat)*cos(lon)
 v(2,i) = vec(1)*cos(lon) - vec(2)*cos(colat)*sin(lon) + vec(3)*sin(colat)*sin(lon)
 v(3,i) =                   vec(2)*sin(colat)          + vec(3)*cos(colat)
end do

!#[2]## set output
vector = v

return
end


!########################################### subroutine xyz2lonlatalt
subroutine xyz2lonlatalt(xyz,lonorigin,latorigin,lonlatalt)
implicit none
real(8),intent(in) :: xyz(3),lonorigin,latorigin
real(8),intent(out) :: lonlatalt(3)
colatorigin = 90.d0 - latorigin
x=xyz(1) ; y=xyz(2) ; z=xyz(3)

lon = x/d2r/earthrad/sin(colatorigin*d2r) + lonorigin
lat = y/d2r/earthrad + latorigin
alt = z ! [km]

lonlatalt(1:3)=(/lon,lat,alt/)

return
end subroutine xyz2lonlatalt

!########################################### subroutine lonlatalt2xyz
subroutine lonlatalt2xyz(lonlatalt,lonorigin,latorigin,xyz)
implicit none
real(8),intent(in) :: lonlatalt(3)
real(8),intent(in) :: lonorigin,latorigin
real(8),intent(out) :: xyz(3)
lon=lonlatalt(1)
lat=lonlatalt(2)
alt=lonlatalt(3)
colatorigin = 90.d0 - latorigin

x = (lon - lonorigin)*d2r*earthrad*sin(colatorigin*d2r) ! [km]
y = (lat - latorigin)*d2r*earthrad ! [km]
z = alt ! [km]

xyz(1:3)=(/x,y,z/)

return
end subroutine lonlatalt2xyz

!########################################### subroutine decar2sphrcl
subroutine lonlatalt2xyzspherical(lonlatalt,xyzspherical)
implicit none
real(8),intent(in) :: lonlatalt(3) ! londitude [deg],latitutde[deg], altitude [km]
real(8),intent(out) :: xyzspherical(3) ! cartesian coordinate [km]
lon =lonlatalt(1) ; lat=lonlatalt(2) ; alt=lonlatalt(3)
colat = 90.d0 - lat ! [deg]

r = earthrad + alt ! [km]
x = r*sin(colat*d2r)*cos(lon*d2r) ! [km]
y = r*sin(colat*d2r)*sin(lon*d2r) ! [km]
z = r*cos(colat*d2r)              ! [km]

xyzspherical(1:3)=(/x,y,z/)

return
end subroutine lonlatalt2xyzspherical

!########################################### subroutine sphrcl2decar
subroutine xyzspherical2lonlatalt(xyzspherical,lonlatalt)
implicit none
real(8),intent(in) :: xyzspherical(3)
! x,y,z corrdinate [km] when the origin is at the center of the earth
real(8),intent(out) :: lonlatalt(3)
x=xyzspherical(1) ! [km]
y=xyzspherical(2) ! [km]
z=xyzspherical(3) ! [km]

r=sqrt(x**2. + y**2. + z**2.) ! [km]
lon = atan2(y,x)*r2d          ! [deg]
lat = asin(z/r)*r2d           ! [deg]
alt = r - earthrad            ! [km]

lonlatalt(1:3)=(/lon,lat,alt/)

return
end subroutine xyzspherical2lonlatalt


end module spherical
