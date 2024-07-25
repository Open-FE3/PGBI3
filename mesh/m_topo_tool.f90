!# coded on 2018.08.27
module topo_tool
implicit none

contains

!##########################################################   readgebco2
subroutine readgebco2(gebcofile,gebco_grd)
use coastline_data
implicit none
character(70),  intent(in)    :: gebcofile
type(grid_data),intent(inout) :: gebco_grd
integer(4)                    :: i
real(8)                       :: alt

!#[1]## read grd data
open(2,file=gebcofile)
 do i=1, gebco_grd%node ! node is set
  read(2,*) gebco_grd%lonlatalt(1:2,i),alt ! alt [m] (upward positive)
  gebco_grd%lonlatalt(3,i) = alt
  if ( alt .eq. 0.d0 ) gebco_grd%lonlatalt(3,i) = 1.d0 ! inserted by T.M. on Oct. 24, 2015
 end do
close(2)

!#[2]## calculate neast and nsouth
do i=1,gebco_grd%node
 if (gebco_grd%lonlatalt(2,i) .ne. gebco_grd%lonlatalt(2,i+1)) goto 99    ! latitude
end do
99 continue
gebco_grd%neast=i
gebco_grd%nsouth=gebco_grd%node/gebco_grd%neast ! neast and nsout are set for h_grd

write(*,*) "### readgebco2 end! ###"
return
end subroutine readgebco2
!########################################################## neasatsouth3
subroutine TRUNCATEGRD3(h_grd,g_grd,g_meshpara)
use constants
use coastline_data
use param_mesh
! ne1 is # of node read from west, ns1 is that from north
implicit none
type(meshpara), intent(in)  :: g_meshpara
type(grid_data),intent(in)  :: h_grd ! whole grid
type(grid_data),intent(out) :: g_grd ! truncated grid
integer(4) :: node, ne1, ne2, ns1, ns2
integer(4) :: i,j,ii,jj
real(8)    :: wlon,elon,slat,nlat,lenout,lonorigin,latorigin
real(8)    :: wlon0,elon0,slat0,nlat0
real(8),allocatable,dimension(:,:) :: lonlatalt ! 2018.08.28

!#[0]## set parameter
lonlatalt = h_grd%lonlatalt  ! lonlatalt(3,h_grd%node)
wlon      = g_meshpara%wlon
elon      = g_meshpara%elon
slat      = g_meshpara%slat
nlat      = g_meshpara%nlat
lenout    = g_meshpara%lenout
lonorigin = g_meshpara%lonorigin
latorigin = g_meshpara%latorigin

!#[1]## calculate wlon0,elon0,slat0,nlat0
 wlon0 = wlon-lenout/(earthrad*d2r*dcos(latorigin*d2r))
 elon0 = elon+lenout/(earthrad*d2r*dcos(latorigin*d2r))
 slat0 = slat-lenout/(earthrad*d2r)
 nlat0 = nlat+lenout/(earthrad*d2r)

!#[2]## create ne1,ne2
!#      wlon0           elon0
!#        |               |
!#  |   | <---- ne2 ------> |            |
!#  1  ne1              (ne1-1)+ne2    neast
!#      1                   ne2
ne1=0;ne2=0
do i=1,h_grd%neast
 if (lonlatalt(1,i) .lt. wlon0) ne1=i              ! start longitude
 if (lonlatalt(1,i) .lt. elon0) ne2=i+1-(ne1-1)  ! end longitude
end do

!#[3]## calculate ns1 and ns2
ns1=0;ns2=0
do i=1,h_grd%nsouth
 if (lonlatalt(2,(i-1)*h_grd%neast +1) .gt. nlat0 ) ns1=i
 if (lonlatalt(2,(i-1)*h_grd%neast +1) .gt. slat0 ) ns2=i+1-(ns1-1)
end do

!#[4]## generate new grd
node=ne2*ns2
CALL ALLOCATEGRD(g_grd,node)
g_grd%neast=ne2
g_grd%nsouth=ns2
jj=0
do j=ns1, (ns1-1)+ns2
 do i=ne1,(ne1-1)+ne2
  ii= (j-1)*h_grd%neast + i
  jj=jj+1
  g_grd%lonlatalt(1:3,jj)=h_grd%lonlatalt(1:3,ii)
  g_grd%xyz(1:3,jj)      =h_grd%xyz(1:3,ii)
 end do
end do
write(*,*) "ns1=",ns1,"ne1=",ne1
write(*,*) "ns2=",ns2,"ne2=",ne2,"truncated node=",g_grd%node


write(*,*) "### TRUNCATEGRD3 end!###"

return
end subroutine TRUNCATEGRD3
end module
