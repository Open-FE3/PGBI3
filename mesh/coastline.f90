!## Modified on 2016.09.06
!## Coded by T. Minami on 2014/01/15
!## last modification date is 2014/04/25
!## May 22, 2014, zlabel is introduced to identify zero-depth nodes.
!## Modified on Oct. 24 and 25, 2015, to avoid h=0 in gebcofile, see readgebco2
!## Please search for the modification with the words "Oct. 24, 2015" or "Oct. 25, 2015"
!## This program
!## 1 extract coastline data from gebco input file, assuming there are no h=0 data in gebco file
!## 2 generate polygon data from coastline nodes, assuming the north east corner is in the ocean
program coastline
use addcorner11m
use bgelem
use constants ! added on 2016.09.03
use coastline_data
use param_mesh ! see m_param_mesh.f90, added on Sep. 7, 2016
use topo_tool ! 2018.08.27
implicit none
type(meshpara)       :: g_meshpara ! see m_param_mesh.f90
type(grid_data)      :: gebco_grd  ! see m_coastline_data.f90
type(grid_data)      :: g_grd      ! truncated grid lon,lat,alt
type(coast_data)     :: g_coast, h_coast
type(poly_data)      :: g_poly     ! see m_coastline_data.f90
type(poly_data)      :: h_poly     ! smoothened polygon data
type(poly_data)      :: i_poly     ! polygon after add corner
type(bound_data)     :: g_bound    ! ncmax, zlabel(ncmax), ibelong(ncmax)
integer(4)           :: node,ncmax,node0,lpmax,i
! node0 is the # of nodes in read file, node is the # of nodes used in calculation
integer(4),allocatable,dimension(:) :: iflgsealand
integer(4),            dimension(3) :: iflg,ilflg
!  ilflg(i) is polygon # whithin nclose land polygons
!  where i-th corner node should be included for land
real(8),dimension(3) :: xcorner,ycorner
character(70)        :: head
character(70)        :: geofile, geofileki, gebcofile
character(70)        :: coastfile="coast.dat", zlfile  ="zlabel.dat"
character(70)        :: posfile  ="bgmesh.pos",pos5file="pos5.dat"
type(bgele)          :: pos    ! 2018.08.28
! type bgele is defined in "bgele.f90", and pos stores the info regarding background mesh file

!#[0]## read parameters
CALL READMESHPARA(g_meshpara)
CALL CALCARTESIANBOUND(g_meshpara)
     gebcofile = g_meshpara%topofile
     head      = g_meshpara%head
     geofile   = head(1:len_trim(head))//".geo"
     geofileki = head(1:len_trim(head))//"ki.geo"

!#[1]## read topography file
 CALL COUNTNODE1(gebcofile,node0)
 CALL ALLOCATEGRD(gebco_grd,node0)    !allocate (lon(node0),lat(node0),h0(node0))
 CALL READGEBCO2(gebcofile,gebco_grd) ! lon [deg], lat [deg], alt [m] (upward positive) are obtained

!#[2]## truncate topography data based on given parameters
 CALL TRUNCATEGRD3(gebco_grd, g_grd, g_meshpara)
 CALL DEALLOCATEGRD(gebco_grd)        ! see m_coastline_data.f90

!#[3]## calculate cartesian coordinate XYZ ! x (eastward) [km], y (northward) [km], z(upward) [m]
 CALL LONLATTOXY4(g_grd,g_meshpara)              ! calculate xyz from lon,lat, alt
      !g_grd%xyz(3,:)=-8000.d0
      !do i=1,g_grd%node
	!  if (sqrt(g_grd%xyz(1,i)**2. + g_grd%xyz(2,i)**2.) .lt. 300.d0 ) g_grd%xyz(3,i)=200.d0
      !end do

!#[4]## output background mesh
 CALL OUTBGMESH14(g_grd,posfile,pos,g_meshpara) ! make posfile
 call outposgeo(5,pos,pos5file) ! 5 is the # for center bgmesh info, see bgele.f90
 CALL outbgmesh14_obs(g_meshpara,posfile,pos)   ! pos generation is added on 2018.08.28

     node  = g_grd%node
     ncmax = node/50 ! ncmax is the max # of points on coastline

!#[5]## make coast
 call allocatecoast(g_coast,ncmax) ! ncmax is max # of points on coastlines
 call makecoast5(g_grd,g_coast)    ! make coast nodes and output coastfile
! call coastout6(g_coast,coastfile) ! commented out 2018.11.09

      lpmax=500 ! max total number of polygons

!#[6]## compose polygon for land and oceans
 call allocatepoly(g_poly,lpmax,ncmax,lpmax) ! nclose for g_poly is unknow here
 call makepolygon7(g_coast,g_grd,g_poly)     ! lpoly0= # of polygons
! call outpolygon8(g_poly,g_coast,head,0) ! commented out on 2018.11.09

     allocate (iflgsealand(lpmax)) ! iflgsealand : 0 -> not determined, 1 -> in ocean, 2 -> in land

!#[7]## choose and exclude land lake polygons
 ! judge if polygons are in land or oceans
 call judgeocean9(iflgsealand,g_coast,g_poly,g_grd,lpmax) ! 2018.08.28
 call landpolygon10(iflgsealand,g_poly,g_coast,lpmax) ! ind_r is modified
! call outpolygon8(g_poly,g_coast,head,1) ! commented out on 2018.11.09

!#[8]## smoothen the coastlines by spline and background size information
 call smoothen10_5(g_poly,h_poly,g_coast,h_coast,g_meshpara,500.d0,pos) ! h_poly is generated
     call deallocatecoast(g_coast)
     call deallocatepoly(g_poly)
     call outpolygon12(h_poly,head,2) ! commented out 2018.11.09

!#[9]## add corner nodes and integrate nclose polygons
 call allocatebound(g_bound,ncmax)
 call addcorner11(h_poly,i_poly,g_grd,g_bound,iflg,ilflg,xcorner,ycorner) ! i_poly is generated
! call outpolygon12(i_poly,head,3) ! commented out 2018.11.09

!#[10]## output .geo file
call outpolygeo13(i_poly,geofile,geofileki,g_bound,zlfile,iflg,ilflg,xcorner,ycorner,h_poly,g_meshpara) ! g_meshpara is added on 2017.06.28

end program coastline
!##
!##########################################################   countnode1
subroutine countnode1(gebcofile,node)
implicit none
integer(4) :: node
character(70) :: gebcofile
node=0
open(1,file=gebcofile)
do while (node .ge. 0)
  read(1,*,end=99)
  node=node+1
end do
99 continue
close(1)
write(*,'(a,a,a,i10)') "# of poins of ", trim(gebcofile)," is", node  ! inserted on Oct. 24, 2015
write(*,*) "### countnode1 end! ###"
return
end subroutine countnode1


!########################################################### lonlattoxy4
subroutine lonlattoxy4(g_grd,g_meshpara)
use constants
use coastline_data
use param_mesh
implicit none
type(meshpara),    intent(in)      :: g_meshpara
type(grid_data),   intent(inout)   :: g_grd
real(8)                            :: lonorigin,latorigin
integer(4)                         :: i,node
real(8),allocatable,dimension(:,:) :: lonlatalt ! 2018.08.28

!#[0]## set
 node      = g_grd%node           ! 2108.08.28
 allocate( lonlatalt(3,node))     ! 2018.08.28
 lonlatalt = g_grd%lonlatalt      ! lonlatalt(3,g_grd%node) 2018.08.28
 lonorigin = g_meshpara%lonorigin
 latorigin = g_meshpara%latorigin

!#[1]## calculate xyz
do i=1,g_grd%node
  g_grd%xyz(1,i) = earthrad*dcos(latorigin*d2r)*(lonlatalt(1,i)-lonorigin)*d2r ! eastward [km]
  g_grd%xyz(2,i) = earthrad*(lonlatalt(2,i)-latorigin)*d2r      ! northward [km]
  g_grd%xyz(3,i) = lonlatalt(3,i)
end do

!#[2]## output topo.xyz
if ( .true. ) then
open(1,file="topo.xyz")
write(1,'(3g15.7)') (g_grd%xyz(1,i),g_grd%xyz(2,i),g_grd%xyz(3,i),i=1,g_grd%node)
close(1)
end if

write(*,*) "### lonlattoxy4 end!###"
return
end
!########################################################### makecoast5
subroutine makecoast5(g_grd,g_coast)
use coastline_data
implicit none     ! make points on coastlines
type(grid_data),       intent(in)     :: g_grd
type(coast_data),      intent(inout)  :: g_coast
integer(4)                            :: node,neast,nsouth,ncmax
integer(4)                            :: i,ii,j,jj,i1,i2,i3,nbound,ncoast
real(8),   allocatable,dimension(:)   :: x,y,h ! x (northward), y(eastward), z(upward)
real(8),   allocatable,dimension(:)   :: cx,cy
integer(4),allocatable,dimension(:,:) :: ind0
real(8),               dimension(2)   :: cxy

!#[1]## set input
neast     = g_grd%neast
nsouth    = g_grd%nsouth
node      = g_grd%node
ncmax     = g_coast%ncmax
allocate( x(node),y(node),h(node))
allocate( cx(ncmax),cy(ncmax)    )
allocate( ind0(ncmax,2)         ) ! 2018.08.28
x(1:node) = g_grd%xyz(2,1:node)    ! north    [km]
y(1:node) = g_grd%xyz(1,1:node)    ! east     [km]
h(1:node) = g_grd%xyz(3,1:node)    ! altitude [m]

jj=0
!!## NOTE THAT the top right corner of the calculation space is in the OCEAN ##!!!!!
!# right boundary
j=neast
do i=1,nsouth-1 ! downward
  ii=neast*(i-1)+j
  i1=ii;i3=ii+neast
  if (h(i1)*h(i3) .lt. 0) then
    jj=jj+1
    call interpolate(h(i1),h(i3),x(i1),x(i3),y(i1),y(i3),cxy(1:2))
!    write(2,*) cxy(1), cxy(2),ii," 2"
    cx(jj)=cxy(1)
    cy(jj)=cxy(2)
    ind0(jj,1:2)=(/ii,2/)
    end if
end do
!# bottom boundary
i=nsouth ! leftward
do j=neast-1,1,-1
 ii=neast*(i-1)+j
 i1=ii;i2=ii+1
 if (h(i1)*h(i2) .lt. 0) then
  jj=jj+1
  call interpolate(h(i1),h(i2),x(i1),x(i2),y(i1),y(i2),cxy(1:2))
  !write(2,*) cxy(1), cxy(2),ii," 1"
  cx(jj)=cxy(1)
  cy(jj)=cxy(2)
  ind0(jj,1:2)=(/ii,1/)
 end if
end do
!# left boundary
j=1
do i=nsouth-1,1,-1 ! upward
 ii=neast*(i-1)+j
 i1=ii;i3=ii+neast
 if (h(i1)*h(i3) .lt. 0) then
  jj=jj+1
  call interpolate(h(i1),h(i3),x(i1),x(i3),y(i1),y(i3),cxy(1:2))
  !write(2,*) cxy(1), cxy(2),ii," 2"
  cx(jj)=cxy(1)
  cy(jj)=cxy(2)
  ind0(jj,1:2)=(/ii,2/)
 end if
end do
!# top boundary
i=1 ! rightward
do j=1,neast-1
 ii=neast*(i-1)+j
 i1=ii;i2=ii+1
 if (h(i1)*h(i2) .lt. 0) then
  jj=jj+1
  call interpolate(h(i1),h(i2),x(i1),x(i2),y(i1),y(i2),cxy(1:2))
  !write(2,*) cxy(1), cxy(2),ii," 1"
  cx(jj)=cxy(1)
  cy(jj)=cxy(2)
  ind0(jj,1:2)=(/ii,1/)
 end if
end do
nbound=jj
!# inside the boundaries
!# "1"
do i=2,nsouth-1
 do j=1,neast-1
  ii=neast*(i-1)+j
  i1=ii;i2=ii+1
! 1
  if (h(i1)*h(i2) .lt. 0) then
   jj=jj+1
   call interpolate(h(i1),h(i2),x(i1),x(i2),y(i1),y(i2),cxy(1:2))
   !write(2,*) cxy(1), cxy(2),ii," 1"
   cx(jj)=cxy(1)
   cy(jj)=cxy(2)
   ind0(jj,1:2)=(/ii,1/)
  end if
 end do
end do
!# "1"
do i=1,nsouth-1
 do j=2,neast-1
  ii=neast*(i-1)+j
  i1=ii;i3=ii+neast
  ! 2
  if (h(i1)*h(i3) .lt. 0) then
   jj=jj+1
   call interpolate(h(i1),h(i3),x(i1),x(i3),y(i1),y(i3),cxy(1:2))
   !write(2,*) cxy(1), cxy(2),ii," 2"
   cx(jj)=cxy(1)
   cy(jj)=cxy(2)
   ind0(jj,1:2)=(/ii,2/)
  end if
 end do
end do
ncoast=jj
write(*,*) "ncoast (# of nodes on coastline is)",ncoast
if ( g_coast%ncmax .lt. ncoast +4   ) then ! +4 is required for addcorner.f90
 write(*,*) "GEGEGE! g_coast%ncmax=",g_coast%ncmax,"is smaller than ncoast+4",ncoast+4
 stop
end if
write(*,*) "### makecoast5 end! ###"
!close(2)

!#[3]## set output
 g_coast%ind=ind0
 g_coast%cxy(1,1:ncmax) = cy(1:ncmax) ! east
 g_coast%cxy(2,1:ncmax) = cx(1:ncmax) ! north
 g_coast%ncoast = ncoast
 g_coast%nbound = nbound

return
end subroutine makecoast5
!########################################################### coastout6
subroutine coastout6(g_coast,coastfile)
use coastline_data
implicit none
type(coast_data),intent(in) :: g_coast
character(70),   intent(in) :: coastfile
integer(4)                  :: i,j

! ind(i,1) : grid node # of original gebco file on left or up side
! ind(i,2) : 1 -> ind(1) indicates left node, 2 -> ind(2) is top side node
open(3,file=coastfile)
 write(3,'(2g15.7,2i8)') (g_coast%cxy(2,i),g_coast%cxy(1,i),(g_coast%ind(i,j),j=1,2),i=1,g_coast%ncoast)
close(3)
write(*,*) "### coastout6 end! ###"

return
end subroutine coastout6
!########################################################### makepolygon7
subroutine makepolygon7(g_coast,g_grd,g_poly) ! lpoly0= # of polygons
!###            make polygons and output "polygon.dat"           ###
!label(num_coastnode,1) : polygon number, default = 0, meaning not belonging to any polygon
!label(num_coastnode,2) : node number in n-th polygon nodes, which begins with 1 in each polygon
!nuwd : up wind direction, 1 -> from top or left, 2 -> from bottom or right
use coastline_data
implicit none
type(coast_data),      intent(inout)   :: g_coast ! label is calculated here
type(grid_data),       intent(in)      :: g_grd
type(poly_data),       intent(inout)   :: g_poly
integer(4)                             :: lpoly0
integer(4),allocatable,dimension(:)    :: lpoly
real(8),   allocatable,dimension(:,:)  :: xpoly,ypoly
integer(4),allocatable,dimension(:,:)  :: ind0, label, ind
real(8),   allocatable,dimension(:)    :: cx,cy
integer(4)                             :: nbound,node,neast,i,l,lpmax, ncmax, ncoast
real(8),   allocatable,dimension(:)    :: h

!#[0]## set input
 ! g_coast
 ncoast = g_coast%ncoast
 ncmax  = g_coast%ncmax
 node   = g_grd%node
 lpmax  = g_poly%lpmax
 allocate( cx(ncmax),cy(ncmax))
 allocate( ind0(ncmax,2), label(ncmax,2), ind(ncmax,2))
 allocate( h(node) )
 allocate( xpoly(lpmax,ncmax) )
 allocate( ypoly(lpmax,ncmax) )
 allocate( lpoly(lpmax)       )
 nbound = g_coast%nbound
 cx(:)  = g_coast%cxy(2,:) ! north
 cy(:)  = g_coast%cxy(1,:) ! east
 ind0   = g_coast%ind
 neast  = g_grd%neast
 h(:)   = g_grd%xyz(3,:) ! altitude [m]
 label  = 0
 lpoly0 = 0

!#[1]## not closed polygon, which touches boundaries
 do i=1,nbound
  if (label(i,1) .eq. 0) then
   call loopelement(i,lpoly0,node,neast,ncmax,ncoast,h,label,ind0,cx,cy,lpoly,xpoly,ypoly,ind,lpmax)
  end if
 end do

!#[2]## including inside the boundaries
 do i=nbound+1,ncoast
  if (label(i,1) .eq. 0) then
   call loopelement(i,lpoly0,node,neast,ncmax,ncoast,h,label,ind0,cx,cy,lpoly,xpoly,ypoly,ind,lpmax)
  end if
 end do

 !#[3]## set label to g_coast
  g_coast%label = label
  g_coast%ind_r = ind
  g_poly%lpoly0 = lpoly0
  g_poly%lpoly = lpoly
  g_poly%xypoly(1,:,:) = ypoly
  g_poly%xypoly(2,:,:) = xpoly

 write(*,*) "Total number of polygon is",lpoly0
 write(*,*) "### makepolygon7 end!###"

return
end subroutine makepolygon7
!############################################### outpolygon8
subroutine outpolygon8(g_poly,g_coast,head,inum)
use coastline_data
implicit none
type(poly_data),       intent(in)     :: g_poly
type(coast_data),      intent(in)     :: g_coast
character(70),         intent(in)     :: head
integer(4),            intent(in)     :: inum
integer(4)                            :: ncoast, nbound, ncmax
integer(4)                            :: nclose,iclose,ii,jj,ij,i,j
integer(4),allocatable,dimension(:,:) :: ind
character(1)                          :: num
character(70)                         :: polygonfile

write(num,'(i1.1)') inum
polygonfile=head(1:len_trim(head))//num//".dat"
open(1,file=polygonfile)

!#[1]## set input
  ncmax  = g_coast%ncmax
  ncoast = g_coast%ncoast
  nbound = g_coast%nbound
  allocate( ind(ncmax,2)) ! 2018.08.28
  ind    = g_coast%ind_r

!# polygon.dat read
nclose=nbound/2
write(*,*) "nclose=",nclose
ii=0
do i=1,g_poly%lpoly0
   iclose=1
   if ( i .le. nclose ) iclose=0
   !# iclose=0 -> polygon is not closed because reaching calculation boundaries
   write(1,*) i,g_poly%lpoly(i),iclose
!  if ( iclose .eq. 0 ) nclose=nclose+1
   do j=1,g_poly%lpoly(i)
     ii=ii+1      !# ii : new ncoast node number
      write(1,*) j, g_poly%xypoly(2,i,j),g_poly%xypoly(1,i,j),ind(ii,1),ind(ii,2)
   end do
end do
close(1)
!# Now the ncoast coastline nodes are ordered 
!    to be {1 2 3},{4,5,6,7,8,9}, where {} means one polygon
if (ncoast .ne. ii) then
  write(*,*) "GEGEGE!!! ncoast is not equal to ii!!ncoast=",ncoast,"ii=",ii
  stop
end if
!# polygon.dat read end!
if ( nbound/2 .ne. nclose) then
  write(*,*) "GEGEGE nbound/2 is not equal to nclose!!! nclose=",nclose,"nbound=",nbound
  stop
end if
write(*,*) "# of coastal nodes on the calculation boundaries is",nbound
write(*,*) "# of boundary polygon (nclose) is",nclose
write(*,*) "# out put",polygonfile
write(*,*) "### outpolygon8 end!###"

return

end subroutine outpolygon8
!########################################################### outpolygon12
subroutine outpolygon12(g_poly,head,inum)
use coastline_data
implicit none
type(poly_data),intent(in)  :: g_poly
character(70),  intent(in)  :: head
integer(4),     intent(in)  :: inum
integer(4)                  :: i,k,iclose
character(1)                :: num
character(70)               :: polygonfile
!# iclose : 0 -> not closed polygon, 1 -> closed polygon

!#[1]## setting
write(num,'(i1)') inum
polygonfile=head(1:len_trim(head))//num//".dat"

!#[2]## output
open(1,file=polygonfile)
do i=1,g_poly%lpoly0
 iclose=1 ! closed
 if ( i .le. g_poly%nclose) iclose=0 ! unclosed
 write(1,'(3i7)') i,g_poly%lpoly(i),iclose
 write(1,'(i7,2g15.7)')(k,g_poly%xypoly(2,i,k),g_poly%xypoly(1,i,k),k=1,g_poly%lpoly(i))
end do
close(1)

write(*,*) "### outpolygon12 end! ###"
return
end subroutine outpolygon12

!########################################################### judgeocean9
subroutine judgeocean9(iflgsealand,g_coast,g_poly,g_grd,lpmax)
!###        judge if polygons are in land or oceans        ###
! iflgsealand : 0 -> not determined, 1 -> in ocean, 2 -> in land
! iflgcoastpresence : 0 -> no coastline, i -> i-th polygon coastline exists
! iflgcoastpresence(i,1) : right side of the node, (i,2) : downside of the node
use coastline_data
implicit none
type(poly_data),       intent(inout)  :: g_poly
type(coast_data),      intent(in)     :: g_coast
type(grid_data),       intent(in)     :: g_grd
integer(4),            intent(in)     :: lpmax
integer(4),            intent(inout)  :: iflgsealand(lpmax)
integer(4)                            :: k,i,j,istart,ii,jj,ij,lpoly0,ncmax,neast,node
integer(4),allocatable,dimension(:)   :: lpoly
integer(4),allocatable,dimension(:,:) :: ind,iflgcoastpresence

!#[1]## set input
 allocate(lpoly(lpmax)) ! 2018.08.28
 node   = g_grd%node
 ncmax  = g_coast%ncmax ! 2018.08.28
 allocate( iflgcoastpresence(node,2) ) ! 2018.08.28
 allocate( ind(ncmax,2))             ! 2018.08.28
 neast  = g_grd%neast
 ind    = g_coast%ind_r
 lpoly  = g_poly%lpoly
 lpoly0 = g_poly%lpoly0

!# generate iflgcoastpresence
iflgcoastpresence(:,:)=0
k=0
do i=1,lpoly0
 do j=1,lpoly(i)
  k=k+1
  if (ind(k,2) .eq. 1) iflgcoastpresence(ind(k,1),1)=i
  if (ind(k,2) .eq. 2) iflgcoastpresence(ind(k,1),2)=i
 end do
end do
!# generate iflgsealand
iflgsealand(:)=0
j=0
do i=1,lpoly0
! j is the last coast line node in i-th polygon
j=j+lpoly(i)
k=0
istart=ind(j,1)
ij=mod(istart,neast)
if ( istart .lt. neast ) goto 86 ! top boundary
do ii=istart/neast,1,-1 ! upward
jj=neast*(ii-1)+ij
if ( iflgcoastpresence(jj,2) .ne. i .and. iflgcoastpresence(jj,2) .ne. 0) k=k+1 ! how many times coastlines are crossed
end do
86 continue
if ( ij .eq. 0) goto 87
do jj=ij,neast ! rightward
if ( iflgcoastpresence(jj,1) .ne. i .and. iflgcoastpresence(jj,1) .ne. 0) k=k+1 ! how many times coastlines are crossed
end do
87 continue
iflgsealand(i)=mod(k,2)+1
!write(*,*) "polygon #",i,"iflgsealand=",iflgsealand(i)
end do
write(*,*)"### judgeocean9 end!  ###"
return

end subroutine judgeocean9

!########################################################### landpolygon10
subroutine landpolygon10(iflgsealand,g_poly,g_coast,lpmax)
!###       exclude polygons in lands (iflgsealand=2)       ###
use coastline_data
implicit none
type(coast_data),      intent(inout)    :: g_coast
type(poly_data),       intent(inout)    :: g_poly
integer(4),            intent(in)       :: lpmax
integer(4),            intent(in)       :: iflgsealand(lpmax)
integer(4)                              :: lpoly0, lpoly0n, nboundn, ncoastn
integer(4)                              :: ncoast,nclose,nbound,ncmax
integer(4)                              :: i,ii,jj,k,kk
integer(4),allocatable,dimension(:,:)   :: ind,indn
integer(4),            dimension(lpmax) :: lpolyn,lpoly
real(8),   allocatable,dimension(:,:)   :: xpoly, ypoly, xpolyn, ypolyn

!#[1]## set input
 ncmax  = g_coast%ncmax ! 2018.08.28
 allocate( ind(ncmax,2), indn(ncmax,2))              ! 2018.08.28
 allocate( xpoly( lpmax,ncmax), ypoly( lpmax,ncmax)) ! 2018.08.28
 allocate( xpolyn(lpmax,ncmax), ypolyn(lpmax,ncmax)) ! 2018.08.28
 ncoast = g_coast%ncoast
 nbound = g_coast%nbound
 ind    = g_coast%ind_r
 lpoly0 = g_poly%lpoly0
 lpoly  = g_poly%lpoly
 xpoly  = g_poly%xypoly(2,:,:) ! north
 ypoly  = g_poly%xypoly(1,:,:) ! east

!#[2]##
nclose=nbound/2
nboundn=nbound
xpolyn(:,:)=0.d0;ypolyn(:,:)=0.d0
indn(:,:)=0;lpolyn(:)=0;ncoastn=0
ii=0;jj=0;kk=0
write(*,*) "lpoly0=",lpoly0
do i=1,lpoly0
 !write(*,*) "lpoly0=",lpoly0,"i=",i
 if (iflgsealand(i) .eq. 1) then ! polygon is in the ocean
  jj=jj+1
  !xpolyn(jj,1:ncmax)=xpoly(i,1:ncmax)
  !ypolyn(jj,1:ncmax)=ypoly(i,1:ncmax)
  xpolyn(jj,1:lpoly(i))=xpoly(i,1:lpoly(i)) ! July 15, 2015, the above 2 lines are replaced by these 2 lines
  ypolyn(jj,1:lpoly(i))=ypoly(i,1:lpoly(i))
  indn(ii+1:ii+lpoly(i),1:2)=ind(kk+1:kk+lpoly(i),1:2)
  lpolyn(jj)=lpoly(i)
  lpoly0n=jj
  ncoastn=ncoastn+lpoly(i)
  !write(*,*) "check3"
  ii=ii+lpoly(i)
 else if ( iflgsealand(i) .eq. 2 .and. i*2 .le. nbound) then
  nclose=nclose-1
  nboundn=nboundn-1
 end if
 kk=kk+lpoly(i)
end do

! set new polygon data
 g_coast%ind_r(:,:)=indn(:,:)
 g_coast%nbound=nboundn
 g_coast%ncoast=ncoastn
!
 g_poly%lpoly0        = lpoly0n
 g_poly%xypoly(2,:,:) = xpolyn(:,:) ! north
 g_poly%xypoly(1,:,:) = ypolyn(:,:) ! east
 g_poly%lpoly(:)      = lpolyn(:)
 g_poly%nclose        = nclose             ! nclose is set

write(*,*) "Total number of polygon is reduced to",lpoly0
write(*,*) "# of coastal node is reduced to",ncoast
write(*,*) "# of coastal nodes on the calculation boundaries is",nbound
write(*,*) "# of boundary polygon (nclose) is",nclose
write(*,*) "### landpolygon10 end!###"
return

end subroutine landpolygon10
!########################################################### smoothen10_5
subroutine smoothen10_5(g_poly,h_poly,g_coast,h_coast,g_meshpara,scale,pos)
use bgelem
use coastline_data
use param_mesh
implicit none
type(meshpara),        intent(in)       :: g_meshpara
type(coast_data),      intent(inout)    :: g_coast
type(poly_data),       intent(in)       :: g_poly
type(poly_data),       intent(out)      :: h_poly
type(coast_data),      intent(out)      :: h_coast
integer(4)                              :: lpoly0,nclose,nclose2, ncmax
integer(4),allocatable,dimension(:)     :: lpoly,lpoly2
integer(4),allocatable,dimension(:,:,:) :: ind2
real(8),   allocatable,dimension(:,:)   :: xpoly,ypoly,xpoly2,ypoly2
integer(4)                              :: i,j,k,lpn,lpn1,lp,lp1,lp0,lp0n
integer(4)                              :: iclose,lpoly0pre,inum
integer(4)                              :: ii,iind,ncoast
! iind is the accumulated number of nodes for ind modify
real(8)                                 :: scale,space,tn,dist,prespace,x1,y1,alpha
real(8),   allocatable,dimension(:)     :: t,X,X2,Y,Y2
real(8),   allocatable,dimension(:)     :: tnn
real(8)                                 :: sizecoastratio
type(bgele)                             :: pos
integer(4),allocatable,dimension(:,:)   :: ind

!#[1]## set input
  ncoast = g_coast%ncoast
  ncmax  = g_coast%ncmax
  nclose = g_poly%nclose
  lpoly0 = g_poly%lpoly0
  allocate( tnn(ncoast),ind(ncmax,2),ind2(nclose,2,2) )
  allocate( lpoly(lpoly0), lpoly2(lpoly0))
  allocate( xpoly( lpoly0,ncmax), ypoly( lpoly0,ncmax))
  allocate( xpoly2(lpoly0,ncmax), ypoly2(lpoly0,ncmax))
  ind    = g_coast%ind_r
  xpoly  = g_poly%xypoly(2,:,:)
  ypoly  = g_poly%xypoly(1,:,:)
  lpoly  = g_poly%lpoly
  sizecoastratio = g_meshpara%sizecoastratio

!#
! scale : length from the origin that permits sparce meshes
xpoly2(:,:)=0.d0 ; ypoly2(:,:)=0.d0 ; lpoly2(:)=0 ; ind2(:,:,:)=0.d0
lp0=lpoly0
!write(*,*) "lp0=lpoly0,",lpoly0
lp0n=0 ! new lpoly0, which excludes small islands
nclose2=nclose ! # of unclosed polygons
iind=0 ! iind is for accumulated lpoly
do i=1,lp0 ! polygon loop start!
  iind=iind+lpoly(i)
  !write(*,*) "investigate",i,"-th polygon"
  !write(*,'(a,i6,a,i6)') ("ind(1)=",ind(k,1),"ind(k,2)",ind(k,2),k=iind-lpoly(i)+1,iind)
  !do i=1,5 ! polygon loop start!
  !  write(*,*) "i=",i,"loop start! lpoly(i)=",lpoly(i)
  lp0n=lp0n+1
  lp=lpoly(i)
  lp1=lp ! if polygon is not closed, lp1 should be lp
  if ( i .gt. nclose )lp1=lp+1  ! if pooygon is closed, count the start node twice and refer it to both the start and the end
  allocate(t(lp1),Y(lp1),Y2(lp1),X(lp1),X2(lp1))
  X(1:lp)=xpoly(i,1:lp)
  Y(1:lp)=ypoly(i,1:lp)
  if (i .gt. nclose) then ! when the polygon is closed
    X(lp+1)=X(1)
    Y(lp+1)=Y(1)
  end if
  t(:)=0.d0
  t(1)=0.d0
  do j=2,lp1  ! set t as a distance
    t(j)=t(j-1)+dsqrt((X(j)-X(j-1))**2.d0+(Y(j)-Y(j-1))**2.d0)
  end do
  call SPLINE(t, X,lp1,X2) !make spline coefficients
  call SPLINE(t, Y,lp1,Y2)
  !  write(*,'(a2,f15.7,1x,a4,2f15.7)')("t=",t(j),"X,Y=",X(j),Y(j),j=1,lp1)
  lpn1=1
  tnn(1)=0.d0
  ! make tnn
  do while (tnn(lpn1) .lt. t(lp1))
    lpn1=lpn1+1
    if (ncoast .lt. lpn1) then
      write(*,*) "GEGEGE ncoast=",ncoast,"is less than lpn=",lpn
      stop
    end if
    call SPLINT(t,X,X2,lp1,tnn(lpn1-1),x1) ! interpolate
    call SPLINT(t,Y,Y2,lp1,tnn(lpn1-1),y1) ! interpolate
    call sizebgele(x1,y1,pos,space)
    tnn(lpn1)=tnn(lpn1-1)+space*sizecoastratio ! smaller node spacing on coastlines against the bgmesh lead to not adding extra nodes on coastlines
    !  if (lpn1 .eq. 2) tnn(lpn1)=space*0.7 ! in order to constrain the size of polygon in the vicinity of calculation boundaries, commented out on Oct. 25, 2015
    if (lpn1 .eq. 2 .and. i .le. nclose) tnn(lpn1)=space*0.7 ! inserted on Oct. 25, 2015
  end do
  alpha=t(lp1)/tnn(lpn1) ! shrink the length of tnn into that of t, delete comment out on Oct. 24, 2015
  tnn(1:lpn1)=alpha*tnn(1:lpn1) ! tnn is created, delete comment out on Oct. 24, 2015
  !tnn(lpn1)=t(lp1)   !    comment out on Oct. 24, 2015
  !tnn(lpn1-1)=(tnn(lpn1-2)+tnn(lpn1))/2.d0  ! comment out on Oct. 24, 2015
  lpn=lpn1                                 ! when the polygon is not closed, lpn should be lpn1
  if ( i .gt. nclose ) lpn=lpn1-1 ! when the polygon is closed, lpn should exclude the last point (identical to start node)
  lpoly2(lp0n)=lpn
  ! if (lpn .le. 3) then ! Do not add islands consisting of less than 4 points, commented out on Oct. 25, 2015
  if (lpn .le. 4) then ! Do not add islands consisting of less than 4 points, inserted on Oct. 25, 2015
    xpoly2(lp0n,1:lpn)=0.d0
    ypoly2(lp0n,1:lpn)=0.d0
    lpoly2(lp0n)=0
    lp0n=lp0n-1
    if ( i .le. nclose )   nclose2=nclose2-1 ! reduce the # of unclosed polygon
  else
    ! make xpoly2 and ypoly2
    !write(*,*) "poygon # i=",i,"alpha=",alpha
    do j=1,lpn
      call SPLINT(t,X,X2,lp1,tnn(j),xpoly2(lp0n,j)) ! interpolate
      call SPLINT(t,Y,Y2,lp1,tnn(j),ypoly2(lp0n,j))
    end do
    if ( i .le. nclose ) then !! ASSIGN ind if ( i .le. nclose)
      ! between starting and end nodes, we cannot define ind.
      ind2(lp0n,1,:)=ind(iind-lpoly(i)+1,:) ! starting node of lp0n-th polygon
      ind2(lp0n,2,:)=ind(iind,:) ! and end node of unclosed polygons, we can still assign ind()
      write(*,*) "lp0n=",lp0n
      !   write(*,*) "ind2(lp0n,1,1-2)=",ind2(lp0n,1,1),ind2(lp0n,1,2)
      !   write(*,*) "ind2(lp0n,1,1-2)=",ind2(lp0n,2,1),ind2(lp0n,2,2)
    end if
    !      write(*,*) "lpn=",lpn,"is too small to create new polygon."
  end if
  deallocate(t,X,X2,Y,Y2)
end do  ! i loop end

!#[3]## set output h_poly and h_coast
!#[3-1]## h_poly
call allocatepoly(h_poly,lp0n+1,g_poly%ncmax,nclose2) ! +1 is for when nclose = 0
write(*,*) "check3 lp0n=",lp0n
h_poly%lpoly0              = lp0n
h_poly%xypoly(2,1:lp0n,:)  = xpoly2(1:lp0n,1:g_poly%ncmax)
h_poly%xypoly(1,1:lp0n,:)  = ypoly2(1:lp0n,1:g_poly%ncmax)
h_poly%lpoly(1:lp0n)       = lpoly2(1:lp0n)
h_poly%nclose              = nclose2
h_poly%ind2(1:nclose2,:,:) = ind2(1:nclose2,:,:) ! nclose2 < nclose
lpoly0pre                  = lpoly0

!#[3-2]## h_coast
ncoast=0
do i=1,lp0n
 ncoast = ncoast + h_poly%lpoly(i)
end do
call allocatecoast(h_coast,ncoast)
h_coast%ncoast = ncoast
write(*,*) "# of points on coast is", h_coast%ncoast
ii=0
do i=1,lp0n
 do j=1,h_poly%lpoly(i)
   ii=ii+1
   h_coast%cxy(1:2,ii)=h_poly%xypoly(1:2,i,j)
 end do
end do

write(*,*) "new lpoly0 is reduced to",lp0n, "from",lpoly0pre
write(*,*) "nclose",nclose,"is reduced to nclose2",nclose2

write(*,*) "### smoothen10_5 end!###"
return
end subroutine smoothen10_5
!########################################################### outpolygeo13
!# g_meshpara is included as input on 2017.06.28
subroutine outpolygeo13(i_poly,geofile,geofileki,g_bound,zlfile,&
&                       iflg,ilflg,xcorner,ycorner,h_poly,g_meshpara)
!####           output polygon2.geo               ####
use coastline_data
use param_mesh ! see m_param_mesh.f90, added on 2017.06.28
implicit none
type(poly_data),       intent(in)     :: i_poly ! after nclose polygons are combined
type(poly_data),       intent(in)     :: h_poly ! before nclose polygons are combined
type(bound_data),      intent(in)     :: g_bound
type(meshpara),        intent(in)     :: g_meshpara ! see m_param_mesh.f90, 2017.06.28
integer(4)                            :: nclose, ncmax
integer(4),allocatable,dimension(:)   :: lpoly_origin
integer(4),allocatable,dimension(:)   :: zlabel,ibelong
integer(4)                            :: lpoly0,i,j,jj,k,is
integer(4),allocatable,dimension(:)   :: lpoly
real(8),               dimension(3)   :: xcorner,ycorner
integer(4),            dimension(3)   :: iflg,ilflg,ncorner
integer(4),allocatable,dimension(:)   :: istart,iend
!! istart(i) and iend(i) are starting and ending node # of i-th polygon, respectively.
real(8),   allocatable,dimension(:,:) :: xpoly,ypoly
! nline is the number of lines on calculation boundaries necesary
! for closing land polygons, nline4(i,1:4) store lines of i-th polygon
! 1st to 4th boundaries are left, bottom, right, and top calculation boundaries, respectively.
integer(4),allocatable,dimension(:)   :: nline
integer(4),allocatable,dimension(:,:) :: nline4
character(70)                         :: geofile,zlfile,geofileki
real(8)                               :: sizebo

!#[0]## set input
 ncmax        = g_bound%ncmax                   ! 2018.08.28
 lpoly0       = i_poly%lpoly0
 nclose       = h_poly%nclose
 allocate(      lpoly_origin(nclose)         )  ! 2018.08.28
 allocate(      zlabel(ncmax), ibelong(ncmax))  ! 2018.08.28
 allocate(      xpoly(lpoly0,ncmax)          )  ! 2018.08.28
 allocate(      ypoly(lpoly0,ncmax)          )  ! 2018.08.28
 allocate(      nline(nclose), nline4(nclose,4))! 2018.08.28
 lpoly        = i_poly%lpoly                    ! 2018.08.28
 lpoly_origin = h_poly%lpoly                    ! 2018.08.28
 xpoly        = i_poly%xypoly(1,:,:)   ! east
 ypoly        = i_poly%xypoly(2,:,:)   ! north
 sizebo       = g_meshpara%sizebo      ! 2017.06.28
 zlabel       = g_bound%zlabel
 ibelong      = g_bound%ibelong
 allocate(      istart(nclose),iend(nclose)) ! 2018.08.28

!#[1]### open files
open(1,file=geofile)
open(3,file=geofileki)
open(2,file=zlfile)
open(4,file="mshki2ocean.ctl")
jj=0
write(1,*) "lc=",sizebo,";" ! lc is the characteristic length for element size
write(3,*) "lc=",sizebo,";" ! lc is the characteristic length for element size

!# [1] ### generate points
do i=1,lpoly0
 do j=1,lpoly(i)
  jj=jj+1
  write(1,*)"Point(",jj,")={",xpoly(i,j),",",ypoly(i,j),",0.0,lc} ;"
  write(3,*)"Point(",jj,")={",xpoly(i,j),",",ypoly(i,j),",0.0,lc} ;"
  write(2,*) zlabel(jj), ibelong(jj)
 end do
end do
close(2)
write(4,'(a20,i8)') "#### nodeseageo ####",jj
!# add node necessary for land polygon
do i=1,3
 if ( ilflg(i) .ne. 0 ) then
  jj=jj+1 !　ilflg(i) is polygon # in nclose unclosed polygons which i-th corner nodes is included
  write(3,*)"Point(",jj,")={",ycorner(i),",",xcorner(i),",0.0,lc} ;"
  ncorner(i)=jj ! jj is the node # which i-th corner corresponds to
 end if
end do
write(4,'(a20,i8)') "#### nodeseageo ####",jj

!# [2] ### generate lines
jj=0;is=0
do i=1,lpoly0
is=is+lpoly(i)
do j=1,lpoly(i)-1 ! lpoly(i) does not include starting node
jj=jj+1
write(1,*) "Line(",jj,")={",jj,",",jj+1,"} ;"
write(3,*) "Line(",jj,")={",jj,",",jj+1,"} ;"
end do
jj=jj+1
write(1,*) "Line(",jj,")={",jj,",",is-lpoly(i)+1,"} ;"
write(3,*) "Line(",jj,")={",jj,",",is-lpoly(i)+1,"} ;"
end do

!# add lines necessary for land polygons
istart(:)=0;iend(:)=0
nline(:)=0
is=1 ! the north east corner node must be added in subroutine addcorner11
do i=1,nclose
 do j=1,3
  if ( iflg(j) .eq. i-1 ) is=is+1 ! node oreder is already arranged by including corner nodes for first polygon
 end do
 is=is+lpoly_origin(i)
 iend(i)=is
 istart(i)=is-lpoly_origin(i)+1 !iend(i) and istart is end and start node number in nclose unclosed polygons
 !write(*,*) "i=",i,"lpoly_origin(i)=",lpoly_origin(i),"istart=",istart(i),"iend(i)=",iend(i)
 jj=jj+1
 nline(i)=1
 nline4(i,nline(i))=jj
 write(3,*) "Line(",jj,")={",iend(i),","
 do j=1,3
  if ( ilflg(j) .eq. i)then
   write(3,*) ncorner(j),"} ;"
   jj=jj+1
   nline(i)=nline(i)+1
   nline4(i,nline(i))=jj
   write(3,*)"Line(",jj,")={",ncorner(j),","
  end if
 end do
 write(3,*) istart(i),"};"
end do

!# [3] ### generate Line loop
is=0;k=0
do i=1,lpoly0 ! line loop mainly for ocean surface
is=is+lpoly(i);k=is-lpoly(i)
write(1,*) "Line Loop(",i,")={",(k+j,",",j=1,lpoly(i)-1),k+lpoly(i),"} ;"
write(3,*) "Line Loop(",i,")={",(k+j,",",j=1,lpoly(i)-1),k+lpoly(i),"} ;"
end do
!# line loop for land polygon
jj=lpoly0
!write(*,*) "lpoly0=",lpoly0
!write(*,*) "nclose=",nclose
do i=1,nclose
jj=jj+1
write(3,*) "Line Loop(",jj,")={",(j,",",j=istart(i),iend(i)-2),iend(i)-1 ! iend - istart lines indicate i-th unclosed polygons
write(3,*) (",",nline4(i,j),j=1,nline(i)),"} ;" ! i-th unclosed polygon is closed
end do
!# [4] ### generate Surface for delauney meshing
write(1,*) "Plane Surface(1)={",(i,",",i=1,lpoly0-1),lpoly0,"} ;"
write(3,*) "Plane Surface(1)={",(i,",",i=1,lpoly0-1),lpoly0,"} ;"
do i=2,lpoly0
write(3,*) "Plane Surface(",i,")={",i,"} ;" ! island polygons
end do
do i=lpoly0+1,lpoly0+nclose
write(3,*) "Plane Surface(",i,")={",i,"} ;" ! continent, land touching boundaries, polygons
end do
!# [5] ### difine physical entity
!write(3,*) "oceansurface = 10"
!write(3,*) "Physical Surface(oceansurface)={1};"
!write(3,*) "landsurface = 20"
!write(3,*) "Physical Surface(landsurface)={2",(",",i,i=3,lpoly0+nclose), "};"
close(1)
close(3)
close(4)
write(*,*) "### outpolygeo13 end! ###"
return
end subroutine outpolygeo13
!########################################################### outbgmesh14
subroutine outbgmesh14(g_grd,posfile,pos,g_meshpara)
use bgelem
use constants
use coastline_data
use param_mesh
implicit none
type(meshpara), intent(in)  :: g_meshpara
type(grid_data),intent(in)  :: g_grd
character(70),  intent(in)  :: posfile
type(bgele),    intent(out) :: pos
real(8)                     :: wlon,elon,slat,nlat,lonorigin,latorigin,sizein,sizebo
integer(4)                  :: node,ifile
real(8)                     :: x1,y1,x2,y2,si,sb,r2,x01,x02,y01,y02

!#[0]## set
 x01  = g_meshpara%xbound(1) -20.
 x1   = g_meshpara%xbound(2)
 x2   = g_meshpara%xbound(3)
 x02  = g_meshpara%xbound(4) +20.
 y01  = g_meshpara%ybound(1) -20.
 y1   = g_meshpara%ybound(2)
 y2   = g_meshpara%ybound(3)
 y02  = g_meshpara%ybound(4) +20.
 si   = g_meshpara%sizein
 sb   = g_meshpara%sizebo
 node = g_grd%node

!#
ifile=1
r2=dsqrt(2.d0)
call initbgele(9,pos) ! see bgele.f90 2018.08.28

open(ifile,file=posfile)
write(ifile,*)'View "backgraound mesh"{'
! top left
call SQWRITE(ifile,x01,x1,y2,y02,(/sb,sb,si,sb/))
call setbgele(pos,1,x01,x1,y2,y02,(/sb,sb,si,sb/))
! top center
call SQWRITE(ifile,x1,x2,y2,y02,(/sb,si,si,sb/))
call setbgele(pos,2,x1,x2,y2,y02,(/sb,si,si,sb/))
! top right
call SQWRITE(ifile,x2,x02,y2,y02,(/sb,si,sb,sb/))
call setbgele(pos,3,x2,x02,y2,y02,(/sb,si,sb,sb/))
! middle left
call SQWRITE(ifile,x01,x1,y1,y2,(/sb,sb,si,si/))
call setbgele(pos,4,x01,x1,y1,y2,(/sb,sb,si,si/))
! middle center
call SQWRITE(ifile,x1,x2,y1,y2,(/si,si,si,si/))
call setbgele(pos,5,x1,x2,y1,y2,(/si,si,si,si/))
! middle right
call SQWRITE(ifile,x2,x02,y1,y2,(/si,si,sb,sb/))
call setbgele(pos,6,x2,x02,y1,y2,(/si,si,sb,sb/))
! bottom left
call SQWRITE(ifile,x01,x1,y01,y1,(/sb,sb,sb,si/))
call setbgele(pos,7,x01,x1,y01,y1,(/sb,sb,sb,si/))
! bottom center
call SQWRITE(ifile,x1,x2,y01,y1,(/si,sb,sb,si/))
call setbgele(pos,8,x1,x2,y01,y1,(/si,sb,sb,si/))
! bottom right
call SQWRITE(ifile,x2,x02,y01,y1,(/si,sb,sb,sb/))
call setbgele(pos,9,x2,x02,y01,y1,(/si,sb,sb,sb/))
write(ifile,*)"};"
close(ifile)
write(*,*) "### outbgmesh14 end! ###"
return
end subroutine outbgmesh14
!######################################################### SQWRITE
subroutine SQWRITE(ifile,x1,x2,y1,y2,v)
! x is the horizontal axis, y is the vertical axis
implicit none
integer(4)           :: ifile,i,j
real(8)              :: x1,x2,y1,y2
real(8),dimension(4) :: x,y,z,v

x(1:4)=(/x1,x1,x2,x2/)
y(1:4)=(/y2,y1,y1,y2/)
z(:)=0.d0
write(ifile,*) "SQ(",x(1),",",y(1),",",z(1)
write(ifile,*)  (",",x(j),",",y(j),",",z(j),j=2,4)
write(ifile,*) "){",v(1)
write(ifile,*)  (",",v(i),i=2,4)
write(ifile,*)"};"

return
end subroutine SQWRITE
!######################################################### interpolate
subroutine interpolate(h1,h2,x1,x2,y1,y2,cxy)
implicit none
real(8)              :: h1,h2,x1,x2,y1,y2,a
real(8),dimension(2) :: cxy

a=dabs(h1)/(dabs(h1)+dabs(h2))
cxy(1)=x1*(1.d0-a)+x2*a
cxy(2)=y1*(1.d0-a)+y2*a
if (cxy(2) > 1.d+4) then
write(*,*) "a,1-a",a,1.-a
write(*,*) "h1,h2=",h1,h2
write(*,*) "x1,x2=",x1,x2
write(*,*) "y1,y2=",y1,y2
write(*,*) "cxy(1,2)=", cxy(1),cxy(2)
stop
end if

return
end
!######################################################### loopelement
! iii : start node # of polygon
! label(num_coastnode,1) : polygon number, default = 0, meaning not belonging to any polygon
! label(num_coastnode,2) : node number in n-th polygon nodes, which begins with 1
! iclose : 0 -> polygon is not closed, i.e. touching boundaries, 1 -> polygon is closed
subroutine loopelement(iii,l,node,neast,ncmax,ncoast,h,label,ind0,&
&                      cx,cy,lpoly,xpoly,ypoly,ind,lpmax)
implicit none
integer(4)                                      :: i,ii,iii,l,j,k,node,ncmax,ncoast,neast
integer(4)                                      :: nuwd,nuwd2,npoly,iclose,lpmax
integer(4),dimension(ncmax,2),    intent(in)    :: ind0
integer(4),dimension(ncmax,2),    intent(inout) :: label,ind
integer(4),dimension(lpmax),      intent(inout) :: lpoly
real(8),   dimension(lpmax,ncmax),intent(inout) :: xpoly,ypoly
integer(4),dimension(ncoast)                    :: polygon
real(8),   dimension(node)                      :: h
real(8),   dimension(ncmax)                     :: cx,cy

 i=iii  ! i is the starting node
 l=l+1 ! increase polygon number
 j=1
! i-th node on coast nods
 label(i,1)=l ! set polygon number
 label(i,2)=j ! set node number in l-th polygon
 polygon(1)=i
!write(*,*) "polygon #=",l,"start!!"
! set first nuwd
 nuwd=0
 if (ind0(i,2) .eq. 1) then
    nuwd=1 ! from top
    if ( ind0(i,1) .gt. node - neast ) nuwd = 3 ! from bottom
 end if
 if (ind0(i,2) .eq. 2) then
    nuwd=2 ! from left
    if ( mod(ind0(i,1),neast) .eq. 0 ) nuwd = 4 ! from right
 end if
!   write(*,*) "ind(i,1)=",ind(i,1),"ind(i,2)=",ind(i,2),"nuwd=",nuwd
!##  polygon loop start
100 continue
 call findnextnode(i,nuwd,ncoast,ncmax,node,neast,h,ind0,ii,nuwd2)
!write(*,*) "i,ind(i,1),ind(i,2),nuwd=",i,ind(i,1),ind(i,2)
!write(*,*) "ii,ind(ii,1),ind(ii,2),nuwd2=",ii,ind(ii,1),ind(ii,2),nuwd2
!### if the circle is closed ###
if ( label(ii,1) .eq. l ) then
    !   write(*,*) "polygon is closed"
    iclose=1
    goto 200
end if
 label(ii,1)=l
 j=j+1
 label(ii,2)=j
 polygon(j)=ii ! ii is the next node of i in j-th polygon
!### if reach calculation boundaries ###
 if ( (ind0(ii,1) .le. neast .and. ind0(ii,2) .eq. 1 ) .or. &  !top boundary
    &     (ind0(ii,1) .gt. node-neast .and. ind0(ii,2) .eq. 1 ) .or. & !bottom boundary
    &     (mod(ind0(ii,1),neast) .eq. 0 .and. ind0(ii,2) .eq. 2 ) .or. & !right boundary
    &     (mod(ind0(ii,1),neast) .eq. 1 .and. ind0(ii,2) .eq. 2 )) then ! left boundary
!	write(*,*) "reach calculation boundary"
    iclose=0
   goto 200
 end if
! setting for next step
 i=ii
 nuwd=nuwd2
 goto 100
!## polygon loop end
200 continue
npoly=j ! npoly is the # of nodes which constitute l-th polygon
!write(ifile,'(3i7)') l,npoly,iclose
!write(ifile,'(i7,2g15.7,2i7)') (j,cx(polygon(j)),cy(polygon(j)),&
!&                       ind0(polygon(j),1),   ind0(polygon(j),2),j=1,npoly)
!write(*,*) "polygon #",l,"end! # of corner is",npoly
k=0
do j=1,l-1 ! l is the polygon #
   k=k+lpoly(j) ! k is summation of nodes of polygons from 1-st to l-th polygons
end do
lpoly(l)=npoly ! # of nodes in l-th polygon is npoly
do j=1,npoly ! npoly is the # of nodes which constitue l-th polygon
   xpoly(l,j)=cx(polygon(j)) ! polygon(j) is the node # of j-th node of l-th polygon
   ypoly(l,j)=cy(polygon(j))
   k=k+1
   ind(k,1:2)=ind0(polygon(j),1:2)
   ! convert ind0 to ind, which is ordered in polygon node order
end do
return
end
!
!############################################## findnextnode
! findnextnode return the next node on coastline
!# nuwd stands for up wind direction in looking for on-coastline nodes
! nuwd = 1 : from the top or left side
! nuwd = 2 : from the bottom or right side
! then we can make tabele as follows
!  nuwd    1         2              3             4
!         from top  from left      from bottom   from right
!# i is the i-th node in ncoast nodes on coastlines
!# h are altitude data from gebco original file [m]
!# ind(ncoast,2) points to left/top original gebco grid nodes
! ind(i,1) : grid node # of original gebco file on left or up side
! ind(i,2) : 1 -> ind(1) indicates left node, 2 -> ind(2) is top side node
subroutine findnextnode(i,nuwd,ncoast,ncmax,node,neast,h,ind0,ii,nuwd2)
implicit none
integer(4)                    :: i,nuwd,neast,ncoast,node,ii,j,ij,jj,nuwd2,ncmax
real(8),   dimension(node)    :: h
integer(4),dimension(ncmax,2) :: ind0
integer(4),dimension(3,2)     :: ind1
integer(4),dimension(3)       :: id3,iflag
integer(4),dimension(3,4)     :: nuwd34
data nuwd34(1,1:4) /4,3,4,3/
data nuwd34(2,1:4) /2,1,2,1/

data nuwd34(3,1:4) /1,2,3,4/
! look for on-coastline nodes on the other three edges
! determine three (ind(1), ind(2)) pairs
! set ind1
ind1(:,:)=0
if ( nuwd .eq. 1) then
ind1(1,1)=ind0(i,1)       ; ind1(1,2)=2
ind1(2,1)=ind0(i,1)+1     ; ind1(2,2)=2
ind1(3,1)=ind0(i,1)+neast ; ind1(3,2)=1
end if
if ( nuwd .eq. 2) then
ind1(1,1)=ind0(i,1)       ; ind1(1,2)=1
ind1(2,1)=ind0(i,1)+neast ; ind1(2,2)=1
ind1(3,1)=ind0(i,1)+1     ; ind1(3,2)=2
end if
if ( nuwd .eq. 3) then
ind1(1,1)=ind0(i,1)-neast   ; ind1(1,2)=2
ind1(2,1)=ind0(i,1)-neast+1 ; ind1(2,2)=2
ind1(3,1)=ind0(i,1)-neast   ; ind1(3,2)=1
end if
if ( nuwd .eq. 4) then
ind1(1,1)=ind0(i,1)-1       ; ind1(1,2)=1
ind1(2,1)=ind0(i,1)-1+neast ; ind1(2,2)=1
ind1(3,1)=ind0(i,1)-1       ; ind1(3,2)=2
end if
if ( ind1(1,1) .eq. 0 ) then
write(*,*) "GEGEGE nuwd was not set from 1 to 4!!"
stop
end if
!##       ind1 set end      ####
!
!## look for next node on coastline ##
ij=0
id3(1:3)=0
iflag(1:3)=0
do j=1,ncoast
! write(*,*) "j,ind(j,1),ind(j,2)",j,ind(j,1),ind(j,2)
 do jj=1,3
  if (ind0(j,1) .eq. ind1(jj,1) .and. ind0(j,2) .eq. ind1(jj,2)) then
  ij=ij+1
  id3(jj)=j
  iflag(jj)=1
  end if
 end do
end do
if (ij .eq. 0 ) then
write(*,*) "GEGEGE there are no next coastline nodes!! ij=",ij
write(*,*) h(ind0(i,1)-1),h(ind0(i,1)),h(ind0(i,1)+1)
write(*,*) h(ind0(i,1)-1+neast),h(ind0(i,1)+neast),h(ind0(i,1)+1+neast)
stop
end if
!## when ij = 1, 2, and 3 ##
if ( ij .eq. 1 ) then
 do j=1,3
  if ( iflag(j) .eq. 1) then
  ii=id3(j)
  nuwd2=nuwd34(j,nuwd)
  end if
 end do
end if
!## ij = 2
if ( ij .eq. 2 ) then
 write(*,*) "GEGEGE there are two next coast line nodes!! node # on coasatline is",i
 write(*,*) "nuwd=",nuwd
 if (ind0(i,1)-1-neast .ge.1) then
  write(*,*) h(ind0(i,1)-1-neast),h(ind0(i,1)-neast),h(ind0(i,1)+1-neast)
 end if
 write(*,*) h(ind0(i,1)-1),h(ind0(i,1)),h(ind0(i,1)+1)
 if (ind0(i,1)+1+neast .le. node) then
  write(*,*) h(ind0(i,1)-1+neast),h(ind0(i,1)+neast),h(ind0(i,1)+1+neast)
 end if
stop
end if
if ( ij .eq. 3 ) then
 if ( h(ind0(i,1)) > 0 ) then
 ii=id3(1)
 nuwd2=nuwd34(1,nuwd)
 else
 ii=id3(2)
 nuwd2=nuwd34(2,nuwd)
 end if
end if
return
end
!############################################################## SPLINT
SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
IMPLICIT REAL*8 (A-H,O-Z)
DIMENSION XA(N), YA(N), Y2A(N)
KLO = 1
KHI = N
1 IF ( KHI - KLO .GT. 1 ) THEN
K = ( KHI + KLO )/2
IF ( XA(K) .GT. X ) THEN
KHI = K
ELSE
KLO = K
END IF
GO TO 1
END IF
H = XA(KHI) - XA(KLO)
IF ( H .EQ. 0.D0 ) PRINT *,'BAD XA INPUT'
A = ( XA(KHI) - X )/H
B = ( X - XA(KLO) )/H
Y = A*YA(KLO) + B*YA(KHI) + ( ( A**3 - A )*Y2A(KLO) + ( B**3 - B ) &
&    *Y2A(KHI) )*(H**2)/6.D0
RETURN
END
!############################################################## SPLINE
SUBROUTINE SPLINE(X,Y,N,Y2)
IMPLICIT REAL*8 (A-H,O-Z)
!DIMENSION X(N), Y(N), Y2(N), U(2048)
DIMENSION X(N), Y(N), Y2(N), U(N)
! Interpolation by Natural Spline: Y2(1)=U(1)=QN=UN=0.d0
Y2(1) = 0.D0
U(1) = 0.D0
DO 10 I = 2, N - 1
SIG = ( X(I) - X(I-1) )/( X(I+1) - X(I-1) )
P = SIG*Y2(I-1) + 2.D0
Y2(I) = ( SIG - 1.D0 )/P
U(I) = ( 6.D0*( ( Y(I+1) - Y(I) )/( X(I+1) - X(I) ) - ( Y(I) - &
& Y(I-1) )/( X(I) - X(I-1) ) )/( X(I+1) - X(I-1) ) - SIG*U(I-1) )/P
10 CONTINUE
!
QN = 0.D0
UN = 0.D0
!
Y2(N) = ( UN - QN*U(N-1) )/( QN*Y2(N-1) + 1.D0 )
!
DO 20 K = N - 1, 1, -1
Y2(K) = Y2(K)*Y2(K+1) + U(K)
20 CONTINUE
!
RETURN
END



