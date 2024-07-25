program extrude
!# calztopo is modified on 2017.06.27
!## This program is coded by T.MINAMI on May 21, 2014
!## to generate 3D mesh by extruding 2-D mesh, polygon.msh
!## modified to apply minimum depth along coastlines on July 18, 2014
!## element node order correction based on volume is applied on July 22, 2014
use constants ! added on 2016.09.03
use param_mesh
use mesh_type
implicit none
type(meshpara)                        :: g_meshpara
type(mesh)                            :: h_mesh
integer(4)                            :: node,ecount,inib,lcount,ns3,neles
integer(4),allocatable,dimension(:)   :: nz ! # of vertical nodes for each horizontal column
real(8),   allocatable,dimension(:)   :: x,y,z
integer(4),allocatable,dimension(:)   :: zlabel,ibelflg,nzptr,nelecount,nhnptr,lbelong
integer(4),allocatable,dimension(:,:) :: neleptr, n2,n3
real(8),   allocatable,dimension(:)   :: zs3
real(8),               parameter      :: coastlinedepth = 0.d0 ! only 0 is accepted in this program on 2016.09.07
! file
character ( 70 )                      :: head,mshfile,polyzfile,allnodefile,gebcofile
character ( 70 )                      :: zlfile="zlabel.dat"
character ( 70 )                      :: oceanfile="ocean.msh"
! zlabel(i) = 0 -> node is on coastlines, 1 -> node is not coastline, necessary to calculate depth
! ibelflg(i) = 0 : ith initial line  is on coastlines, no need to calculate depths
!                    1 : ith initial line is not on coastline, it is necessary to calculate depth
! lbelong(j)  : the initial line number that j-th line belongs to, which is same as the number of the starting initial node in .geo file
! ns3 : # of nodes in the ocean, necessary for tsunami hydrodynamic simulations

!#[1]## read parameters
 CALL READMESHPARA(g_meshpara)
 CALL CALCARTESIANBOUND(g_meshpara)
      gebcofile=g_meshpara%topofile
      head=g_meshpara%head
      mshfile=head(1:len_trim(head))//"_ki.msh"
      polyzfile=head(1:len_trim(head))//"z.msh"
      allnodefile=head(1:len_trim(head))//"allnode.geo"

!#[2]## read zlabel.dat
 call countinib01(zlfile,inib)
     allocate(zlabel(inib),ibelflg(inib)) ! inib is the # of nodes in .geo file, namely # of initial boundary nodes
 call readzlabel02(zlfile,inib,zlabel,ibelflg)

!#[3]## read head_ki.msh
 call READMESH_TOTAL(h_mesh,mshfile)
      node    = h_mesh%node
      ecount  = h_mesh%ntri
      lcount  = h_mesh%nlin
      allocate(x(node),y(node),z(node),nz(node),n3(ecount,3),n2(lcount,2),lbelong(lcount))
      x       = h_mesh%xyz(1,:)
      y       = h_mesh%xyz(2,:)
      n3      = h_mesh%n3
      n2      = h_mesh%n2
      lbelong = h_mesh%n2flag(:,2)

call calztopo2(gebcofile, h_mesh, zlabel,ibelflg,inib,coastlinedepth,g_meshpara)
     deallocate (zlabel,ibelflg)
     z = h_mesh%xyz(3,:)
!     z = z * 30.d0 ! only for visualize

call polygonz3(mshfile,polyzfile,z,node)
call countznode4(ns3,nz,h_mesh,g_meshpara,node) !!ns3 is the number of nodes in the ocean 2018.08.28
allocate(zs3(ns3), nhnptr(ns3))  ! nhnptr(i) indicates node column # which the i-th oceanic node belongs to
call mkznode5(zs3,z,ns3,node,nz,nhnptr,node) ! nz(i) is the number of nodes in i-th column
call polygonallnode6(allnodefile,x,y,node,zs3,ns3,nz,node)
allocate (nzptr(ns3)) ! nzptr(i) : node number of i-th deep node
call sortzorder7(zs3,ns3,nzptr)
allocate (neleptr(node,10),nelecount(node))        ! neles is the total # of 3D ocean mesh, maybe 10 is maximum in delauney mesh generating
call countele8(n3(1:ecount,1:3),ecount,node,neleptr,nelecount,neles,ns3,nz(1:node)) ! neleptr(i:node,1:7) indicates which element column i-th horizontal node is included
call maketetra9(n3(1:ecount,1:3),ecount,node,nzptr,neleptr,nhnptr,nelecount,ns3,nz(1:node),neles,x(1:node),y(1:node),zs3,oceanfile)
end program extrude
!############################################# countinib01
subroutine countinib01(zlfile,inib)
implicit none
integer(4)    :: i,inib ! inib stands for # of initial boundary nodes
character(70) :: zlfile

open(1,file=zlfile)
inib=0
do while ( i .ge. 0)
read(1,*,end=99)
inib=inib+1
end do
99 continue
close(1)
write(*,*) "# of initial boundary nodes 'inib' is", inib
write(*,*) "### COUNTINIB01 END!! ###"

return
end subroutine countinib01
!############################################# readzlabel02
subroutine readzlabel02(zlfile,inib,zlabel,ibelflg)
implicit none
integer(4)    :: i,inib,zlabel(inib),ibelflg(inib)
character(70) :: zlfile
open(1,file=zlfile)
read(1,*) (zlabel(i),ibelflg(i),i=1,inib)
close(1)
!write(*,'(2i5)')(zlabel(i),ibelflg(i),i=1,100)
write(*,*) "### READZLABEL02 END!! ###"
return
end subroutine readzlabel02
!############################################ calztopo2
subroutine calztopo2(gebcofile,h_mesh,zlabel,ibelflg,inib,coastlinedepth,g_meshpara)
use constants
use param_mesh
use mesh_type
implicit none
type(mesh),            intent(inout)  :: h_mesh
type(meshpara),        intent(in)     :: g_meshpara
integer(4),            intent(in)     :: inib, zlabel(inib), ibelflg(inib)
character(70),         intent(in)     :: gebcofile
real(8),               intent(in)     :: coastlinedepth
real(8),   allocatable,dimension(:)   :: x,y,z
integer(4),allocatable,dimension(:)   :: lbelong,zlabelful ! 2018.08.28
integer(4),allocatable,dimension(:,:) :: n2                ! 2018.08.28
integer(4)                            :: lcount, nlin           ! 2018.08.28
integer(4)                            :: ntopo,node,neast,nsouth,i,j,ii,k,neast1,nsouth1
real(8)                               :: lonorigin,latorigin, calz
real(8),   allocatable,dimension(:)   :: lon,lat,zt,x1,y1
real(8),               dimension(4)   :: xbound,ybound ! 2017.06.27

!#[0]## set input
 node      = h_mesh%node
 nlin      = h_mesh%nlin
 allocate(   zlabelful(node)) ! 2018.08.28
 allocate(   lbelong(nlin)  ) ! 2018.08.28
 allocate(   n2(nlin,2))
 allocate(   x(node),y(node),z(node)) ! 2018.08.28
 x         = h_mesh%xyz(1,:)
 y         = h_mesh%xyz(2,:)
 n2        = h_mesh%n2
 lcount    = h_mesh%nlin
 lbelong   = h_mesh%n2flag(:,2)
 lonorigin = g_meshpara%lonorigin
 latorigin = g_meshpara%latorigin
 xbound    = g_meshpara%xbound  ! 2017.06.27
 ybound    = g_meshpara%ybound  ! 2017.06.27

!# [1] ### count lines of gebcofile
 ntopo=0
 open(1,file=gebcofile)
 do while (ntopo .ge. 0)
  read(1,*,end=99)
  ntopo=ntopo+1 ! ntopo is # of lines in gebcofile
 end do
99 continue

!# [2] ### read coordinates in gebcofile
 allocate(lon(ntopo),lat(ntopo),zt(ntopo),x1(ntopo),y1(ntopo))
 rewind(1)
 do i=1,ntopo
  read(1,*) lon(i),lat(i),zt(i)
 end do
 close(1)

!# [3] ### measure # of nodes in horizontal and vertical directions
 do i=1,ntopo
  if ( lon(i+1) .lt. lon(i) ) then
   neast=i; goto 98 ! neast is # of nodes in horizontal direction
  end if
 end do
98 continue  ! nsouth is # of nodes in vertical direction
 nsouth=ntopo/neast
 write(*,*) "nsouth=",nsouth,"neast=",neast,"ntopo=",ntopo

!# [4] ### convert londitude and  latitude to cartecian coordinates about gebco data
 y1(:) = earthrad*(lat(:)-latorigin)*d2r ! [km] vertical cooridinate
 x1(:) = earthrad*dcos(latorigin*d2r)*(lon(:)-lonorigin)*d2r ! [km] horizontal coordinate
 zt(:) = zt(:)/1000.d0 ! [m] -> [km]
 if (.false.) then ! 2018.11.09
 open(1,file="gebcoz.dat")
  write(1,'(3g15.7)') (x1(i),y1(i),zt(i), i=1,ntopo)
 close(1)
 end if ! 2018.11.09
 deallocate(lon,lat)

!# [5] ### calculate zlabel for all nodes, zlabelful
 zlabelful(1:inib)=zlabel(1:inib) ! zlabelful(i) tells if the i-th node is on coastlines (0) or in the ocean (1).
 do i=1,lcount ! lcount is # of line elements of polygonki.msh
  if ( n2(i,2) .gt. inib ) then ! if n2(i,2), which is the end node of i-th line, is larger than inib, the line is added by gmsh -2
  zlabelful(n2(i,2))=ibelflg(lbelong(i)) ! ibelflg tells if the initial line is on the coastlines (0) or not (1).
  end if
 end do
 zlabelful(lcount+1:node)=1 ! When i is larger than lcount, i-th node must be in the ocean (1)

!# [6] ### calculate depth "z" at all the horizontal nodes
 z(:)=0.d0
 do k=1,node ! node is total # of nodes in ocean surface horizontal mesh, polygon_ki.msh
 !if ( zlabelful(k) .eq. 0 ) cycle ! Patern 1 : if the k-th node is on the coastline, we don't have to calculate z
  if ( zlabelful(k) .eq. 0 ) then !  Patern 2:  if the k-th node is on the coastline, we set z=-0.01d0 [km]
   z(k)=coastlinedepth
   cycle
  end if
  neast1=0; nsouth1=0
  do i=1,nsouth-1
   if ( y1((i-1)*neast+1) .ge. y(k) .and. y(k) .ge. y1(i*neast+1) ) nsouth1=i
  end do
  do j=1,neast-1
   if ( x1(j) .le. x(k) .and. x(k) .le. x1(j+1)) neast1=j
  end do
  if ( neast1 .eq. 0 .or. nsouth1 .eq. 0 ) goto 999
  ii=(nsouth1-1)*neast+neast1
  z(k)=calz(x(k),y(k),x1(ii),x1(ii+1),y1(ii),y1(ii+neast),zt(ii),zt(ii+1),zt(ii+neast),zt(ii+neast+1))
  if (z(k)  .gt. 0.d0 ) z(k)=-0.1d0

!######################################################################## 2017.06.27 modified
if ( (x(k) .lt. xbound(2) .or. xbound(3) .lt. x(k) ) .or. &
&    (y(k) .lt. ybound(2) .or. ybound(3) .lt. y(k) )) then ! 2017.06.27
 if ( z(k) .gt. - 0.1d0 ) z(k) = -0.1d0 ! 2017.06.27 to avoid too small dihedral angles
else         ! 2017.06.27
 if ( z(k) .gt. -0.015d0 ) z(k) = -0.015d0 ! inside the focus area
end if
!######################################################################## 2017.06.27 modified

!write(*,*) "x1(ii),x(ii+1),y(ii),y(ii+neast),zt(ii),zt(ii+1),zt(ii+neast),zt(ii+neast+1)"
!write(*,*) x1(ii),x1(ii+1),y1(ii),y1(ii+neast),zt(ii),zt(ii+1),zt(ii+neast),zt(ii+neast+1)
!write(*,*) "nsouth1=",nsouth1,"neast1=",neast1
!write(*,*) "k=",k,"x,y=",x(k),y(k),"z(k)=",z(k)
!stop
end do
! output topo file
if (.false.) then ! 2018.11.09
open(1,file="polygonz.dat")
write(1,'(3g15.7)') (x(i),y(i),z(i), i=1,node)
close(1)
end if ! 2018.11.09
!# set output
h_mesh%xyz(3,:) = z

write(*,*) "### CALZTOPO2 END!! ###"
return
999 continue
write(*,*) "GEGEGE when node # is k=",k,"x(k),y(k)=",x(k),y(k)
write(*,*) "neast1=",neast1,"nsouth1=",nsouth1
stop
end subroutine calztopo2
!######################################################  function calz
function calz(x,y,x1,x2,y1,y2,z1,z2,z3,z4)
implicit none
real(8) :: calz,x,y,x1,x2,y1,y2,z1,z2,z3,z4,aa,a1,a2,a3,a4
aa=(x2-x1)*(y1-y2)
a1=(x2-x)*(y-y2)/aa
a2=(x-x1)*(y-y2)/aa
a3=(x2-x)*(y1-y)/aa
a4=(x-x1)*(y1-y)/aa
calz=a1*z1+a2*z2+a3*z3+a4*z4
end function calz
!##################################################### polygonz3
subroutine polygonz3(mshfile,polyzfile,z,nmax) ! output polyzfile (mshfile)
implicit none
real(8),  dimension(nmax) :: z
integer(4)                :: i,j,node,nmax
real(8)                   :: x1,y1,z1
character(70)             ::  mshfile,polyzfile
character(80)             :: line

open(1,file=mshfile)
open(2,file=polyzfile)
do i=1,4
read(1,'(a80)') line
write(2,'(a80)') line
end do
read(1,*) node
write(2,*) node
do i=1,node
read(1,*) j,x1,y1,z1
write(2,*) j,x1,y1,z(i)
end do
do while ( i .ge. 0)
read(1,'(a80)',end=99) line
write(2,'(a80)') line
end do
99 continue
close(1)
close(2)

write(*,*) "### POLYGONZ3 END!! ###"
return
end subroutine polygonz3
!#################################################### countznode4
!# modified on 2017.06.27
!# modified on 2016.09.13
subroutine  countznode4(ns3,nz,h_mesh,g_meshpara,node) ! 2018.08.28
use param_mesh
use mesh_type
use horizontalresolution ! for value
implicit none
type(mesh),            intent(in)     :: h_mesh
type(meshpara),        intent(in)     :: g_meshpara
integer(4),            intent(in)     :: node
integer(4),            intent(out)    :: nz(node),ns3 ! # of vertical nodes for each horizontal column
real(8),   allocatable,dimension(:,:) :: xy
real(8),   allocatable,dimension(:)   :: z, depth_nlayer
integer(4)                            :: i,j,iflag,nlayer_max

!#[1]##
 nlayer_max   = g_meshpara%nlayer_max
 allocate(xy(2,node),z(node))
 allocate(depth_nlayer(nlayer_max-1))
 xy(1:2,:)    = h_mesh%xyz(1:2,:)
 z(:)         = h_mesh%xyz(3,:)
 iflag        = g_meshpara%iflag_extrude ! 0 for nothing, 1 : vertical interval is less than horizontal
 depth_nlayer = g_meshpara%depth_nlayer
 write(*,*) "g_meshpara%iflag_extrude=",iflag

!#[2]##
ns3   = 0 ! number of nodes in the ocean
nz(:) = 0
do i=1,node
 if (z(i) .eq. 0.d0 ) nz(i) =1 ! only coastline nodes
 if ( z(i) .lt. 0.d0 ) then
!  nz(i)= 2 + int((abs(z(i)) - 0.001 )/dz_base) ! z(i) is negative, namely upward positive
! if ( -2.d0 .lt. z(i) .and. z(i) .lt. -0.05d0 ) nz(i)=max(3,nz(i))
! !#[2-2]## vertical resolution is reduced so less than horizontal, inserted on 2016.09.13
  nz(i) = 2
  do j = 1,nlayer_max -1
   if ( z(i) .lt. depth_nlayer(j) ) nz(i) = j + 2
  end do
 end if
 if ( iflag .eq. 1 .and. nz(i) .gt. 1 ) then
!  dz  = abs(z(i))/float(nz(i)-1)
!  dz2 = value(xy(1,i),xy(2,i),g_meshpara)
!  if ( dz2 .lt. dz ) then
!   !write(*,*) "dz=",dz,"dz2=",dz2,"nz(i)=",nz(i),"->", int(abs(z(i)) / dz2 + 1)
!   nz(i) = abs(z(i)) / dz2 + 1
!  end if

 else if ( iflag .ge. 2 ) then
  goto 99
 end if

 ns3=ns3+nz(i)
 !if (z(i) .le. 100 ) ns=ns+2 ! not insert nodes, a node on the surface and bottom
 !if (100 .le. z(i) .and. z(i) .le. 1000) ns=ns+3 ! insert one node
end do
write(*,*) "Total number of node in the ocean (ns3)=",ns3
write(*,*) "### COUNTNODE4 END!! ###"
return
!
99 continue
write(*,*) "GEGEGE! g_meshpara%iflag_extrude",iflag ,">= 2 is not supported."
stop
end subroutine countznode4
!################################################### mkznode5
subroutine mkznode5(zs3,z,ns3,node,nz,nhnptr,nmax)
implicit none
integer(4), intent(in)    :: nmax,node,ns3 ! 2018.08.28
real(8),    intent(in)    :: z(nmax)
integer(4), intent(in)    :: nz(nmax)      ! # of vertical nodes for each horizontal column
real(8),    intent(inout) :: zs3(ns3)
integer(4), intent(inout) :: nhnptr(ns3)   ! indicates node columun # which the node in the ocean belongs to
integer(4)                :: i,j,k,ii,jj

ii=0
do i=1,node
ii=ii+1
nhnptr(ii)=i
zs3(ii)=0.d0
if (z(i) .eq. 0.d0) cycle
jj=nz(i)-2 ! jj=z(i)/1000.
do k=1,jj+1 ! from 1 to nz - 1; second to bottom node
ii=ii+1
nhnptr(ii)=i
zs3(ii)=z(i)*float(k)/float(jj+1) ! uniform vertical interval is adopted
end do
end do
if ( ns3 .ne. ii) goto 99
write(*,*) "### MKZNODE5 END!! ###"
return

99 write(*,*) "GEGEGE ii is not equal to ns3; ii and ns3=",ii,ns3
stop
end subroutine mkznode5
!###################################################### polygonallnode6
subroutine polygonallnode6(allnodefile,x,y,node,zs3,ns3,nz,nmax)
! output allnodefile (geofile)
implicit none
integer(4),   intent(in) :: node,ns3,nmax
real(8),      intent(in) :: x(nmax),y(nmax)
integer(4),   intent(in) :: nz(nmax) ! # of vertical nodes for each horizontal column
real(8),      intent(in) :: zs3(ns3)
character(70),intent(in) :: allnodefile
integer(4)               :: i, j, ii

open(1,file=allnodefile)
ii=0
!write(1,'(a11)') "$MeshFormat"
!write(1,'(a7)') "2.2 0 8"
!write(1,'(a14)') "$EndMeshFormat"
!write(1,'(a6)') "$Nodes"
!write(1,*) ns3
do i=1,node
do j=1,nz(i)
ii=ii+1
!write(1,*) ii,x(i),y(i),zs3(ii)
write(1,'(a)') "lc=500.0;"
write(1,'(a,i10,a,3(g15.7,a))') "Point(",ii,")={",x(i),",",y(i),",",zs3(ii),",lc};"
end do
end do
!write(1,'(a9)')"$EndNodes"
close(1)
write(*,*) "### POLYGONALLNODEZ6 END!! ###"

return
end subroutine polygonallnode6
!###################################################### sortzorder7
subroutine sortzorder7(zs3,ns3,nzptr)
implicit none
integer(4) :: i,j,u1,ns3
real(8),dimension(ns3) :: zs3
integer(4),dimension(ns3) :: nzptr
nzptr(:)=0
nzptr(1)=1 ! this program assumes that z of first node = 0.d0, one of shallowest nodes
do i=2,ns3 ! ns3 is the # of total node for 3-D ocean mesh
! seek where is the suit place for i-th node
do j=1,i-1
  if (zs3(nzptr(j)) .ge. zs3(i) )then
    u1=j       ! j is the nearest upper node, and i-th node should be u1+1-th shallow point
  else
    goto 99
  end if
end do
99 continue
if ( u1 .eq. i -1) nzptr(i) = i
if ( u1 .lt. i  ) then
nzptr(u1+2:i)=nzptr(u1+1:i-1)
nzptr(u1+1)=i
end if
end do
!write(*,*) (zs3(nzptr(i)),i=1,ns3)
write(*,*) "### SORTZORDER7 END!! ###"
return
end subroutine sortzorder7
!####################################################   countele8
subroutine countele8(n3,ecount,node,neleptr,nelecount,neles,ns3,nz) ! neleptr(i:node,1:10) indicates which element column i-th horizontal node is included
implicit none ! node is # of horizontal oceanic surface node, ns3 is total # of 3D oceanic nodes
integer(4) ::i,j,ii,ecount,node,neles,ns3 ! ecount is # of triangle elements; neles is total # of 3D ocean elements
integer(4),dimension(node) :: nelecount ! nelecount(i) is # of elements which i-th horizontal node belongs to
integer(4),dimension(node,10) :: neleptr ! neleptr(i) is element # which i-th horizontal node belongs to
integer(4),dimension(ecount,3) :: n3
integer(4),dimension(node)  :: nz
nelecount(:)=0
neleptr(:,:)=0
neles=0
do i=1,ecount ! make nelecount(1:node) node is # of horizontal surface nodes
do j=1,3
!write(*,*) "i=",i,"nelecount(13260)=",nelecount(13260),"ii=n3(i,j)=",n3(i,j)
ii=n3(i,j)
nelecount(ii)=nelecount(ii)+1
!if ( ii .eq. 13260 .or. i .eq. 25907 )then!
!write(*,*) "i,j=",i,j
!write(*,*) "n3(i,j)=",n3(i,j)
!write(*,*) "new nelecount(n3(i,j))=",nelecount(ii)
!end if
!write(*,*) "i,j=",i,j,"nelecount(n3(i,j))=",nelecount(n3(i,j))
neleptr(ii,nelecount(ii))=i ! set the element # to neleptr
neles=neles+(nz(ii)-1)
end do
end do
write(*,*) "Total # of tetrahedron elements (neles) =",neles
write(*,*) "### COUNTELE8 END!! ###"
return
end subroutine countele8
!########################################################################   maketetra9
subroutine maketetra9(n3,ecount,node,nzptr,neleptr,nhnptr,nelecount,ns3,nz,neles,x,y,zs3,oceanfile)
implicit none
integer(4) :: i,j,k,l,ii,ecount,node,neles,ns3,newnode,unode,hnode,nn,hnode2,hnode3
! ns3 is total # of tetrahedrons
integer(4),dimension(ecount,3) :: n3 ! ecount is # of triangle elements, n3(i,j) is horizontal node number which is j-th node of i-th element
integer(4),dimension(ns3)      :: nzptr,nhnptr
integer(4),dimension(node,10)  :: neleptr ! node is # of horizontal node, neleptr(i,j) is horizontal elemtnt # which is j-th element of i-th horizontal node
integer(4),dimension(node)     :: nelecount,nz,hnns3ptr,hanging_node ! hanging_node(i) is the accumulated # of nodes at i-th node column
integer(4),dimension(4)        :: n4
real(8),   dimension(node)     :: x,y
real(8),   dimension(ns3)      :: zs3 ! zs3(i) is the z coordinate of i-th 3-D oceanic node
real(8),   dimension(4)        :: x1,y1,z1
real(8) :: tetvol_old,vol
character(70) :: oceanfile
!# [0] ### make hnns3ptr : hnns3ptr(i) is the ocean node # which corresponds to the i-th horizontal surface node
hnns3ptr(:)=0;hnns3ptr(1)=1
do i=2,node
hnns3ptr(i)=hnns3ptr(i-1)+nz(i-1)
!write(*,*) "hnns3ptr(i),nz(i)=",hnns3ptr(i),nz(i)
end do
! output n3(1:ecount,1:3)
open(1,file="n3ocean.dat")
do i=1,ecount
write(1,*) i,(n3(i,j),j=1,3)
end do
close(1)
open(1,file=oceanfile)
! out oceanfile ### node
write(1,'(a11)') "$MeshFormat"
write(1,'(a7)') "2.2 0 8"
write(1,'(a14)') "$EndMeshFormat"
write(1,'(a6)') "$Nodes"
write(1,*) ns3
ii=0
do i=1,node
!write(*,*) "Node i=",i
do j=1,nz(i)
ii=ii+1
write(1,'(i10,3g15.7)') ii,x(i),y(i),zs3(ii) ! 2018.08.28
end do
end do
write(1,'(a9)') "$EndNodes"
! out oceanfile ### tetrahedron elements
write(1,'(a9)') "$Elements"
write(1,'(i10)') neles+node
ii=0
do i=1,node ! ocean surface node # are recorded
ii=ii+1
write(1,'(i10,a,i10)') ii,"  15   2   0   1 ",hnns3ptr(i) ! Point elements
end do
hanging_node(:)=1 ! At all node column, the starting node is on the surface
do i=node+1,ns3
!write(*,*) "Element i=",i
!# [1] choose first two nodes for i-th tetrahedrons element ###
newnode=nzptr(i) ! in this loop new node is added and create tetrahedron
unode=newnode-1 ! node # of just above the newnode
!write(*,*) "newnode,unode=",newnode,unode
if ( zs3(unode) .le. zs3(newnode)) goto 99
hnode=nhnptr(newnode) ! nhnptr(i) is horizontal node # corresponds to i-th 3-D node
!write(*,*) "hnode=",hnode,"nelecount(hnode)=",nelecount(hnode)
do j=1,nelecount(hnode) ! nelecount(i) is # of element column surrounding i-th horizontal node
! loop for element columns which surrounding the nhnptr(i) node column
!# [2] choose rest two nodes for the tetrahedron ###
!write(*,*) "j=",j,"neleptr(hnode,j)=",neleptr(hnode,j)
hnode2=0;hnode3=0 ! hnode means horizontal surface node
do k=1,3
nn=n3(neleptr(hnode,j),k)
! n3(i,j) is horizontal node # which is j-th node of i-th element
! neleptr(i,j) is horizontal element # which is j-th elements of i-th horizontal node
if (hnode .ne. nn ) hnode2=nn ! if hnode2 is not given, set hnode2
end do
do k=1,3
nn=n3(neleptr(hnode,j),k)
if (hnode .ne. nn .and. hnode2 .ne. nn) hnode3=nn ! if hnode2 is already given, set hnode3
end do
if ( hnode2 .eq. hnode3 )then
write(*,*) "GEGEGE hnode2=",hnode2,"hnode3=",hnode3
stop
end if
!write(*,*) "hnode2=",hnode2,"hnode3=",hnode3
!# [3] ### set four node #
n4(1:4)=(/unode,hanging_node(hnode2)+hnns3ptr(hnode2)-1,&
&                          hanging_node(hnode3)+hnns3ptr(hnode3)-1,&
&               newnode/)
!if ( n4(2) .eq. 112826 .or. n4(3) .eq. 112826) then
x1(1)=x(hnode);x1(2)=x(hnode2);x1(3)=x(hnode3);x1(4)=x(hnode)
y1(1)=y(hnode);y1(2)=y(hnode2);y1(3)=y(hnode3);y1(4)=y(hnode)
z1(1)=zs3(n4(1));z1(2)=zs3(n4(2));z1(3)=zs3(n4(3));z1(4)=zs3(n4(4))
vol=tetvol_old(x1,y1,z1)
if ( vol < 0.d0) then
l=n4(2)
n4(2)=n4(3)
n4(3)=l
!write(*,*) "element node order changed! element # is", ii+1
end if
ii=ii+1
write(1,'(i10,a,4i10)') ii,"  4   2   0   1",(n4(l),l=1,4)
end do ! create element loop end
hanging_node(hnode)=hanging_node(hnode)+1 ! add 1 to the accumulated number of hanging_node at node columun nhnptr(i)
!if ( hnode .eq. 38534) then
!write(*,*) "i=",i,"newnode=nzptr(i)=",newnode,"hnode=",hnode
!write(*,*) hanging_node(hnode)
!write(*,*) "zs3(newnode)=",zs3(newnode)
!end if
end do ! node loop end
if ( ii .ne. neles + node) goto 98
write(1,'(a12)') "$EndElements"
close(1)
write(*,*) "### MAKETETRA9 END!! ###"
return
99 continue
write(*,*) "GEGEGE i=",i,"newnode=",newnode,"unode=",unode
write(*,*) "z_unode=",zs3(unode)," is deeper than z_newnode=",zs3(newnode)
stop
98 continue
write(*,*) "GEGEGE Total number of oceanic elements (neles) =",neles,"is not equal to created elements (ii)=",ii
end subroutine maketetra9
!#####################################  function tetvol_old
!This program is from Field_base.f90 in Fluidity
function tetvol_old( x, y, z )
implicit none
real(8) tetvol_old
real(8) x(4), y(4), z(4)
real(8) vol, x12, x13, x14, y12, y13, y14, z12, z13, z14
!
x12 = x(2) - x(1)
x13 = x(3) - x(1)
x14 = x(4) - x(1)
y12 = y(2) - y(1)
y13 = y(3) - y(1)
y14 = y(4) - y(1)
z12 = z(2) - z(1)
z13 = z(3) - z(1)
z14 = z(4) - z(1)
!
vol = x12*( y13*z14 - y14*z13 )  &
+ x13*( y14*z12 - y12*z14 ) &
+ x14*( y12*z13 - y13*z12 )
!
tetvol_old = vol/6
!
return
end function tetvol_old
