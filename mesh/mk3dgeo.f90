! line 29 is added on May 30, 2016
! to decrease jumps in flatness of elements
!================================
!<node configuration>
! polygonki.msh                                     ->   ocean.msh
! [ ocean node by polygonki.geo (nodeseageo)    ]        [ ocean node by polygonki.geo      ]
! [ corner node for land        (naddedcorner)  ]        [ added node on oceanbound         ]
! [ added node on oceanbound in polygonki.msh   ]        [ added node in below ocean surface]
! [ added node on landbound in polygonki.msh    ]
! [ node in ocean in polygonki.msh              ]
! [ node in land in polygonki.msh               ]
!==
program mk3dgeo
use constants  ! added on 2016.09.06
use mesh_type  ! added on 2016.09.07
use param_mesh ! see m_param_mesh.f90, added on Sep. 7, 2016
use mesh_relation
use line_type
use topo_tool      ! 2018.08.27
use coastline_data ! 2018.08.27
implicit none
type(mesh)      :: ocean_mesh, ki_mesh
type(meshpara)  :: g_meshpara
type(relate)    :: g_relate
type(line_info) :: l_line
integer(4)      :: n2ds ! ocean.msh,n3ds is 3D nodes in the ocean including boundary nodes
integer(4)      :: ntri ! 2018.11.18
integer(4)      :: nodek,ntrik,ntriall ! for polygonki.msh
real(8),   allocatable,dimension(:,:) :: xyz,xyz3d
integer(4),allocatable,dimension(:,:) :: n3all
integer(4),allocatable,dimension(:)   :: nzs,n3top,n3bot
integer(4),allocatable,dimension(:,:) :: ki23dptr ! ki23dptr(i,j) is i=1, surface, i=2 bottom
integer(4)      :: n3d, n3dn ! n3d: total number of 3d geo file, n3dn= n3d + 8 (= 3d corner)
integer(4)      :: nlinbry(4),linbry_sb(2,4,100) ! nlinbry(i) is # of lines on i-th calculation boundary, linbry(i,j) is line number of j-th line of i-th boundary
integer(4)      :: m4(4,2)
character(70)   :: oceanfile="ocean.msh",head,kifile ! input files
character(70)   :: premsh="pre3d.msh" ! subroutine premsh6_5
character(70)   :: pregeo="pre3d.geo" ! subroutine outpregeo8
character(70)   :: bgmeshfile="bgmesh3d.pos"
character(70)   :: ctlfile="mshki2ocean.ctl",ki23dfile="ki23dptr.dat"
!#
character(70)   :: geofile,geofileki,gebcofile ! 2018.08.27
integer(4)      :: node0                       ! 2018.08.27
type(grid_data) :: gebco_grd                   ! 2018.08.27 see m_coastline_data.f90

!#[0]## read meshparameters for head
 CALL READMESHPARA(g_meshpara)
 CALL calcartesianbound(g_meshpara)
     head=g_meshpara%head
     kifile=head(1:len_trim(head))//"ki.msh"
 CALL READGEORELATION(g_relate,ctlfile) ! read nodeseageo, nodegeo, naddcorner

!#[1]## read mesh
 call READMESH_TOTAL(ocean_mesh, oceanfile)
 call READMESH_TOTAL(ki_mesh, kifile)
  n2ds  = ocean_mesh%npoi  ! n2ds is # of 2d ocean nodes
  nodek = ki_mesh%node
  ntrik = ki_mesh%ntri

!#[2]## read topography and reflect topography in ki_mesh
 CALL calzkiland(g_meshpara,g_relate,ki_mesh)  ! 2018.08.27

!#[3]## make n3d
 call findnodenum3(ki_mesh,ocean_mesh,g_relate,n3d)

!#[4]## make nzs
     allocate (nzs(n2ds))
 call makenzs5(nzs,ocean_mesh) ! nzs(i) is # of ocean nodes in the i-th node column

!#[5]## make xyz for pre3d.geo
     allocate( ki23dptr(2,nodek),xyz(3,n3d))
 call makenode4(ki23dptr,g_relate,xyz,n3d,ki_mesh,ocean_mesh,nzs)
 call outki23dptr(nodek,ki23dptr,ki23dfile)

!#[6]## make n3top and n3bot, which store the node # in n3d nodes
     allocate(n3top(ntrik),n3bot(ntrik),n3all(ntrik*2,3))
 call maketopbot6(n3top,n3bot,n3all,ntriall,ki23dptr,ki_mesh)
 call premsh6_5(xyz,n3d,n3all(1:ntriall,1:3),ntriall,premsh) ! output pre3d.msh

!#[7]## find horizontal boundary lines
 call MKLINE_V2(l_line, n3d, ntriall,3,n3all(1:ntriall,1:3)) ! make l_line%line(nline,2) see m_line_type.f90
 call MKN3  (l_line, n3d, ntriall, n3all(1:ntriall,1:3)) ! make l_line%n3line(ntriall,3)
     deallocate(n3all)
 call findbryline7(ki_mesh,nlinbry,linbry_sb,l_line,ki23dptr,xyz,n3d,m4)

!#[8]## make last node set
       n3dn=n3d+8
     allocate(xyz3d(3,n3dn)) ! n3dn is total node # of output pre3d.geo
 call prepare7_5(xyz,n3d,xyz3d,n3dn,g_meshpara)

!#[9]## output pregeo file
 ntri = l_line%ntri ! 2018.11.18
 call outpregeo8(l_line,n3top,n3bot,xyz3d,n3d,pregeo,nlinbry,&
 &               linbry_sb,ki_mesh,g_meshpara,m4,ntri) ! 2018.11.18
 call outbgmesh3d(bgmeshfile,g_meshpara) ! see outbgmesh3d.f90

end program mk3dgeo
!########################################## calzkailand
!# 2018.08.28
!# copy the part of calztopo2 in extrude.f90
subroutine calzkiland(g_meshpara,g_relate,ki_mesh)
use param_mesh
use mesh_type
use constants
use mesh_relation
implicit none
type(meshpara), intent(in)         :: g_meshpara
type(mesh),     intent(inout)      :: ki_mesh
type(relate),   intent(in)         :: g_relate
integer(4)                         :: ntopo
real(8),allocatable,dimension(:)   :: lon,lat,zt,x1,y1
real(8),allocatable,dimension(:,:) :: xyz
real(8)                            :: lonorigin,latorigin
integer(4)                         :: i, j, ii
integer(4)                         :: nsouth,neast,nsouth1,neast1
integer(4)                         :: nodek, maxseabry
integer(4)                         :: nodegeo,nodeseageo,maxlandbry,maxinsea
character(50)                      :: gebcofile

!# [0] ### set
 gebcofile  = g_meshpara%topofile
 lonorigin  = g_meshpara%lonorigin
 latorigin  = g_meshpara%latorigin
 nodegeo    = g_relate%nodegeo
 nodeseageo = g_relate%nodeseageo
 maxseabry  = g_relate%maxseabry
 maxlandbry = g_relate%maxlandbry
 maxinsea   = g_relate%maxinsea
 nodek      = ki_mesh%node
 xyz        = ki_mesh%xyz

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

!# [4] ### convert lon and lat to cartecian coordinates about gebco data
 y1(:) = earthrad*(lat(:)-latorigin)*d2r ! [km] vertical cooridinate
 x1(:) = earthrad*dcos(latorigin*d2r)*(lon(:)-lonorigin)*d2r ! [km] horizontal coordinate
 zt(:) = zt(:)/1000.d0 ! [m] -> [km]
 deallocate(lon,lat)

!# [5] ## cal z from topo
 do i=nodeseageo+1,nodegeo ! land boundary geo nodes [LAND]
  call findcalz(ntopo,nsouth,neast,x1,y1,zt,xyz(1:2,i),xyz(3,i))
 end do
 do i=maxseabry+1,maxlandbry ! land boudary nodes added by gmsh [LAND]
  call findcalz(ntopo,nsouth,neast,x1,y1,zt,xyz(1:2,i),xyz(3,i))
 end do
 do i=maxinsea+1,nodek ! nodes in land added by gmsh [LAND]
  call findcalz(ntopo,nsouth,neast,x1,y1,zt,xyz(1:2,i),xyz(3,i))
 end do

!# [6] ## output
 ki_mesh%xyz = xyz

 write(*,*) "### CALZKILAND END!! ###"
 return
end
!####################################################### find and calz
!# 2018.08.28
subroutine findcalz(ntopo,nsouth,neast,x1,y1,zt,xy,zout)
implicit none
integer(4),intent(in)  :: ntopo,nsouth,neast
real(8),   intent(out) :: zout
real(8),   intent(in)  :: xy(2)
real(8),   intent(in)  :: x1(ntopo),y1(ntopo),zt(ntopo)
real(8)                :: calz
integer(4)             :: i, j, neast1, nsouth1, ii

  do i=1,nsouth-1
   if ( y1((i-1)*neast+1) .ge. xy(2) .and. xy(2) .ge. y1(i*neast+1) ) nsouth1=i
  end do
  do j=1,neast-1
   if ( x1(j) .le. xy(1) .and. xy(1) .le. x1(j+1)) neast1=j
  end do
  if ( neast1 .eq. 0 .or. nsouth1 .eq. 0 ) goto 999
  ii=(nsouth1-1)*neast+neast1
  zout=calz(xy(1),xy(2),x1(ii),x1(ii+1),y1(ii),y1(ii+neast),&
&           zt(ii),zt(ii+1),zt(ii+neast),zt(ii+neast+1))

  return

999 continue
write(*,*) "GEGEGE! No corresponding coordinate; xy(2)=",xy(1),xy(2)
write(*,*) "neast1=",neast1,"nsouth1=",nsouth1
stop
end
!######################################################  function calz
!# 2018.08.28
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
!########################################## out ki23dfile
subroutine outki23dptr(nodek,ki23dptr,ki23dfile)
implicit none
integer(4),intent(in) :: nodek,ki23dptr(2,nodek)
character(70),intent(in) :: ki23dfile
integer(4) :: i

open(1,file=ki23dfile)
 write(1,*) nodek
 do i=1,nodek
  write(1,*) ki23dptr(1:2,i)
 end do
close(1)

return
end
!########################################## findnodenum3
subroutine findnodenum3(ki_mesh,ocean_mesh,g_relate,n3d)
use mesh_type
use mesh_relation
implicit none
type(mesh),intent(in) :: ki_mesh,ocean_mesh
type(relate),intent(inout) :: g_relate
integer(4),intent(out) :: n3d
integer(4) :: nodegeo,nodeseageo,maxseabry,maxlandbry,maxinsea
integer(4) :: i,nlineinsea,neleinsea
integer(4) :: nodek,ntrik,n2ds,n3ds,nlink
integer(4) :: n3k(ki_mesh%ntri,4),n2k(ki_mesh%nlin,3)
!#[1]## set input
n2ds       = ocean_mesh%npoi
n3ds       = ocean_mesh%node
nodek      = ki_mesh%node
nlink      = ki_mesh%nlin
ntrik      = ki_mesh%ntri
n2k(:,1  ) = ki_mesh%n2flag(:,2)
n2k(:,2:3) = ki_mesh%n2(:,1:2)
n3k(:,1  ) = ki_mesh%n3flag(:,2)
n3k(:,2:4) = ki_mesh%n3(:,1:3)
nodegeo    = g_relate%nodegeo
nodeseageo = g_relate%nodeseageo

maxlandbry = 0
maxseabry  = 0
nlineinsea = 0

do i=1,nlink
 if (n2k(i,1) .le. nodeseageo ) then
  maxseabry=max(maxseabry,n2k(i,2),n2k(i,3))
  nlineinsea=nlineinsea+1
 else if (n2k(i,1) .gt. nodeseageo ) then
  !write(*,*) "maxlandbry, n2k(i,2:3)", maxlandbry,n2k(i,2),n2k(i,3)
  !write(*,*) "nodeseageo=",nodeseageo,"n2k(i,1)=",n2k(i,1)
  maxlandbry=max(maxlandbry,n2k(i,2),n2k(i,3))
 end if
end do

if ( maxlandbry .lt. maxseabry ) maxlandbry=maxseabry
! if no nodes are added by gmsh -2 for land baoundaries
write(*,*) "nodeseageo=", nodeseageo
write(*,*) "nodegeo=",    nodegeo
write(*,*) "maxseabry=",  maxseabry
write(*,*) "maxlandbry=", maxlandbry
write(*,*) "nlineinsea=", nlineinsea
! find maxinsea
maxinsea=0
neleinsea=0
do i=1,ntrik
 if (n3k(i,1) .eq. 1 ) then ! n3k(i,1) = 1 indicates the element is in the ocean
  maxinsea=max(maxinsea,n3k(i,2),n3k(i,3),n3k(i,4))
  neleinsea=neleinsea+1
 end if
end do
write(*,*)"maxinsea=",maxinsea
write(*,*)"neleinsea=",neleinsea

! calculate n3d
n3d=n3ds+( nodegeo - nodeseageo ) + ( maxlandbry - maxseabry )+ &
& ( nodek -maxinsea)
write(*,*) "# of total land ocean nodes (n3d) is",n3d

! calculate neles
!neles=nodeseageo+nlineinsea+neleinsea
!write(*,*) "neles=",neles
if (n2ds .ne. nodeseageo + (maxseabry-nodegeo)+(maxinsea-maxlandbry) ) then
write(*,*) "GEGEGE n2ds=",n2ds,"is not equal to"
write(*,*) "nodeseageo + (maxseabry-nodegeo)+(maxinsea-maxlandbry)=",&
&                nodeseageo + (maxseabry-nodegeo)+(maxinsea-maxlandbry)
stop
end if

!# set output
g_relate%maxlandbry = maxlandbry
g_relate%maxseabry = maxseabry
g_relate%maxinsea = maxinsea

write(*,*) "### FINDNODENUM3 END!! ###"
return
end subroutine findnodenum3
!##################################################  makenode4
subroutine makenode4(ki23dptr,g_relate,xyz,n3d,ki_mesh,ocean_mesh,nzs)
use mesh_relation
use mesh_type
implicit none
type(relate),intent(in) :: g_relate
type(mesh),intent(in)  :: ki_mesh, ocean_mesh
integer(4),intent(in)  :: n3d,nzs(ocean_mesh%npoi)
real(8)   ,intent(out) :: xyz(3,n3d)
integer(4),intent(out) :: ki23dptr(2,ki_mesh%node)
integer(4) :: nodegeo,nodeseageo,maxseabry,maxlandbry,maxinsea
integer(4) :: nodek,n2ds,n3ds,n1s(ocean_mesh%npoi,1)
integer(4) :: i,ii,j

!#[1]## set input
 nodegeo    = g_relate%nodegeo
 nodeseageo = g_relate%nodeseageo
 maxseabry  = g_relate%maxseabry
 maxlandbry = g_relate%maxlandbry
 maxinsea   = g_relate%maxinsea
 nodek      = ki_mesh%node
 n3ds       = ocean_mesh%node
 n2ds       = ocean_mesh%npoi
 n1s        = ocean_mesh%n1


!#[2]## set ocean nodes to xyz
xyz(1:3,1:n3ds) = ocean_mesh%xyz(1:3,1:n3ds) ! ocean.msh nodes are taken after by pre3d.geo

!write(*,*)"n3ds=",n3ds
ii=0 ! ii is node # in n2ds
j=0 ! j is extra node # after n3ds
do i=1,nodeseageo ! ocean boundary geo nodes [OCEAN]
 ii=ii+1
 ki23dptr(1,i)=n1s(ii,1)             ! surface
 ki23dptr(2,i)=n1s(ii,1) + nzs(ii)-1 ! bottom
end do

!write(*,*) "2 x,y,z(1)=",x(1),y(1),z(1)
do i=nodeseageo+1,nodegeo ! land boundary geo nodes [LAND]
 j=j+1 ! j + n3ds is the node number of pre3d.geo nodes
 ki23dptr(1,i)  = n3ds+j                ! surface
 ki23dptr(2,i)  = n3ds+j                ! bottom
 xyz(1:3,n3ds+j)= ki_mesh%xyz(1:3,i)
end do
!write(*,*) "3 x,y,z(1)=",x(1),y(1),z(1)
do i=nodegeo+1,maxseabry ! ocean boundary nodes added by gmsh [OCEAN]
 ii=ii+1
 ki23dptr(1,i)=n1s(ii,1)               ! surface
 ki23dptr(2,i)=n1s(ii,1) + nzs(ii) - 1 ! bottom
end do
!write(*,*) "4 x,y,z(1)=",x(1),y(1),z(1)
do i=maxseabry+1,maxlandbry ! land boudary nodes added by gmsh [LAND]
 j=j+1
 ki23dptr(1,i)=n3ds+j                 ! surface
 ki23dptr(2,i)=n3ds+j                 ! bottom
 xyz(1:3,n3ds+j)=ki_mesh%xyz(1:3,i)
end do
!write(*,*) "5 x,y,z(1)=",x(1),y(1),z(1)
do i=maxlandbry+1,maxinsea ! nodes in ocean added by gmsh [OCEAN]
ii=ii+1
ki23dptr(1,i)=n1s(ii,1)               ! surface
ki23dptr(2,i)=n1s(ii,1) + nzs(ii) - 1 ! bottom
end do
!write(*,*) "6 x,y,z(1)=",x(1),y(1),z(1)
write(*,*) "nodeseageo=",nodeseageo
write(*,*) "nodegeo=",nodegeo
write(*,*) "maxseabry=",maxseabry
write(*,*) "maxlandbry=",maxlandbry
write(*,*) "maxinsea=",maxinsea
write(*,*) "nodek=",nodek
write(*,*) "n3ds=",n3ds
write(*,*) "n3d=",n3d
do i=maxinsea+1,nodek ! nodes in land added by gmsh [LAND]
 j=j+1
 ki23dptr(1,i)=n3ds+j   ! surface
 ki23dptr(2,i)=n3ds+j   ! bottom
 xyz(1:3,n3ds+j)=ki_mesh%xyz(1:3,i)
 !write(*,*) "nodek,i=",nodek,i,"j=",j,"n3d,n3ds+j=",n3d,n3ds+j
 !write(*,*) " x,y,z(1)=",x(1),y(1),z(1)
end do

!write(*,*) "7 x,y,z(1)=",x(1),y(1),z(1)
return
end subroutine makenode4
!########################################## make nzs5
subroutine makenzs5(nzs,ocean_mesh)
use mesh_type
implicit none
type(mesh),intent(in) :: ocean_mesh
integer(4),intent(out) :: nzs(ocean_mesh%npoi)
integer(4) :: n2ds,n3ds,n1s(ocean_mesh%npoi,1)
integer(4) :: i
!#[1]## set input
n3ds = ocean_mesh%node
n2ds = ocean_mesh%npoi
n1s = ocean_mesh%n1

!#[1]##
do i=1,n2ds-1
 nzs(i)=n1s(i+1,1) - n1s(i,1)
end do
nzs(n2ds)=n3ds - n1s(n2ds,1) +1

write(*,*) "### MAKENZS5 END!! ###"
return
end subroutine makenzs5
!########################################## maketopbot6
!# modified on 2016.09.08
subroutine maketopbot6(n3top,n3bot,n3all,ntriall,ki23dptr,ki_mesh)
use mesh_type
implicit none
type(mesh),intent(in)  :: ki_mesh
integer(4),intent(in)  :: ki23dptr(2,ki_mesh%node)
integer(4),intent(out) :: n3top(ki_mesh%ntri),n3bot(ki_mesh%ntri),n3all(ki_mesh%ntri*2,3),ntriall
integer(4) :: ntrik,nodek,n2ds ! n2ds is total # of horizontal 2d ocean nodes; n3ds is # of nodes extruded 3d ocean mesh
integer(4)  :: n3k(ki_mesh%ntri,4),m3top(ki_mesh%ntri,3)
integer(4) :: i,j,n1,n2,n3,m3(3)
!#[1]## set input
ntrik = ki_mesh%ntri
nodek = ki_mesh%node
n3k(:,1) = ki_mesh%n3flag(:,2)
n3k(:,2:4) = ki_mesh%n3(:,1:3)

!#[2]## top surface
do i=1,ntrik !ntrik is # of land & ocean surface triangles
 n1=n3k(i,2);n2=n3k(i,3);n3=n3k(i,4)
 n3all(i,1:3) = (/ ki23dptr(1,n1),ki23dptr(1,n2),ki23dptr(1,n3) /)
 n3top(i) = i ! face id
end do

!#[3]## bottom surface
j=ntrik
do i=1,ntrik
 n1=n3k(i,2);n2=n3k(i,3);n3=n3k(i,4)
 if ( ki23dptr(1,n1) .ne. ki23dptr(2,n1) .or. ki23dptr(1,n2) .ne. ki23dptr(2,n2) .or. &
   &  ki23dptr(1,n3) .ne. ki23dptr(2,n3) ) then !if one or more points are not the same as surface node
  j=j+1
  n3all(j,1:3) = (/ ki23dptr(2,n1),ki23dptr(2,n2),ki23dptr(2,n3)/)
  n3bot(i) = j ! triangle id
 else
  n3bot(i) = i         ! triangle id
 end if
end do
ntriall = j

!write(*,*) "i=",i,"n3k(i,1)=",n3k(i,1),"n3k(i,2:4)=",n3k(i,2:4),"nzs(1:3)=",nzs(n1),nzs(n2),nzs(n3)
!write(*,*) "n3top(i,1:3)=",n3top(i,1:3)
!write(*,*) "n3bot(i,1:3)=",n3bot(i,1:3)
write(*,*) "### MAKETOPBOT6 END!! ###"
return
end subroutine maketopbot6
!########################################## premsh6_5
subroutine premsh6_5(xyz,n3d,n3all,ntriall,premsh)
implicit none
integer(4),intent(in) :: ntriall,n3d,n3all(ntriall,3)
real(8),intent(in)    :: xyz(3,n3d)
character(70),intent(in) :: premsh
integer(4) :: i,j
open(1,file=premsh)
write(1,'(a11)') "$MeshFormat"
write(1,'(a7)') "2.2 0 8"
write(1,'(a14)') "$EndMeshFormat"
write(1,'(a6)') "$Nodes"
write(1,*) n3d
do i=1,n3d
write(1,*) i,xyz(1:3,i)
end do
write(1,'(a9)') "$EndNodes"
! out oceanfile ### tetrahedron elements
write(1,'(a9)') "$Elements"
write(1,*) ntriall !ntrik
do i=1,ntriall
write(1,*) i,"  2   2   0   1 ",(n3all(i,j),j=1,3)
end do
write(1,'(a12)') "$EndElements"
close(1)
write(*,*) "### PREMSH6_5 END!! ###"
return
end subroutine premsh6_5
!########################################## findbryline7
subroutine findbryline7(ki_mesh,nlinbry,linbry_sb,l_line,ki23dptr,xyz,n3d,m4)
use mesh_type
use line_type
implicit none
type(mesh),     intent(in)  :: ki_mesh
type(line_info),intent(in)  :: l_line
integer(4),     intent(in)  :: ki23dptr(2,ki_mesh%node),n3d
real(8),        intent(in)  :: xyz(3,n3d)
integer(4),     intent(out) :: nlinbry(4),linbry_sb(2,4,100),m4(4,2)
integer(4)                  :: linbry(4,100)
! n2k(j,1) is the belonging old line # for line j, n2k(j,2:3) are the start and end node #
integer(4)                        :: nodek,nlink
integer(4)                        :: n2k(ki_mesh%nlin,3)
real(8),dimension(ki_mesh%node,2) :: xyk
real(8),dimension(4) :: d=(/1.d0, -1.d0, -1.d0, -1.d0/)
! nlinbry(i) is # of lines on i-th calculation boundary
! linbry(i,j) is line number of j-th line of i-th boundary
! 1 - 4 calculation boundaries are east, south, west, and north boundaries, respectively.
integer(4) :: line,i,j,k,l,ii,jj,n1,n2,n3,is,isnext,ilast,idir

!#[1]## set input
nodek      = ki_mesh%node
nlink      = ki_mesh%nlin
n2k(:,1  ) = ki_mesh%n2flag(:,2)
n2k(:,2:3) = ki_mesh%n2(:,1:2)
xyk(:,1)=ki_mesh%xyz(1,:) ! x, y of polygonki.msh
xyk(:,2)=ki_mesh%xyz(2,:)

!#[2]## make nlinbry and linebry
nlinbry(:)=0
ilast=0
is=1 ! starting node number, which is north east corner node
do i=1,4 ! 1: right; 2: bottom; 3: right; 4 : top boundary
 if ( i .eq. 1 .or. i .eq. 3 ) ii = 1 !   ( east, west)
 if ( i .eq. 2 .or. i .eq. 4 ) ii = 2 !   ( bottom, top)
 10 continue
 do j=1,nlink ! line loop  ; nlink is # of boundary lines in polygonki.msh
  n1=n2k(j,2);n2=n2k(j,3) ! set the start and end node of j-th line
  if ( (n1 .eq. is .or. n2 .eq. is ) .and. abs(xyk(n1,ii)-xyk(n2,ii)) .lt. 1.d-4 .and. ilast .ne. j ) then
  ! search line start from is
  nlinbry(i)=nlinbry(i)+1
  linbry(i,nlinbry(i))=j
  isnext=n2
  if ( is .eq. n2)      linbry(i,nlinbry(i))=-j
  if ( is .eq. n2)      isnext=n1
  is=isnext
  ilast=j
!write(*,*) "i=",i,"nlinbry(i)=",nlinbry(i)
!if ( linbry(i,nlinbry(i)) .gt. 0) write(*,*) n1,"->",n2
!if ( linbry(i,nlinbry(i)) .lt. 0) write(*,*) n2,"->",n1
!write(*,*) "xyk(n1,ii),xyk(n2,ii), delt=",xyk(n1,ii),xyk(n2,ii),(xyk(n1,ii)-xyk(n2,ii)),1.d-5
goto 10
end if
!if (n2k(j,1) .eq. 186)then
!write(*,*) "n2k(j,2:3)",n2k(j,2),n2k(j,3)
!write(*,*) "xyk(n1,ii),xyk(n2,ii), delt=",xyk(n1,ii),xyk(n2,ii),(xyk(n1,ii)-xyk(n2,ii))
!end if
end do
! if no lines selected for next linbry, move on next calculation boundary
end do

!#[3]## find up and down corner node
do i=1,4
 line=abs(linbry(i,1)) ; idir = line/linbry(i,1)
 n1 = n2k(line,(1-idir)/2 + 2) ! (1-idir)/2 = 0 for (idir = 1), = 1 for (idir = -1)
 m4(i,1)=ki23dptr(1,n1)
 m4(i,2)=ki23dptr(2,n1)
end do

!#[4]## find line id from l_line
  do i=1,4
    do j=1,nlinbry(i)
      do l=1,2 ! l=1: surface, l=2: bottom
      line = abs(linbry(i,j)) ; idir = linbry(i,j)/line
      n1=ki23dptr(l,n2k(line,2)) ; n2=ki23dptr(l,n2k(line,3))
      if ( n2 .lt. n1 ) then
       n3 = n1
	 n1 = n2
	 n2 = n3
	 idir = idir * (-1.)
      end if
      do k=l_line%line_stack(n1-1)+1,l_line%line_stack(n1)
       if (l_line%line_item(k) .eq. n2 ) goto 30
      end do
      goto 98
      30 continue
      linbry_sb(l,i,j) = k*idir ! line number for n3d line index is set
     end do
    end do
  end do


!#[5]## write nodes
if (.false.) then ! 2018.11.09
open(1,file="pre3dbry.geo")
write(1,*) "lc=500.0 ;"
 do i=1,n3d
 write(1,*)"Point(",i,")={",xyz(1,i),",",xyz(2,i),",",xyz(3,i),",lc} ;"
end do
!### write lines only on boundaries
jj=0
do i=1,4
 do j=1,nlinbry(i)
 jj=jj+1
 line=linbry_sb(1,i,j)
!  n1=n2k(abs(line),2);n2=n2k(abs(line),3)
!  write(*,*) "i=",i,"j=",j,"line=",line,"x1,y1=",(xyk(n1,k),k=1,2),"x2,y2=",(xyk(n2,k),k=1,2)
if ( line  .gt. 0) then
 write(1,*) "Line(",line,")={",l_line%line(1,line),",",l_line%line(2,line),"} ;"
else
 write(1,*) "Line(",-line,")={",l_line%line(1,-line),",",l_line%line(2,-line),"} ;"
end if
end do
end do
close(1)
end if ! 2018.11.09
if ( is .ne. 1) goto 99
write(*,*) "### FINDBRYLINE7 END!! ###"
return
99 write(*,*) "GEGEGE boundary round loop cannot be achieved."
write(*,*) "The end of node is", isnext
write(*,*) "x,y=",xyk(isnext,1),xyk(isnext,2)
stop
98 write(*,*) "GEGEGE corresponding boundary line in n3d was not found!!"
write(*,*) "line in n2k is",n1,"->",n2,"idir=",idir
stop
end subroutine findbryline7
!########################################## prepare7_5
subroutine prepare7_5(xyz,n3d,xyz3d,n3dn,g_meshpara)
use param_mesh
use mesh_type
implicit none
type(meshpara),intent(in) :: g_meshpara
integer(4),    intent(in) :: n3d,n3dn
real(8),       intent(in) :: xyz(3,n3d)
real(8),      intent(out) :: xyz3d(3,n3dn)
real(8)    :: zmax,zmin,x01,x02,y01,y02

!#[1]## set input
zmin = g_meshpara%zmin
zmax = g_meshpara%zmax
x01  = g_meshpara%xbound(1)
x02  = g_meshpara%xbound(4)
y01  = g_meshpara%ybound(1)
y02  = g_meshpara%ybound(4)

!! x3d, y3d, and z3d
xyz3d(1:3,1:n3d)=xyz(1:3,1:n3d) ! x,y,z are set in subroutine mknode4

xyz3d(1:3,n3d+1) = (/x02,y02,zmax/)
xyz3d(1:3,n3d+2) = (/x02,y01,zmax/)
xyz3d(1:3,n3d+3) = (/x01,y01,zmax/)
xyz3d(1:3,n3d+4) = (/x01,y02,zmax/)
xyz3d(1:2,n3d+5:n3d+8) = xyz3d(1:2,n3d+1:n3d+4)
xyz3d(3,  n3d+5:n3d+8) = zmin

write(*,*) "### PREPARE7_5 END!! ###"
return
end subroutine prepare7_5
!######################################### outpregeo8
subroutine outpregeo8(l_line,n3top,n3bot,xyz3d,n3d,pregeo,nlinbry,linbry_sb,&
&                     ki_mesh,g_meshpara,m4,ntri) ! 2018.11.18 ntri is added
use constants
use mesh_type
use param_mesh
use line_type
implicit none
type(line_info) ,intent(in)    :: l_line
type(mesh),      intent(in)    :: ki_mesh
type(meshpara),  intent(in)    :: g_meshpara
integer(4),      intent(in)    :: n3d
real(8),         intent(in)    :: xyz3d(3,n3d+8)
integer(4),      intent(in)    :: nlinbry(4),m4(4,2)
integer(4),      intent(in)    :: ntri   ! 2018.11.18
integer(4),      intent(inout) :: linbry_sb(2,4,100)
integer(4),      intent(inout) :: n3top(ntri),n3bot(ntri) ! 2018.11.18
! n23dlinbry(i,j,k) is start node(k=1) or end node(k=2) of j-th line on the i-th boundary
! linbry(i,j) is positive/negative line # corresponding to j-th line
! along the i-th calculation boundary;
! 1 to 4-th boundaries are right, bottom, left, and top boundary
! negative # indicates opposite direction of the original line
integer(4)               :: ntrik,nline,ntriall,nl
integer(4)               :: i,j,k,l,jj
integer(4),allocatable,dimension(:,:) :: n3line ! 2018.11.18
integer(4)               :: line(2,16),n43(2,4,3)
real(8),dimension(n3d+8) :: x3d,y3d,z3d
character(70)            :: pregeo
character(3)             :: lc
real(8)                  :: x1,x2,y1,y2
real(8)                  ::  sizein,sizebo ! 2017.06.28

!#[1]## set input
  ntrik   = ki_mesh%ntri
  x3d(:)  = xyz3d(1,:)
  y3d(:)  = xyz3d(2,:)
  z3d(:)  = xyz3d(3,:)
  x1      = g_meshpara%xbound(2) ! west
  x2      = g_meshpara%xbound(3) ! east
  y1      = g_meshpara%xbound(2) ! south
  y2      = g_meshpara%xbound(3) ! north
  nline   = l_line%nline
  ntriall = l_line%ntri
  allocate(n3line(ntriall,3)) ! 2018.11.18
  n3line  = l_line%n3line
  sizein  = g_meshpara%sizein
  sizebo  = g_meshpara%sizebo

!#[2]## Output points
  open(1,file=pregeo)

! Point Output
  write(1,*) "lc1=",sizein,";" ! 2017.06.28
  write(1,*) "lc2=",sizebo,";" ! 2017.06.28
  do i=1,n3d+8 ! n3d coordinates inputdo i=1,nodes
   lc="lc2"
   if ( x1 .lt. x3d(i) .and. x3d(i) .lt. x2 .and. y1 .lt. y3d(i) .and. y3d(i) .lt. y2 .and. i .le. n3d ) lc="lc1"
   write(1,*) "Point(",i,")={",x3d(i),",",y3d(i),",",z3d(i),","//lc(1:3)//"};"
  end do

!#[3]## Output lines
! additional 16 line                  !    4  top
  nl=nline
  line(1:2,1)  = (/ n3d + 1, n3d + 2 /) !   ---
  line(1:2,2)  = (/ n3d + 2, n3d + 3 /) ! 1 | | 3
  line(1:2,3)  = (/ n3d + 3, n3d + 4 /) !   ---
  line(1:2,4)  = (/ n3d + 4, n3d + 1 /) !    2
  line(1:2,5)  = (/ n3d + 1, m4(1,1) /) ! (x,y) = top    left  vertical line
  line(1:2,6)  = (/ n3d + 2, m4(2,1) /) ! (x,y) = bottom left  vertical line
  line(1:2,7)  = (/ n3d + 3, m4(3,1) /) ! (x,y) = bottom right vertical line
  line(1:2,8)  = (/ n3d + 4, m4(4,1) /) ! (x,y) = top    right vertical line
  line(1:2,9)  = (/ m4(1,2), n3d + 5 /)
  line(1:2,10) = (/ m4(2,2), n3d + 6 /)
  line(1:2,11) = (/ m4(3,2), n3d + 7 /)
  line(1:2,12) = (/ m4(4,2), n3d + 8 /)
  line(1:2,13) = (/ n3d + 5, n3d + 6 /)
  line(1:2,14) = (/ n3d + 6, n3d + 7 /)
  line(1:2,15) = (/ n3d + 7, n3d + 8 /)
  line(1:2,16) = (/ n3d + 8, n3d + 5 /)
  do i=1,nline
   write(1,*) "Line(", i,")={",l_line%line(1,i),",",l_line%line(2,i),"};" ! node id is not changed
  end do
  do i=1,16
   write(1,*) "Line(", nl+i,")={",line(1,i),",",line(2,i),"};"
  end do

!#[4]## output lineloop
!#[4-1]## triangles
  do i=1,ntriall
  write(1,*) "Line Loop(",i ,")={",n3line(i,1),",",n3line(i,2),",",n3line(i,3),"};"
  write(1,*) "Plane Surface(",i ,")={",i ,"};"
  end do

!#[4-1]## upper 5 and lower 5 faces
! for upper line loop
  nl=nline
  n43(1,1,1:3)=(/ -(nl+6),-(nl+1),nl+5/)
  n43(1,2,1:3)=(/ -(nl+7),-(nl+2),nl+6/)
  n43(1,3,1:3)=(/ -(nl+8),-(nl+3),nl+7/)
  n43(1,4,1:3)=(/ -(nl+5),-(nl+4),nl+8/)

! for lower line loop
  n43(2,1,1:3)=(/ nl+10,-(nl+13),- (nl+9)/)
  n43(2,2,1:3)=(/ nl+11,-(nl+14),-(nl+10)/)
  n43(2,3,1:3)=(/ nl+12,-(nl+15),-(nl+11)/)
  n43(2,4,1:3)=(/ nl+ 9,-(nl+16),-(nl+12)/)

! upmost face
  write(1,*) "Line Loop(",ntriall +1,")={",nl+1,",",nl+2,",",nl+3,",",nl+4,"};"
  write(1,*) "Plane Surface(",ntriall+1,")={",ntriall +1,"};"
!
  jj= ntriall + 1 ! line loop
  do l=1,2 ! upper and lower
  do i=1,4
  write(1,*) "Line Loop(",i+jj+(l-1)*4,")={",(n43(l,i,k),",",k=1,3) ! comma is the last
  write(1,*) (linbry_sb(l,i,k),",",k=1,nlinbry(i)-1),linbry_sb(l,i,nlinbry(i)),"};"
  write(1,*) "Plane Surface(",i+jj+(l-1)*4,")={",i+jj+(l-1)*4,"};"
  end do
  end do
! downmost face
  write(1,*) "Line Loop(",ntriall +10,")={",nl+13,",",nl+14,",",nl+15,",",nl+16,"};"
  write(1,*) "Plane Surface(",ntriall+10,")={",ntriall+10,"};"

!#[5]## output surface loop
  jj=ntriall
  write(1,*) "Surface loop(1)={",jj+1,",",jj+2,",",jj+3,",",jj+4,",",jj+5
  write(1,*)                  (",",n3top(i),i=1,ntrik),"};"
  write(1,*) "Volume(1)={1};"

  write(1,*) "Surface Loop(2)={",jj+6,",",jj+7,",",jj+8,",",jj+9,",",jj+10
  write(1,*)                  (",",n3bot(i),i=1,ntrik),"};"
  write(1,*) "Volume(2)={2};"



  100 continue
  close(1)
  write(*,*) "### OUTPREGEO8 END!! ###"
return
end subroutine outpregeo8

