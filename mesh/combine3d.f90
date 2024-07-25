! Modified on 2016.09.10
! "meshpara.h" -> m_param_mesh.f90
!#
! line 29 is added on May 30, 2016
! to decrease jumps in flatness of elements
!==
! Coded by T. MINAMI on 2014.07.31
! This program combine ocean.msh (tet4) and pre3d.msh (tet4) to em3d.msh (tet4)
! This program assumed that pre3d.geo includes the whole ocean.msh nodes
! prior to the nodes for air and land
program combine3d
use param_mesh
use mesh_type
implicit none
!include "meshpara.h" !inserted on may 30, 2016, only to use zratio
type(meshpara)                        :: g_meshpara
type(mesh)                            :: pre3d_mesh,ocean_mesh
integer(4)                            :: nodes,ntets,ntris,nlins,npois ! ocean.msh
integer(4)                            :: nodea,nteta,ntria
integer(4)                            :: nlina,npoia ! for pre3d.msh
real(8),   allocatable,dimension(:,:) :: xyzs, xyza
integer(4),allocatable,dimension(:,:) :: n4s,n3s,n2s
integer(4),allocatable,dimension(:,:) :: n4a,n3a,n2a
integer(4),allocatable,dimension(:)   :: n1s, n1a
character(70)                         :: oceanfile="ocean.msh"
character(70)                         :: pre3dfile="pre3d.msh"
character(70)                         :: em3dfile="em3d.msh"
integer(4)                            :: i ! added on May 30, 2016
real(8)                               :: a,b ! added on May 31, 2016
!#####################

!#[0]##  read parameters
 call READMESHPARA(g_meshpara)

!#[1]##  read two mesh files
 call READMESH_TOTAL(ocean_mesh, oceanfile)
 call READMESH_TOTAL(pre3d_mesh, pre3dfile)

!#[2]##     read airlandfile
!======================= the following line is to decrease jumps in element
!                        flatness between ocean layer and the other regions
!do i=1,nodea
!    if (  l_up .le. xyza(3,i) ) then
!      xyza(3,i) = xyza(3,i) - (l_up -l0_up)
!    else if ( 0.d0 .lt. xyza(3,i) .and.  xyza(3,i) .le. l_up )then
!      a=-1./l0_up*(l_up/l0_up -1.d0) ; b=1.d0 -2.*a*l0_up
!      xyza(3,i)= - b/2./a - dsqrt(xyza(3,i)/a +(b/2./a)**2.d0)
!    else if ( l_down .lt. xyza(3,i) .and.  xyza(3,i) .le. 0.d0 )then
!      a=-1./l0_down*(l_down/l0_down -1.d0) ; b=1.d0 -2.*a*l0_down
!      xyza(3,i)= - b/2./a + dsqrt(xyza(3,i)/a +(b/2./a)**2.d0)
!    else if ( xyza(3,i) .le. l_down )  then
!      xyza(3,i)=xyza(3,i) - (l_down - l0_down)
!    else
!    write(*,*) "GEGEGE there are z not corresponding any zgrid xyza(3,i)=",xyza(3,i)
!    stop
!end if
!end do
!#####################

!#[3]##     out put em3dfile
 CALL OUTEM3D(em3dfile,pre3d_mesh,ocean_mesh)
end program combine3d
!##

!############################################  OUTEM3d
subroutine OUTEM3D(em3dfile,pre3d_mesh,ocean_mesh)
use mesh_type
implicit none
character(70),         intent(in)     :: em3dfile
type(mesh),            intent(in)     :: pre3d_mesh,ocean_mesh
integer(4)                            :: nodea,nteta,ntets,npois
integer(4),allocatable,dimension(:,:) :: n4a
integer(4),allocatable,dimension(:,:) :: n4s
integer(4),allocatable,dimension(:)   :: n1s
real(8),   allocatable,dimension(:,:) :: xyza
integer(4)                            :: i, j, l
integer(4)                            :: ishift

!#[1]## set input
nodea      = pre3d_mesh%node
nteta      = pre3d_mesh%ntet
npois      = ocean_mesh%npoi
ntets      = ocean_mesh%ntet ! 2018.08.28
allocate( n4s(ntets,5) )     ! 2018.08.28
allocate( n4a(nteta,5) )     ! 2018.08.28
allocate( n1s(npois  ) )     ! 2018.08.28
allocate( xyza(3,nodea))     ! 2018.08.28
xyza       = pre3d_mesh%xyz
n4a(:,1  ) = pre3d_mesh%n4flag(:,2)
n4a(:,2:5) = pre3d_mesh%n4(:,1:4)
ntets = ocean_mesh%ntet
n1s(:)     = ocean_mesh%n1(:,1)
n4s(:,1  ) = ocean_mesh%n4flag(:,2)
n4s(:,2:5) = ocean_mesh%n4(:,1:4)

!#[2]## output
open(1,file=em3dfile)
write(1,'(a11)') "$MeshFormat"
write(1,'(a7)') "2.2 0 8"
write(1,'(a14)') "$EndMeshFormat"
write(1,'(a6)') "$Nodes"
write(1,'(i10)') nodea
do i=1,nodea
write(1,'(i10,3g15.7)') i, (xyza(j,i),j=1,3) ! 2018.08.28
end do
write(1,'(a9)') "$EndNodes"
! out oceanfile ### tetrahedron elements
write(1,'(a9)') "$Elements"
write(1,'(i10)') npois+ntets+nteta  !ntets+nteta
do i=1,npois
write(1,'(i10,a,i10)') i,"  15   2   0   1", n1s(i) ! 1 for ocean
end do
ishift=npois
do i=1,ntets
write(1,'(i10,a,4i10)') ishift+i,"  4   2   0   2", (n4s(i,l),l=2,5) ! 2 for ocean
end do
ishift=npois+ntets
do i=1,nteta
if (n4a(i,1) .eq. 2 ) n4a(i,1) = 3 ! 1 for air, 3 for land < CHECK!!
write(1,'(i10,a,5i10)') i+ishift,"  4   2   0",(n4a(i,l),l=1,5) ! n4a(1)=2 for air, 3 for land
end do
write(1,'(a12)') "$EndElements"
close(1)
write(*,*) "### OUTEM3D END ###"
end subroutine outem3d


