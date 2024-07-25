program msh2spherical

use spherical
use mesh_type
use param_mesh

implicit none
type(meshpara) :: g_meshpara
integer(4) :: i,iout
!character(50) :: infile="em3d.msh", ofile="em3d_spherical.msh"
!character(50) :: infile="ocean.msh", ofile="ocean_spherical.msh"
character(70) :: infile="polygonki.msh", ofile="polygonki_spherical.msh"
type(mesh)    :: em_mesh
real(8) :: lonorigin,latorigin

!#[0]## read meshpara
CALL readmeshpara(g_meshpara)

!#[1]## set lonorigin, latorigin
lonorigin = g_meshpara%lonorigin
latorigin = g_meshpara%latorigin

!#[2]## read mesh
write(*,*) "input file : ", infile
CALL READMESH_TOTAL(em_mesh,infile) ! see m_mesh_type.f90
CALL ALLOCATEMESHSPHERICAL(em_mesh)

!#[3]## calculate xyzspherical
do i=1,em_mesh%node
 call xyz2lonlatalt(em_mesh%xyz(1:3,i),lonorigin,latorigin,em_mesh%lonlatalt(1:3,i))
 call lonlatalt2xyzspherical(em_mesh%lonlatalt(1:3,i),em_mesh%xyzspherical(1:3,i))
end do
write(*,*) "### calculation of xyz spherical end!! ###"

!#[4]## output
write(*,*) "output file : ", ofile
em_mesh%icoordinateflag=3 ! main coordinate is xyzspherical
iout=1
open(iout,file=ofile)
CALL MESHOUT(iout,em_mesh)
close(iout)

end program msh2spherical
