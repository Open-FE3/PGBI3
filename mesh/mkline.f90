! Coded on April 21, 2016
! to generate line file for 3dmesh.msh
program genline
use mesh_type  ! see m_mesh_type.f90
use line_type  ! see m_line_type.f90
implicit none
character(70)   :: filename="em3d.msh"
type(mesh)      :: g_mesh     ! see m_mesh_type.f90
type(line_info) :: g_line     ! see m_line_type.f90
integer(4)      :: nline, ntet,i,j

!#[1]## Mesh READ
  CALL MSHCOUNT1    (filename,g_mesh)  ! npoi = 0 makes no problems
  CALL ALLOCATE_MESH(g_mesh)
  CALL READMSH2     (filename,g_mesh)

!#[2]## Line information
!  CALL MKLINE(g_line, g_mesh%node, g_mesh%ntet, g_mesh%n4) ! make g_line
  CALL MKLINE_V2(g_line, g_mesh%node, g_mesh%ntet,4,g_mesh%n4) ! make g_line
  CALL MKN6  (g_line, g_mesh%node, g_mesh%ntet, g_mesh%n4) ! make g_line%n6line
  nline= g_line%nline
  ntet = g_line%ntet

!#[3]## Output lines and n6line
open(10,file="lineinfo.dat")
write(10,*) nline
write(10,'(2i10)') ((g_line%line(j,i),j=1,2),i=1,nline)
write(10,*) ntet
write(10,'(6i10)') ((g_line%n6line(i,j),j=1,6),i=1,ntet)
close(10)

end program genline
