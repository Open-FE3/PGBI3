! coded on 2016.09.07
module mesh_relation
implicit none

type relate
 integer(4) :: nodeseageo
 integer(4) :: nodegeo
 integer(4) :: naddedcorner
 integer(4) :: maxseabry
 integer(4) :: maxlandbry
 integer(4) :: maxinsea
end type

contains

!##################################################  readgeorelation
subroutine readgeorelation(g_relate,infile)
type(relate),intent(inout) :: g_relate
character(50),intent(in) :: infile

  open(1,file=infile)
   read(1,'(20x,i10)') g_relate%nodeseageo
   read(1,'(20x,i10)') g_relate%nodegeo
  close(1)
  g_relate%naddedcorner = g_relate%nodegeo - g_relate%nodeseageo

end subroutine readgeorelation
end module