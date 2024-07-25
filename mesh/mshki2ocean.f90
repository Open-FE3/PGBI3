program mshki2ocean
use param_mesh ! see m_param_mesh.f90, added on Sep. 7, 2016
use mesh_type
implicit none
type(mesh)                          :: ki_mesh
type(meshpara)                      :: g_meshpara
integer(4)                          :: nodes,ntrik
integer(4)                          :: nodeseageo,nodegeo,naddedcorner,neles
character(70)                       :: ctlfile="mshki2ocean.ctl"
! this file provides nodeseageo and nodegeo as parameters
integer(4)                          :: maxseabry,maxlandbry,maxinsea
! these are max node number of ocean boundary, land boundary, and in the ocean in kifile
integer(4),allocatable,dimension(:) :: ki2oceanptr,ocean2kiptr
character(70)                       :: head,kifile,oceanfile

!#[1]## read meshparameters for head
CALL READMESHPARA(g_meshpara)
     head=g_meshpara%head
     kifile=head(1:len_trim(head))//"ki.msh"
     oceanfile=head(1:len_trim(head))//"_ki.msh"
  open(1,file=ctlfile)
   read(1,'(20x,i10)') nodeseageo
   read(1,'(20x,i10)') nodegeo
  close(1)
  naddedcorner = nodegeo - nodeseageo

!#[2]## read kifile
!call mshcount1(kifile,nodek,ntetk,ntrik,nlink,npoik)
!call readmsh2(kifile,xk,yk,zk,nodek,n4k,n3k,n2k,n1k,ntetk,ntrik,nelek,nlink,npoik)
 call READMESH_TOTAL(ki_mesh,kifile)
      ntrik=ki_mesh%ntri

!#[3]## find nodes for ocean mesh
call findnodenum3(ki_mesh,nodes,nodegeo,nodeseageo,maxseabry,maxlandbry,maxinsea,neles)

!#[4]## make ki2ocean pointer
allocate(ki2oceanptr(ntrik),ocean2kiptr(nodes))
call makenodeptr4(ki2oceanptr,ocean2kiptr,ntrik,nodeseageo,nodegeo,maxseabry,maxlandbry,maxinsea,nodes)

!#[5]## output oceanfile
call outocean(oceanfile,nodes,ocean2kiptr,ki2oceanptr,ki_mesh,neles,nodeseageo,ntrik) ! 2018.08.28

end program
!########################################## makenodeptr4
subroutine makenodeptr4(ki2oceanptr,ocean2kiptr,ntrik,nodeseageo,nodegeo,maxseabry,maxlandbry,maxinsea,nodes)
implicit none
integer(4),intent(in) :: ntrik,nodeseageo,nodegeo,maxseabry,maxlandbry,maxinsea,nodes
integer(4),intent(out) :: ki2oceanptr(ntrik),ocean2kiptr(nodes)
integer(4) :: i,j
ki2oceanptr(:)=99999999
j=0
do i=1,nodeseageo
j=j+1
ki2oceanptr(i)=j;ocean2kiptr(j)=i
end do
do i=nodegeo+1,maxseabry
j=j+1
ki2oceanptr(i)=j;ocean2kiptr(j)=i
end do
do i=maxlandbry+1,maxinsea
j=j+1
ki2oceanptr(i)=j;ocean2kiptr(j)=i
end do
if ( j .ne. nodes) goto 10
write(*,*) "### MAKENODEPTR4 END!! ###"
return
10 write(*,*) "GEGEGE! j =",j,"not equal to nodes=",nodes
stop
end subroutine makenodeptr4
!########################################## findnodenum3
subroutine findnodenum3(ki_mesh,nodes,nodegeo,nodeseageo,maxseabry,maxlandbry,&
&                                         maxinsea,neles)
use mesh_type
implicit none
type(mesh),            intent(in)     :: ki_mesh
integer(4),            intent(in)     :: nodegeo,nodeseageo
integer(4),            intent(out)    :: nodes,maxseabry,maxlandbry,maxinsea,neles
integer(4)                            :: nlink,ntrik
integer(4),allocatable,dimension(:,:) :: n2k,n3k
integer(4)                            :: i,nlineinsea,neleinsea

!#[1]## set input
 ntrik=ki_mesh%ntri
 nlink=ki_mesh%nlin
 allocate(n2k(nlink,3))
 allocate(n3k(ntrik,4))
 n2k(1:nlink,1  ) = ki_mesh%n2flag(1:nlink,2)
 n2k(1:nlink,2:3) = ki_mesh%n2(1:nlink,1:2)
 n3k(1:ntrik,1  ) = ki_mesh%n3flag(1:ntrik,2)
 n3k(1:ntrik,2:4) = ki_mesh%n3(1:ntrik,1:3)
 maxlandbry=0;maxseabry=0
 nlineinsea=0

 do i=1,nlink
  if (n2k(i,1) .le. nodeseageo ) then
   maxseabry=max(maxseabry,n2k(i,2),n2k(i,3))
   nlineinsea=nlineinsea+1
  else if (n2k(i,1) .gt. nodeseageo ) then
   maxlandbry=max(maxlandbry,n2k(i,2),n2k(i,3))
  end if
 end do
 if ( maxlandbry .lt. maxseabry ) maxlandbry=maxseabry
 write(*,*) "nodeseageo=",nodeseageo
 write(*,*) "maxseabry=",maxseabry
 write(*,*) "maxlandbry=",maxlandbry
 write(*,*) "nlineinsea=",nlineinsea

! find maxinsea
 maxinsea=0
 neleinsea=0
 do i=1,ntrik
  if (n3k(i,1) .eq. 1 ) then ! n3k(i,1) = 1 indicates the element is in the ocean
   maxinsea=max(maxinsea,n3k(i,2),n3k(i,3),n3k(i,4))
   !write(*,*) "maxinsea,n3k(i,2),n2k(i,3),n3k(i,4)",maxinsea,n3k(i,2),n3k(i,3),n3k(i,4)
   neleinsea=neleinsea+1
  end if
 end do
 write(*,*)"maxinsea=",maxinsea
 write(*,*)"neleinsea=",neleinsea
! calculate nodes
 nodes=nodeseageo+(maxseabry-nodegeo)+(maxinsea-maxlandbry)
! calculate neles
 neles=nodeseageo+nlineinsea+neleinsea
 write(*,*) "nodes=",nodes
 write(*,*) "neles=",neles
 write(*,*) "### FINDNODENUM3 END!! ###"

return
end subroutine findnodenum3
!########################################## outocean
subroutine outocean(oceanfile,nodes,ocean2kiptr,ki2oceanptr,ki_mesh,&
&                                 neles,nodeseageo,ntrik) ! 2018.08.28
use mesh_type
implicit none
type(mesh),             intent(in)     :: ki_mesh
integer(4),             intent(in)     :: ntrik
integer(4),             intent(in)     :: ki2oceanptr(ntrik)
integer(4),             intent(in)     :: ocean2kiptr(nodes)
integer(4),             intent(in)     :: nodes,neles,nodeseageo
character(70),          intent(in)     :: oceanfile
real(8),    allocatable,dimension(:)   :: xk,yk,zk
integer(4), allocatable,dimension(:,:) :: n2k,n3k
integer(4)                             :: i,j,nlink,nodek

!#[1]## set input
 nodek=ki_mesh%node
 nlink=ki_mesh%nlin
 allocate(xk(nodek),yk(nodek),zk(nodek))
 allocate(n2k(nlink,3))
 allocate(n3k(ntrik,4))
 xk=ki_mesh%xyz(1,:)
 yk=ki_mesh%xyz(2,:)
 zk=ki_mesh%xyz(3,:)
 n2k(1:nlink,1  ) = ki_mesh%n2flag(1:nlink,2)
 n2k(1:nlink,2:3) = ki_mesh%n2(1:nlink,1:2)
 n3k(1:ntrik,1  ) = ki_mesh%n3flag(1:ntrik,2)
 n3k(1:ntrik,2:4) = ki_mesh%n3(1:ntrik,1:3)

!#[2]## output oceanfile
 open(1,file=oceanfile)
 write(1,'(a11)') "$MeshFormat"
 write(1,'(a7)') "2.2 0 8"
 write(1,'(a14)') "$EndMeshFormat"
 write(1,'(a6)') "$Nodes"
 write(1,'(i10)') nodes
 do i=1,nodes
  write(1,'(i10,3g15.7)') i,xk(ocean2kiptr(i)),yk(ocean2kiptr(i)),zk(ocean2kiptr(i)) ! 2018.08.28
 end do
 write(1,'(a6)') "$EndNodes"
 write(1,'(a9)') "$Elements"
 write(1,'(i10)') neles
 do i=1,nodeseageo ! Point elements
  write(1,'(i10,a,2i10)')i," 15  2  0",i,i ! 2018.08.28
 end do
 j=nodeseageo
 do i=1,nlink ! Line elements
  if ( n2k(i,1) .le. nodeseageo) then
   j=j+1
   write(1,'(i10,a,3i10)')j," 1 2 0",n2k(i,1),ki2oceanptr(n2k(i,2)),ki2oceanptr(n2k(i,3))! 2018.08.28
  end if
 end do
 do i=1,ntrik
  if ( n3k(i,1) .eq. 1) then
   j=j+1
   write(1,'(i10,a,3i10)')j," 2 2 0 1",ki2oceanptr(n3k(i,2)),ki2oceanptr(n3k(i,3)),ki2oceanptr(n3k(i,4))! 2018.08.28
  end if
 end do
 if ( j .ne. neles ) goto 10
 write(1,'(a9)') "$EndElements"
 close(1)

!################################# added on May 17, 2016
 open(1,file="ki2oceanptr.dat")
 write(1,*) ntrik
 do i=1,ntrik
  write(1,'(i10)')ki2oceanptr(i)
 end do
 close(1)
!
write(*,*) "### OUTOCEAN END!! ###"
return
10 write(*,*) "GEGEGE j=",j,"and neles=",neles,"are not equal!!!"
stop
end subroutine outocean
