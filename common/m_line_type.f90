! subroutine READLINE is added on April 21, 2016
! Coded on March 3, 2016
module line_type
use fem_edge_util ! added on May 24, 2016
implicit none

type line_info
integer(4) :: nline  ! # of line
integer(4) :: node   ! # of node
integer(4) :: ntet   ! # of tetrahedron
integer(4) :: ntri   ! # of triangle
integer(4),dimension(:),  pointer :: line_stack ! line_stack(0:node)
integer(4),dimension(:,:),pointer :: line       !line(1:2,i) is start and end node
integer(4),dimension(:),  pointer :: line_item
integer(4),dimension(:,:),pointer :: n6line
integer(4),dimension(:,:),pointer :: n3line   ! added on 2016.09.08
     ! line_stack(i) stores total line numbers until i-th node line group,
     !                             where all of lines start with i-th node
     ! line_item(j) is the end node id for j-th line
     ! NOTE that lines stored in line_stack are ordered by the starting node ids, 
     !                             not related to the source current
     ! line(1:2,i) is the start and end node id for i-th line
     ! line_item(i) is the same as line_item_work(i) except size
end type

contains
!######################################## READLINE added on April 21, 2016
subroutine READLINE(filename,g_line)
implicit none
character(70),intent(in) :: filename
type(line_info),intent(out) :: g_line
integer(4) :: nline,ntet,i,j
!
  open(10,file=filename)
   read(10,*) nline
   allocate(g_line%line(2,nline))
   read(10,'(2I10)') ((g_line%line(j,i),j=1,2),i=1,nline)
   read(10,*) ntet
   allocate(g_line%n6line(ntet,6))
   read(10,'(6I10)') ((g_line%n6line(i,j),j=1,6),i=1,ntet)
  close(10)
  g_line%nline=nline
  g_line%ntet =ntet
  write(*,*) "### READLINE END!! ###"
return
end subroutine READLINE

!######################################## MKLINE
! line_stack(i) - line_stack(i-1) corresponds to the number of lines belonging to i-th node
! line_item(j) is the ending node id of j-th line
! NOTE that the line ids are ordered by the starting nodes.
! nline is total number of edge lines, including the edges corresponding to the source cable
subroutine MKLINE(l_line, nodeg, ntetg, n4) ! line_stack,line_item, nline generated
implicit none
type(line_info),intent(out) :: l_line
integer(4),     intent(in)  :: nodeg, ntetg, n4(ntetg,4)
integer(4),     parameter   :: maxnline=3000000
integer(4)                  :: line_item_work(maxnline),line_stack(0:nodeg)
integer(4)                  :: i,j,k,l,n1,n2,nline,icount,n1o,n2o

!#[1]## cal line_item_work, l_line%line_stack
nline=0 ; line_stack(:)=0
do i=1,ntetg ! element loop
  do j=1,6      ! element node loop for starting node
    n1o=n4(i,kl(j,1))
    n2o=n4(i,kl(j,2))
    n1=min(n1o,n2o) ! start node
    n2=max(n1o,n2o) ! end node
      do l=line_stack(n1-1)+1, line_stack(n1)
        if ( n2 .eq. line_item_work(l) ) goto 100  ! check whether the line (n1, n2) exist or not
      end do
      line_item_work(line_stack(n1)+2:nline+1)=line_item_work(line_stack(n1)+1:nline)
      line_stack(n1:nodeg)=line_stack(n1:nodeg)+1
      if ( nline .eq. maxnline ) goto 999
	line_item_work( line_stack(n1) ) = n2 !line_stack(n1) is latest line id
	nline=nline+1
	100 continue
   end do
end do
write(*,*) "nline=",nline

!do i=1,nodeg
!write(*,*) "i=",i,"line_stack(i)=",line_stack(i),"line_item(i)=",(line_item(j),j=line_stack(i-1)+1,line_stack(i))
!end do

!#[2]## set line_item and nline
     allocate( l_line%line(2,nline), l_line%line_item(nline))
     allocate( l_line%line_stack(0:nodeg) )
     l_line%line_stack(0:nodeg)=line_stack(0:nodeg)
     l_line%line_item(1:nline)=line_item_work(1:nline)
     l_line%nline=nline
     l_line%node=nodeg
     l_line%ntet=ntetg

!#[3]## set line
icount=0
do i=1,nodeg
  do j=line_stack(i-1)+1, line_stack(i)
    icount=icount+1
    l_line%line(1:2,icount)=(/ i,l_line%line_item(j) /)
    !write(*,*) "icount=",icount, "line(1:2,icount)=",line(1:2,icount)
  end do
end do
if ( icount .ne. nline ) goto 99

write(*,*) "### MKLINE END ###"
return
999 continue
write(*,*) "GEGEGE! # of line exceeds maxnline =", maxnline
write(*,*) "ntetg=",ntetg,"i=",i
stop
99 write(*,*) "GEGEGE! icount=", icount,"is not equal to nline=",nline
stop

end subroutine MKLINE
!######################################## MKLINE_V2
!# Modified on Sep. 9, 2016, enod is adopted
!# Coded on May 31, 2016
subroutine MKLINE_V2(l_line, nodeg, ne, enod, n4) ! line_stack,line_item, nline generated
implicit none
type(line_info),       intent(out)    :: l_line
integer(4),            intent(in)     :: nodeg, ne, enod
integer(4),            intent(in)     :: n4(ne,enod)
integer(4),            parameter      :: maxnline=3500000
integer(4),            parameter      :: maxnele=100
integer(4),allocatable,dimension(:)   :: line_item_work, line_stack
integer(4),allocatable,dimension(:)   :: nele
integer(4)                            :: i,j,k,l,n1,n2,nline,icount,n1o,n2o,iele
integer(4),allocatable,dimension(:,:) :: node2ele
! node2ele(j,i) is element id of j-th element of i-th node
! nele(i) is # of element including i-th node

!#[0]## gen node-element connectivity data
allocate( node2ele(maxnele,nodeg) ) ! 2018.08.28
allocate( line_item_work(maxnline)) ! 2018.08.28
allocate( line_stack(0:nodeg)     ) ! 2018.08.28
allocate( nele(nodeg)             ) ! 2018.08.28
nele(:)=0
do i=1,ne
 do j=1,enod
  n1=n4(i,j)
  if (maxnele .eq. nele(n1)) goto 998
  nele(n1)=nele(n1)+1     ! increase nele(n1) by 1
  node2ele(nele(n1),n1)=i ! set i-th element to element list of n1 node
 end do
end do
write(*,*) "node-element connectivity end!"

!#[1]## cal line_item_work, l_line%line_stack
nline=0 ; line_stack(:)=0
do n1=1,nodeg ! element loop
  if (mod(n1,100000) .eq. 0) write(*,*) "n1=",n1,"nodeg=",nodeg
  line_stack(n1)=line_stack(n1-1)
  do j=1,nele(n1)     ! element loop for i-th node
    iele=node2ele(j,n1) ! element id
!    write(*,*) "iele=",iele,"n4(iele,1:4)=",n4(iele,1:4)
    do k=1,enod ! node loop
    n2=n4(iele,k)
!    write(*,*) "n1=",n1,"n2=",n2
    if ( n2 .le. n1 ) goto 100 ! below n1 < n2
      do l=line_stack(n1-1)+1, line_stack(n1)
        if ( n2 .eq. line_item_work(l) ) goto 100  ! check whether the line (n1, n2) exist or not
      end do
      if ( nline .eq. maxnline ) goto 999
      line_stack(n1)=line_stack(n1)+1
	line_item_work( line_stack(n1) ) = n2 !line_stack(n1) is latest line id
	nline=nline+1
	100 continue
     end do ! 4 node loop
   end do   ! element loop
!   write(*,*) "node=",n1,"# of lines",&
!       &      line_stack(n1)-line_stack(n1-1),"nele=",nele(n1)
!   write(*,*) "line_item_work=",line_item_work(line_stack(n1-1)+1:line_stack(n1))
!if (n1 .eq. 2) stop
end do      ! node loop
write(*,*) "nline=",nline

!do i=1,nodeg
!write(*,*) "i=",i,"line_stack(i)=",line_stack(i),"line_item(i)=",(line_item(j),j=line_stack(i-1)+1,line_stack(i))
!end do

!#[2]## set line_item and nline
     allocate( l_line%line(2,nline), l_line%line_item(nline))
     allocate( l_line%line_stack(0:nodeg) )
     l_line%line_stack(0:nodeg)=line_stack(0:nodeg)
     l_line%line_item(1:nline)=line_item_work(1:nline)
     l_line%nline=nline
     l_line%node=nodeg
     if ( enod .eq. 4 ) l_line%ntet=ne
     if ( enod .eq. 3 ) l_line%ntri=ne

!#[3]## set line
icount=0
do i=1,nodeg
  do j=line_stack(i-1)+1, line_stack(i)
    icount=icount+1
    l_line%line(1:2,icount)=(/ i,l_line%line_item(j) /)
    !write(*,*) "icount=",icount, "line(1:2,icount)=",line(1:2,icount)
  end do
end do
if ( icount .ne. nline ) goto 99

write(*,*) "### MKLINE_V2 END ###"
return
99 continue
write(*,*) "icount=",icount,"is not equal to nline=",nline
stop
999 continue
write(*,*) "GEGEGE! # of line exceeds maxnline =", maxnline
if ( enod .eq. 4 ) write(*,*) "ntet=",ne,"i=",i
if ( enod .eq. 3 ) write(*,*) "ntri=",ne,"i=",i
stop
998 continue
write(*,*) "GEGEGE! # maxnele=",maxnele,"too small for node-element connectivity data"
stop
end subroutine MKLINE_V2

!######################################## MKN6LINE
! modeified on March 3, 2016
! make n6line, which store 6 lines composing tetrahedron, including direction info
! Coded on 2015.08.15 by T. MINAMI
subroutine MKN6(l_line,nodeg,ntetg,n4g) ! make
use fem_edge_util ! for kl (see m_fem_edge_util.f90)
implicit none
type(line_info),       intent(inout)  :: l_line
integer(4),            intent(in)     :: ntetg,nodeg
integer(4),            intent(in)     :: n4g(ntetg,4)
integer(4),allocatable,dimension(:,:) :: n6line
integer(4)                            ::  nline
integer(4)                            :: i,j,k, icount, lineid,idirection, n11,n12
integer(4),allocatable,dimension(:)   :: line_stack, line_item
integer(4)                            :: node

!#[0] set
 nline      = l_line%nline
 node       = l_line%node       ! 2018.08.28
 allocate( line_stack(0:node) ) ! 2018.08.28
 allocate( line_item(1:nline) ) ! 2018.08.28
 allocate( n6line(ntetg,6)    ) ! 2018.08.28
 line_stack = l_line%line_stack
 line_item  = l_line%line_item

!#[1] cal n6line
do i=1,ntetg ! element loop
 do j=1,6    ! line loop
  !#[1-1] set 2 nodes
  n11=n4g(i, kl(j,1)) ; n12=n4g(i, kl(j,2)) ; idirection=1
  if ( n11 .gt. n12 ) then
   n12=n4g(i, kl(j,1)) ; n11=n4g(i, kl(j,2)) ; idirection=-1
  end if
  icount=line_stack(n11-1)
  !#[1-2] search for the node
  do k=line_stack(n11-1)+1, line_stack(n11)
   icount=icount+1
   if ( line_item(k) .eq. n12 ) goto 100
  end do
  goto 99
  100 continue
  lineid=icount
  n6line(i,j)=idirection*lineid
 end do
 !write(*,*) "i=",i,"ntetg=",ntetg,"n6line(i,1:6)=",n6line(i,1:6)
end do

!#[2]## set n6line
allocate(l_line%n6line(ntetg,6) )
l_line%n6line(1:ntetg,1:6)=n6line(1:ntetg,1:6)

write(*,*) "### MKN6LINE END!! ###"
return
99 write(*,*) "GEGEGE! i=",i,"/ntetg=",ntetg," line n11, n12=",n11,n12, "was not found!"
stop
end subroutine MKN6

!######################################## MKN3LINE
!  Coded on 2016.09.08 by T. MINAMI, based on MKN6
subroutine MKN3(l_line,nodeg,ntri,n3) ! make
use fem_edge_util ! for lm (see m_fem_edge_util.f90)
implicit none
type(line_info),    intent(inout)   :: l_line
integer(4),         intent(in)      :: ntri, n3(ntri, 3), nodeg
integer(4)                          :: n3line(ntri,3), nline
integer(4)                          :: i,j,k, icount, lineid,idirection, n11,n12
integer(4),allocatable,dimension(:) :: line_stack,line_item
integer(4)                          :: node ! 2018.11.18

!#[0] copy
  node  = l_line%node
  nline = l_line%nline ! 2018.11.18
  allocate(line_stack(0:node),line_item(1:nline)) ! 2018.11.18
  line_stack(0:l_line%node)=l_line%line_stack(0:l_line%node)
  line_item(1:nline)       =l_line%line_item(1:nline)

!#[1] cal n3line
do i=1,ntri ! element loop
 do j=1,3    ! line loop
  !#[1-1] set 2 nodes
  n11=n3(i, lm(j,1)) ; n12=n3(i, lm(j,2)) ; idirection=1
  if ( n11 .gt. n12 ) then
   n12=n3(i, lm(j,1)) ; n11=n3(i, lm(j,2)) ; idirection=-1
  end if
  icount=line_stack(n11-1)
  !#[1-2] search for the node
  do k=line_stack(n11-1)+1, line_stack(n11)
   icount=icount+1
   if ( line_item(k) .eq. n12 ) goto 100
  end do
  goto 99
  100 continue
  lineid=icount
  n3line(i,j)=idirection*lineid
 end do
 !write(*,*) "i=",i,"ntetg=",ntetg,"n6line(i,1:6)=",n6line(i,1:6)
end do

!#[2]## set n6line
allocate(l_line%n3line(ntri,3) )
l_line%n3line(1:ntri,1:3)=n3line(1:ntri,1:3)

write(*,*) "### MKN3LINE END!! ###"
return
99 write(*,*) "GEGEGE! i=",i,"/ntri=",ntri," line n11, n12=",n11,n12, "was not found!"
stop
end subroutine MKN3
end module line_type
