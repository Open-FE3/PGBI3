! subroutine READLINE is added on April 21, 2016
! Coded on March 3, 2016
module line_type
implicit none

type line_info
integer(4) :: nline  ! # of line
integer(4) :: node   ! # of node
integer(4) :: ntet   ! # of tetrahedron
integer(4),dimension(:),  pointer :: line_stack ! line_stack(0:node)
integer(4),dimension(:,:),pointer :: line       !line(1:2,i) is start and end node
integer(4),dimension(:),  pointer :: line_item
integer(4),dimension(:,:),pointer :: n6line
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
character(50),intent(in) :: filename
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
integer(4),intent(in) :: nodeg, ntetg, n4(ntetg,4)
integer(4), parameter :: maxnline=3000000
integer(4) :: line_item_work(maxnline),line_stack(0:nodeg)
integer(4) :: i,j,k,l,n1,n2,nline,icount

!#[1]## cal line_item_work, l_line%line_stack
nline=0 ; line_stack(:)=0
do i=1,ntetg ! element loop
  do j=1,4      ! element node loop for starting node
    n1=n4(i,j)
    do k=1,4   ! element node loop for end node
      n2=n4(i,k)
      if ( n1 .lt. n2) then ! deal only with cases that start node id, n1, is less than end node, n2
      do l=line_stack(n1-1)+1, line_stack(n1)
        if ( n2 .eq. line_item_work(l) ) goto 100  ! check whether the line (n1, n2) exist or not
      end do
      line_item_work(line_stack(n1)+2:nline+1)=line_item_work(line_stack(n1)+1:nline)
      line_stack(n1:nodeg)=line_stack(n1:nodeg)+1
      if ( nline .eq. maxnline ) goto 999
        line_item_work( line_stack(n1) ) = n2 !line_stack(n1) is latest line id
        nline=nline+1
        100 continue
      end if
    end do
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

!######################################## MKN6LINE
! modeified on March 3, 2016
! make n6line, which store 6 lines composing tetrahedron, including direction info
! Coded on 2015.08.15 by T. MINAMI
subroutine MKN6(l_line,nodeg,ntetg,n4g) ! make
use fem_edge_util ! for kl (see m_fem_edge_util.f90)
implicit none
type(line_info),intent(inout) :: l_line
integer(4),intent(in) :: ntetg, n4g(ntetg, 4), nodeg
integer(4) :: n6line(ntetg,6), nline
integer(4) :: i,j,k, icount, lineid,idirection, n11,n12
integer(4) :: line_stack(0:l_line%node), line_item(1:l_line%nline)
!#[0] copy
nline=l_line%nline
line_stack(0:l_line%node)=l_line%line_stack(0:l_line%node)
line_item(1:nline)       =l_line%line_item(1:nline)

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

end module line_type