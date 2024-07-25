program ocean_l2q
implicit none
integer(4) :: i,j! nsele is # of triangle corresponding to ocean surfaces
integer(4) :: ntris,nlins,n2ds
integer(4) :: node_num1, node_num2, nodesurf_num
integer(4) :: element_num
real(8),allocatable,dimension(:,:) :: node_xyz1
real(8),allocatable,dimension(:,:) :: node_xyz2
integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node1
integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node2
integer ( kind = 4 ), allocatable, dimension ( :, : ) :: surf_node
integer(4),allocatable,dimension(:,:) :: n3s,n2s
integer(4),allocatable,dimension(:) :: n1s
integer ( kind = 4 ), allocatable, dimension ( :, : ) :: edge_data
integer(4) :: nodesurf_num2
character(50) :: oceanfile="ocean.msh",tet10file='ocean_t10.msh'
!###
! read ocean.msh
CALL MSHCOUNT1(oceanfile, node_num1,element_num, ntris,nlins,n2ds)
          allocate ( node_xyz1(3,node_num1) )
          allocate ( element_node1(4,element_num), n3s(ntris,4),n2s(nlins,3),n1s(n2ds))
          ! npoi = 0 makes no problems
CALL READMSH2(oceanfile, node_xyz1, node_num1, element_node1,&
                           & n3s,n2s,n1s,element_num, ntris,i,nlins,n2ds)
          !## calculate edge data and node_num2 : total number of nodes in tet10 mesh
          allocate ( edge_data(1:5,1:6*element_num) )
CALL tet_mesh_order4_to_order10_size ( element_num, element_node1, node_num1, &
          &  edge_data, node_num2 )
          write ( *, '(a,i8)' ) '  Number of linear nodes = ', node_num1
          write ( *, '(a,i8)' ) '  Number of quadratic nodes = ', node_num2
!##  calculate element_node2
          allocate ( node_xyz2(3,node_num2) )
          allocate ( element_node2(10,element_num) )
CALL tet_mesh_order4_to_order10_compute ( element_num, element_node1, node_num1, &
           & node_xyz1, edge_data, element_node2, node_num2, node_xyz2 )
!## find surface node
          allocate( surf_node(3,node_num2) )
CALL surfacenode3( n2ds,node_num1,node_num2,node_xyz2,n1s,&
	    &   surf_node, nodesurf_num2)
!##  output
CALL writetet10(node_num2,element_num, node_xyz2, element_node2,tet10file,&
&   nodesurf_num2, surf_node(1:3,1:nodesurf_num2))
end program ocean_l2q
!########################################## mshcount1
subroutine mshcount1(mshfile,node,ntet,ntri,nlin,npoi)
implicit none
integer(4),intent(out) :: node,ntet,ntri,nlin,npoi
integer(4) :: i,j,k,itype,ii,jj,kk,etype,nele,n44(4)
character(50),intent(in) :: mshfile
open(1,file=mshfile)
do i=1,4
read(1,*)
end do
read(1,*) node
!write(*,*) "node=",node
do i=1,node+2 ! include "$EndNodes", "$Elements" lines
read(1,*)
end do
read(1,*) nele ! total number of elements
npoi=0;nlin=0;ntri=0;ntet=0
do i=1,nele
read(1,*) j,itype,ii,jj,kk,(n44(k),k=1,etype(itype))
if ( itype .eq. 15) then   ! Point element
npoi=npoi+1
else if ( itype .eq. 1) then ! Line element
nlin=nlin+1
else if ( itype .eq. 2) then   ! Triangle element
ntri=ntri+1
else if ( itype .eq. 4) then ! Tetrahedral element
ntet=ntet+1
end if
end do
close(1)
!write(*,*) "# of Point elements is",npoi
!write(*,*) "# of Line elements is",nlin
!write(*,*) "# of Triangle elements is",ntri
!write(*,*) "# of Tetrahedron elements is",ntet
write(*,*) "### COUNT ", mshfile(1:len_trim(mshfile))," END!! ###"
return
end subroutine mshcount1
!##########################################  readocean2
subroutine readmsh2(mshfile,xyz,node,n4,n3,n2,n1,ntet,ntri,nkk,nlin,npoi)
implicit none
integer(4),intent(in) :: node,ntet,ntri,nlin,npoi
real(8),dimension(3,node),intent(out) :: xyz
integer(4),intent(out) :: n4(4,ntet),n3(ntri,4),n2(nlin,3),n1(npoi)
integer(4) :: itet,itri,ilin,ipoi,nele,inode
integer(4) :: i,j,k,ii,jj,kk,itype,etype,nkk
character(50) :: mshfile
integer(4),dimension(4) :: n44
open(1,file=mshfile)
!# [1] ### skip header
do i=1,4
read(1,*)
end do
!# [2] ### read node
read(1,*) inode
write(*,*) "# of nodes (node) =",inode
if ( inode .gt. node) goto 999
do i=1,node
read(1,*) j,(xyz(j,i),j=1,3)
end do
read(1,*)
!# [3] ### read elements
read(1,*) ! skip the starting line, "$Elements"
read(1,*) nele
write(*,*)"# of elements (nele)=",nele
ipoi=0;ilin=0;itri=0;itet=0;nkk=0
do i=1,nele
read(1,*) j,itype,ii,jj,kk,(n44(k),k=1,etype(itype))
if ( itype .eq. 15) then   ! Point element
ipoi=ipoi+1
n1(ipoi)=n44(1)
else if ( itype .eq. 1) then ! Line element
ilin=ilin+1
n2(ilin,1)=kk
n2(ilin,2:3)=n44(1:2)
else if ( itype .eq. 2 ) then ! Triangle element
itri=itri+1
n3(itri,1)=kk
n3(itri,2:4)=n44(1:3)
if (kk .eq. 1) nkk=nkk+1 ! count triangle group 1
else if ( itype .eq. 4) then ! Tetrahedral element
itet=itet+1
n4(1:4,itet)=n44(1:4)
end if
end do
write(*,*) "# of Point elements is",npoi
write(*,*) "# of Line elements is",nlin
write(*,*) "# of Triangle elements is",ntri," # of traiangle group1 is",nkk
write(*,*) "# of Tetrahedron elements is",ntet
close(1)
write(*,*) "### READ ", mshfile(1:len_trim(mshfile))," END!! ###"
return
999 write(*,*) "GEGEGE node .ne. inode, node=",node,"inode=",inode
stop
end subroutine readmsh2
!###########################################  function etype
function etype(itype)
implicit none
integer(4) :: itype,etype
etype=0
if (itype .eq. 15) etype=1! one point node
if (itype .eq. 1 ) etype=2 ! line
if (itype .eq. 2 ) etype=3 ! triangle
if (itype .eq. 4 ) etype=4 ! tetrahedron
return
end function
!############################################### surfacenode3
subroutine surfacenode3( nodesurf_num1,node_num1,node_num2,node_xyz2,n1s,&
& surf_node, nodesurf_num2)
implicit none
integer(4),intent(in) :: nodesurf_num1
integer(4),intent(in) :: node_num1,node_num2
real(8),dimension(3,node_num2),intent(in) :: node_xyz2
integer(4),dimension(nodesurf_num1),intent(in) :: n1s
integer(4),intent(out) :: nodesurf_num2
integer(4),dimension(3,node_num2),intent(out) :: surf_node
integer(4) :: i, j, k
real(8) :: x, y, z, z2, z3
! surf_node(1,j) is the top node # of j-th horizontal node column
! surf_node(2,j) is node # just beneath the top layer at j-th node column
! surf_node(3,j) is node # at third from the top at j-th horizontal node colmun
!## surf_node(1,j) top surface / initialize the surf_node(2:3,j)=surf_node(1,j)
surf_node(1,1:nodesurf_num1)=n1s(1:nodesurf_num1) ! surface vertices
k=nodesurf_num1
do i=node_num1+1,node_num2
  if ( abs(node_xyz2(3,i)) .lt. 1.d-10 ) then
    k=k+1
    surf_node(1,k)=i
  end if
end do
nodesurf_num2=k
write(*,*) "Nomber of surface node linear tet is",nodesurf_num1
write(*,*) "Nomber of surface node quadratic tet is",nodesurf_num2
!## search for surf_node(2:3,j) : 2nd, 3rd node # of j-th horizontal node column
surf_node(2:3,:)=-1 !
do i=1,node_num2
if ( node_xyz2(3,i) .lt. 0.d0 ) then ! not the surface node pass
x=node_xyz2(1,i)
y=node_xyz2(2,i)
z=node_xyz2(3,i)
do j=1,nodesurf_num2
if ( node_xyz2(1,surf_node(1,j)) .eq. x .and. node_xyz2(2,surf_node(1,j)) .eq. y)then
   if ( surf_node(2,j) .eq. -1 ) then ! in the case that surf_node(2:3) not given
     surf_node(2:3,j)=i
   else ! in the case that surf_node(2:3) are already given
      z2=node_xyz2(3,surf_node(2,j))
      z3=node_xyz2(3,surf_node(3,j))
      if ( z2 .eq. z3 .and. z .gt. z2 ) surf_node(2,j)=i
      if ( z2 .eq. z3 .and. z3 .gt. z ) surf_node(3,j)=i
      if ( z2 .ne. z3 .and. z .gt. z2 ) then
        surf_node(3,j)=surf_node(2,j)
        surf_node(2,j)=i
      end if
      if ( z2 .ne. z3 .and. z2 .gt. z .and. z .gt. z3 ) surf_node(3,j)=i
      if ( node_xyz2(3,surf_node(2,j)) .eq. node_xyz2(3,surf_node(3,j))) goto 99
    end if
end if
end do
end if
end do
write(*,*) "### SURFACENODE3 END ###"
return
99 write(*,*) "GEGEGE correct surf_node is not given when i,j=",i,j
     write(*,*) "x,y,z=",x,y,z
     write(*,'(a,3i7)') "surf_node(1:3,j)=",(surf_node(k,j),k=1,3)
     write(*,*) "x for surf_node(1:3,j) is",(node_xyz2(1,surf_node(k,j)),k=1,3)
     write(*,*) "y for surf_node(1:3,j) is",(node_xyz2(2,surf_node(k,j)),k=1,3)
     write(*,*) "z for surf_node(1:3,j) is",(node_xyz2(3,surf_node(k,j)),k=1,3)
stop
end subroutine surfacenode3
!###############################################   writete10
subroutine writetet10(node_num2,element_num, node_xyz2, element_node2,tet10file,nodesurf_num2, surf_node)
implicit none
integer(4) :: node_num2, element_num
real(8),dimension(3,node_num2) :: node_xyz2
integer(4),dimension(10,element_num) :: element_node2
integer(4),intent(in) :: nodesurf_num2
integer(4),dimension(3,nodesurf_num2),intent(in) :: surf_node
character(50) :: tet10file
integer(4) :: i,j,k,ishift
! out oceanfile ### node
open(1,file=tet10file)
write(1,'(a11)') "$MeshFormat"
write(1,'(a7)') "2.2 0 8"
write(1,'(a14)') "$EndMeshFormat"
write(1,'(a6)') "$Nodes"
write(1,*) node_num2
do i=1,node_num2
write(1,*) i, (node_xyz2(j,i),j=1,3)
end do
write(1,'(a)') "$EndNodes"
! out oceanfile ### tetrahedron elements
write(1,'(a9)') "$Elements"
ishift=3*nodesurf_num2
write(1,*) ishift + element_num
k=0
do j=1,3
do i=1,nodesurf_num2
k=k+1
write(1,'(i7,a,2i7)') k,"  15   2   0 " ,  j, surf_node(j,i) ! Surface node #
end do
end do
do i=1,element_num
write(1,*) i+ishift, "  11   2   0   1 ",(element_node2(j,i),j=1,10) ! TET10 elements
end do
write(1,'(a)') "$EndElements"
close(1)
write(*,'(a)') "### WRITETET10 END ###"
return
end subroutine writetet10
!#####
subroutine ch_cap ( c )
!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
implicit none

character c
integer ( kind = 4 ) itemp

itemp = ichar ( c )

if ( 97 <= itemp .and. itemp <= 122 ) then
c = char ( itemp - 32 )
end if

return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
implicit none

logical ch_eqi
character c1
character c1_cap
character c2
character c2_cap

c1_cap = c1
c2_cap = c2

call ch_cap ( c1_cap )
call ch_cap ( c2_cap )

if ( c1_cap == c2_cap ) then
ch_eqi = .true.
else
ch_eqi = .false.
end if

return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding value.  If C was
!    'illegal', then DIGIT is -1.
!
implicit none

character c
integer ( kind = 4 ) digit

if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

digit = ichar ( c ) - 48

else if ( c == ' ' ) then

digit = 0

else

digit = -1

end if

return
end
subroutine file_column_count ( input_filename, column_num )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Most lines of the file is presumed to consist of COLUMN_NUM words,
!    separated by spaces.  There may also be some blank lines, and some
!    comment lines,
!    which have a "#" in column 1.
!
!    The routine tries to find the first non-comment non-blank line and
!    counts the number of words in that line.
!
!    If all lines are blanks or comments, it goes back and tries to analyze
!    a comment line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the file.
!
!    Output, integer ( kind = 4 ) COLUMN_NUM, the number of columns in the file.
!
implicit none

integer ( kind = 4 ) column_num
logical got_one
character ( len = * ) input_filename
integer ( kind = 4 ) input_unit
integer ( kind = 4 ) ios
character ( len = 255 ) line
!
!  Open the file.
!
call get_unit ( input_unit )

open ( unit = input_unit, file = input_filename, status = 'old', &
form = 'formatted', access = 'sequential', iostat = ios )

if ( ios /= 0 ) then
column_num = -1
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
write ( *, '(a)' ) '  Could not open the file:'
write ( *, '(a)' ) '    ' // trim ( input_filename )
return
end if
!
!  Read one line, but skip blank lines and comment lines.
!
got_one = .false.

do

read ( input_unit, '(a)', iostat = ios ) line

if ( ios /= 0 ) then
exit
end if

if ( len_trim ( line ) == 0 ) then
cycle
end if

if ( line(1:1) == '#' ) then
cycle
end if

got_one = .true.
exit

end do

if ( .not. got_one ) then

rewind ( input_unit )

do

read ( input_unit, '(a)', iostat = ios ) line

if ( ios /= 0 ) then
exit
end if

if ( len_trim ( line ) == 0 ) then
cycle
end if

got_one = .true.
exit

end do

end if

close ( unit = input_unit )

if ( .not. got_one ) then
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
write ( *, '(a)' ) '  The file does not seem to contain any data.'
column_num = -1
return
end if

call s_word_count ( line, column_num )

return
end
subroutine file_row_count ( input_filename, row_num )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of row records in a file.
!
!  Discussion:
!
!    It does not count lines that are blank, or that begin with a
!    comment symbol '#'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) ROW_NUM, the number of rows found.
!
implicit none

integer ( kind = 4 ) bad_num
integer ( kind = 4 ) comment_num
integer ( kind = 4 ) ierror
character ( len = * ) input_filename
integer ( kind = 4 ) input_unit
integer ( kind = 4 ) ios
character ( len = 255 ) line
integer ( kind = 4 ) record_num
integer ( kind = 4 ) row_num

call get_unit ( input_unit )

open ( unit = input_unit, file = input_filename, status = 'old', &
iostat = ios )

if ( ios /= 0 ) then
row_num = -1;
ierror = 1
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
write ( *, '(a)' ) '  Could not open the input file: ' // &
trim ( input_filename )
stop
end if

comment_num = 0
row_num = 0
record_num = 0
bad_num = 0

do

read ( input_unit, '(a)', iostat = ios ) line

if ( ios /= 0 ) then
ierror = record_num
exit
end if

record_num = record_num + 1

if ( line(1:1) == '#' ) then
comment_num = comment_num + 1
cycle
end if

if ( len_trim ( line ) == 0 ) then
comment_num = comment_num + 1
cycle
end if

row_num = row_num + 1

end do

close ( unit = input_unit )

return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
implicit none

integer ( kind = 4 ) i
integer ( kind = 4 ) ios
integer ( kind = 4 ) iunit
logical lopen

iunit = 0

do i = 1, 99

if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

inquire ( unit = i, opened = lopen, iostat = ios )

if ( ios == 0 ) then
if ( .not. lopen ) then
iunit = i
return
end if
end if

end if

end do

return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
implicit none

integer ( kind = 4 ) i
integer ( kind = 4 ) j
integer ( kind = 4 ) k

k = i
i = j
j = k

return
end
subroutine i4col_compare ( m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! I4COL_COMPARE compares columns I and J of an I4COL.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of N columns of vectors
!    of length M.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column J < column I.
!
implicit none

integer ( kind = 4 ) m
integer ( kind = 4 ) n

integer ( kind = 4 ) a(m,n)
integer ( kind = 4 ) i
integer ( kind = 4 ) isgn
integer ( kind = 4 ) j
integer ( kind = 4 ) k
!
!  Check.
!
if ( i < 1 .or. n < i ) then
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
write ( *, '(a)' ) '  Column index I is out of bounds.'
stop
end if

if ( j < 1 .or. n < j ) then
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
write ( *, '(a)' ) '  Column index J is out of bounds.'
stop
end if

isgn = 0

if ( i == j ) then
return
end if

k = 1

do while ( k <= m )

if ( a(k,i) < a(k,j) ) then
isgn = -1
return
else if ( a(k,j) < a(k,i) ) then
isgn = +1
return
end if

k = k + 1

end do

return
end
subroutine i4col_sort_a ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT_A ascending sorts an I4COL.
!
!  Discussion:
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in ascending
!    lexicographic order.
!
implicit none

integer ( kind = 4 ) m
integer ( kind = 4 ) n

integer ( kind = 4 ) a(m,n)
integer ( kind = 4 ) i
integer ( kind = 4 ) indx
integer ( kind = 4 ) isgn
integer ( kind = 4 ) j

if ( m <= 0 ) then
return
end if

if ( n <= 1 ) then
return
end if
!
!  Initialize.
!
i = 0
indx = 0
isgn = 0
j = 0
!
!  Call the external heap sorter.
!
do

call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
if ( 0 < indx ) then

call i4col_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
else if ( indx < 0 ) then

call i4col_compare ( m, n, a, i, j, isgn )

else if ( indx == 0 ) then

exit

end if

end do

return
end
subroutine i4col_sort2_a ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT2_A ascending sorts the elements of each column of an I4COL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A, and the length
!    of a vector of data.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of N columns of M vectors.
!    On output, the elements of each column of A have been sorted in descending
!    order.
!
implicit none

integer ( kind = 4 ) m
integer ( kind = 4 ) n

integer ( kind = 4 ) a(m,n)
integer ( kind = 4 ) col
integer ( kind = 4 ) i
integer ( kind = 4 ) indx
integer ( kind = 4 ) isgn
integer ( kind = 4 ) j

if ( m <= 1 ) then
return
end if

if ( n <= 0 ) then
return
end if
!
!  Initialize.
!
do col = 1, n

i = 0
indx = 0
isgn = 0
j = 0
!
!  Call the external heap sorter.
!
do

call sort_heap_external ( m, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
if ( 0 < indx ) then

call i4_swap ( a(i,col), a(j,col) )
!
!  Compare the I and J objects.
!
else if ( indx < 0 ) then

if ( a(j,col) < a(i,col) ) then
isgn = +1
else
isgn = -1
end if

else if ( indx == 0 ) then

exit

end if

end do

end do

return
end
subroutine i4col_sorted_unique_count ( m, n, a, unique_num )

!*****************************************************************************80
!
!! I4COL_SORTED_UNIQUE_COUNT counts unique elements in an I4COL.
!
!  Discussion:
!
!    The columns of the array may be ascending or descending sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), a sorted array, containing
!    N columns of data.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique columns.
!
implicit none

integer ( kind = 4 ) m
integer ( kind = 4 ) n

integer ( kind = 4 ) a(m,n)
integer ( kind = 4 ) j1
integer ( kind = 4 ) j2
integer ( kind = 4 ) unique_num

if ( n <= 0 ) then
unique_num = 0
return
end if

unique_num = 1
j1 = 1

do j2 = 2, n

if ( any ( a(1:m,j1) /= a(1:m,j2) ) ) then
unique_num = unique_num + 1
j1 = j2
end if

end do

return
end
subroutine i4col_swap ( m, n, a, i, j )

!*****************************************************************************80
!
!! I4COL_SWAP swaps columns I and J of an I4COL.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      A = (
!        1  4  3  2
!        5  8  7  6
!        9 12 11 10 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of N columns of length M.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be swapped.
!
implicit none

integer ( kind = 4 ) m
integer ( kind = 4 ) n

integer ( kind = 4 ) a(m,n)
integer ( kind = 4 ) col(m)
integer ( kind = 4 ) i
integer ( kind = 4 ) j

if ( i < 1 .or. n < i .or. j < 1 .or. n < j ) then

write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'I4COL_SWAP - Fatal error!'
write ( *, '(a)' ) '  I or J is out of bounds.'
write ( *, '(a,i8)' ) '  I =    ', i
write ( *, '(a,i8)' ) '  J =    ', j
write ( *, '(a,i8)' ) '  N =    ', n
stop

end if

if ( i == j ) then
return
end if

col(1:m) = a(1:m,i)
a(1:m,i) = a(1:m,j)
a(1:m,j) = col(1:m)

return
end
subroutine i4i4_sort_a ( i1, i2, j1, j2 )

!*****************************************************************************80
!
!! I4I4_SORT_A ascending sorts a pair of I4's.
!
!  Discussion:
!
!    The program allows the reasonable call:
!
!      call i4i4_sort_a ( i1, i2, i1, i2 )
!
!    and this will return the reasonable result.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I1, I2, the values to sort.
!
!    Output, integer ( kind = 4 ) J1, J2, the sorted values.
!
implicit none

integer ( kind = 4 ) i1
integer ( kind = 4 ) i2
integer ( kind = 4 ) j1
integer ( kind = 4 ) j2
integer ( kind = 4 ) k1
integer ( kind = 4 ) k2
!
!  Copy arguments, so that the user can make "reasonable" calls like:
!
!    call i4i4_sort_a ( i1, i2, i1, i2 )
!
k1 = i1
k2 = i2

j1 = min ( k1, k2 )
j2 = max ( k1, k2 )

return
end
function s_index_last ( s, sub )

!*****************************************************************************80
!
!! S_INDEX_LAST finds the LAST occurrence of a given substring.
!
!  Discussion:
!
!    It returns the location in the string at which the substring SUB is
!    first found, or 0 if the substring does not occur at all.
!
!    The routine is also trailing blank insensitive.  This is very
!    important for those cases where you have stored information in
!    larger variables.  If S is of length 80, and SUB is of
!    length 80, then if S = 'FRED' and SUB = 'RED', a match would
!    not be reported by the standard FORTRAN INDEX, because it treats
!    both variables as being 80 characters long!  This routine assumes that
!    trailing blanks represent garbage!
!
!    This means that this routine cannot be used to find, say, the last
!    occurrence of a substring 'A ', since it assumes the blank space
!    was not specified by the user, but is, rather, padding by the
!    system.  However, as a special case, this routine can properly handle
!    the case where either S or SUB is all blanks.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be searched.
!
!    Input, character ( len = * ) SUB, the substring to search for.
!
!    Output, integer ( kind = 4 ) S_INDEX_LAST.  0 if SUB does not occur in
!    the string.  Otherwise S_INDEX_LAST = I, where S(I:I+LENS-1) = SUB,
!    where LENS is the length of SUB, and is the last place
!    this happens.
!
implicit none

integer ( kind = 4 ) i
integer ( kind = 4 ) j
integer ( kind = 4 ) llen1
integer ( kind = 4 ) llen2
character ( len = * ) s
integer ( kind = 4 ) s_index_last
character ( len = * ) sub

s_index_last = 0

llen1 = len_trim ( s )
llen2 = len_trim ( sub )
!
!  In case S or SUB is blanks, use LEN
!
if ( llen1 == 0 ) then
llen1 = len ( s )
end if

if ( llen2 == 0 ) then
llen2 = len ( sub )
end if

if ( llen1 < llen2 ) then
return
end if

do j = 1, llen1+1-llen2

i = llen1 + 2 - llen2 - j

if ( s(i:i+llen2-1) == sub ) then
s_index_last = i
return
end if

end do

return
end
subroutine s_to_i4 ( s, ival, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters of S used.
!
implicit none

character c
integer ( kind = 4 ) i
integer ( kind = 4 ) ierror
integer ( kind = 4 ) isgn
integer ( kind = 4 ) istate
integer ( kind = 4 ) ival
integer ( kind = 4 ) length
character ( len = * ) s

ierror = 0
istate = 0
isgn = 1
ival = 0

do i = 1, len_trim ( s )

c = s(i:i)
!
!  Haven't read anything.
!
if ( istate == 0 ) then

if ( c == ' ' ) then

else if ( c == '-' ) then
istate = 1
isgn = -1
else if ( c == '+' ) then
istate = 1
isgn = + 1
else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
istate = 2
ival = ichar ( c ) - ichar ( '0' )
else
ierror = 1
return
end if
!
!  Have read the sign, expecting digits.
!
else if ( istate == 1 ) then

if ( c == ' ' ) then

else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
istate = 2
ival = ichar ( c ) - ichar ( '0' )
else
ierror = 1
return
end if
!
!  Have read at least one digit, expecting more.
!
else if ( istate == 2 ) then

if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
ival = 10 * ival + ichar ( c ) - ichar ( '0' )
else
ival = isgn * ival
length = i - 1
return
end if

end if

end do
!
!  If we read all the characters in the string, see if we're OK.
!
if ( istate == 2 ) then
ival = isgn * ival
length = len_trim ( s )
else
ierror = 1
length = 0
end if

return
end
subroutine s_to_i4vec ( s, n, ivec, ierror )

!*****************************************************************************80
!
!! S_TO_I4VEC reads an I4VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer ( kind = 4 ) N, the number of values expected.
!
!    Output, integer ( kind = 4 ) IVEC(N), the values read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
implicit none

integer ( kind = 4 ) n

integer ( kind = 4 ) i
integer ( kind = 4 ) ierror
integer ( kind = 4 ) ilo
integer ( kind = 4 ) ivec(n)
integer ( kind = 4 ) length
character ( len = * ) s

i = 0

ilo = 1

do while ( i < n )

i = i + 1

call s_to_i4 ( s(ilo:), ivec(i), ierror, length )

if ( ierror /= 0 ) then
ierror = -i
exit
end if

ilo = ilo + length

end do

return
end
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8VEC from a string.
!
!  Discussion:
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    S                 DVAL
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) DVAL, the value read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
implicit none

logical ch_eqi
character c
real ( kind = 8 ) dval
integer ( kind = 4 ) ierror
integer ( kind = 4 ) ihave
integer ( kind = 4 ) isgn
integer ( kind = 4 ) iterm
integer ( kind = 4 ) jbot
integer ( kind = 4 ) jsgn
integer ( kind = 4 ) jtop
integer ( kind = 4 ) length
integer ( kind = 4 ) nchar
integer ( kind = 4 ) ndig
real ( kind = 8 ) rbot
real ( kind = 8 ) rexp
real ( kind = 8 ) rtop
character ( len = * ) s

nchar = len_trim ( s )

ierror = 0
dval = 0.0D+00
length = -1
isgn = 1
rtop = 0
rbot = 1
jsgn = 1
jtop = 0
jbot = 1
ihave = 1
iterm = 0

do

length = length + 1

if ( nchar < length+1 ) then
exit
end if

c = s(length+1:length+1)
!
!  Blank character.
!
if ( c == ' ' ) then

if ( ihave == 2 ) then

else if ( ihave == 6 .or. ihave == 7 ) then
iterm = 1
else if ( 1 < ihave ) then
ihave = 11
end if
!
!  Comma.
!
else if ( c == ',' .or. c == ';' ) then

if ( ihave /= 1 ) then
iterm = 1
ihave = 12
length = length + 1
end if
!
!  Minus sign.
!
else if ( c == '-' ) then

if ( ihave == 1 ) then
ihave = 2
isgn = -1
else if ( ihave == 6 ) then
ihave = 7
jsgn = -1
else
iterm = 1
end if
!
!  Plus sign.
!
else if ( c == '+' ) then

if ( ihave == 1 ) then
ihave = 2
else if ( ihave == 6 ) then
ihave = 7
else
iterm = 1
end if
!
!  Decimal point.
!
else if ( c == '.' ) then

if ( ihave < 4 ) then
ihave = 4
else if ( 6 <= ihave .and. ihave <= 8 ) then
ihave = 9
else
iterm = 1
end if
!
!  Scientific notation exponent marker.
!
else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

if ( ihave < 6 ) then
ihave = 6
else
iterm = 1
end if
!
!  Digit.
!
else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

if ( ihave <= 2 ) then
ihave = 3
else if ( ihave == 4 ) then
ihave = 5
else if ( ihave == 6 .or. ihave == 7 ) then
ihave = 8
else if ( ihave == 9 ) then
ihave = 10
end if

call ch_to_digit ( c, ndig )

if ( ihave == 3 ) then
rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
else if ( ihave == 5 ) then
rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
rbot = 10.0D+00 * rbot
else if ( ihave == 8 ) then
jtop = 10 * jtop + ndig
else if ( ihave == 10 ) then
jtop = 10 * jtop + ndig
jbot = 10 * jbot
end if
!
!  Anything else is regarded as a terminator.
!
else
iterm = 1
end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
if ( iterm == 1 ) then
exit
end if

end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to NCHAR.
!
if ( iterm /= 1 .and. length+1 == nchar ) then
length = nchar
end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
ierror = ihave
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
write ( *, '(a)' ) '  Illegal or nonnumeric input:'
write ( *, '(a)' ) '    ' // trim ( s )
return
end if
!
!  Number seems OK.  Form it.
!
if ( jtop == 0 ) then
rexp = 1.0D+00
else
if ( jbot == 1 ) then
rexp = 10.0D+00 ** ( jsgn * jtop )
else
rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = 8 ) &
/ real ( jbot, kind = 8 ) )
end if
end if

dval = real ( isgn, kind = 8 ) * rexp * rtop / rbot

return
end
subroutine s_to_r8vec ( s, n, rvec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer ( kind = 4 ) N, the number of values expected.
!
!    Output, real ( kind = 8 ) RVEC(N), the values read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
implicit none

integer ( kind = 4 ) n

integer ( kind = 4 ) i
integer ( kind = 4 ) ierror
integer ( kind = 4 ) ilo
integer ( kind = 4 ) lchar
real ( kind = 8 ) rvec(n)
character ( len = * ) s

i = 0

ilo = 1

do while ( i < n )

i = i + 1

call s_to_r8 ( s(ilo:), rvec(i), ierror, lchar )

if ( ierror /= 0 ) then
ierror = -i
exit
end if

ilo = ilo + lchar

end do

return
end
subroutine s_word_count ( s, nword )

!*****************************************************************************80
!
!! S_WORD_COUNT counts the number of "words" in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be examined.
!
!    Output, integer ( kind = 4 ) NWORD, the number of "words" in the string.
!    Words are presumed to be separated by one or more blanks.
!
implicit none

logical blank
integer ( kind = 4 ) i
integer ( kind = 4 ) lens
integer ( kind = 4 ) nword
character ( len = * ) s

nword = 0
lens = len ( s )

if ( lens <= 0 ) then
return
end if

blank = .true.

do i = 1, lens

if ( s(i:i) == ' ' ) then
blank = .true.
else if ( blank ) then
nword = nword + 1
blank = .false.
end if

end do

return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, real ( kind = 8 )s, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements
!    I and J.
!    (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
implicit none

integer ( kind = 4 ) i
integer ( kind = 4 ), save :: i_save = 0
integer ( kind = 4 ) indx
integer ( kind = 4 ) isgn
integer ( kind = 4 ) j
integer ( kind = 4 ), save :: j_save = 0
integer ( kind = 4 ), save :: k = 0
integer ( kind = 4 ), save :: k1 = 0
integer ( kind = 4 ) n
integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
if ( indx == 0 ) then

i_save = 0
j_save = 0
k = n / 2
k1 = k
n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
else if ( indx < 0 ) then

if ( indx == -2 ) then

if ( isgn < 0 ) then
i_save = i_save + 1
end if

j_save = k1
k1 = i_save
indx = -1
i = i_save
j = j_save
return

end if

if ( 0 < isgn ) then
indx = 2
i = i_save
j = j_save
return
end if

if ( k <= 1 ) then

if ( n1 == 1 ) then
i_save = 0
j_save = 0
indx = 0
else
i_save = n1
n1 = n1 - 1
j_save = 1
indx = 1
end if

i = i_save
j = j_save
return

end if

k = k - 1
k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
else if ( indx == 1 ) then

k1 = k

end if

do

i_save = 2 * k1

if ( i_save == n1 ) then
j_save = k1
k1 = i_save
indx = -1
i = i_save
j = j_save
return
else if ( i_save <= n1 ) then
j_save = i_save + 1
indx = -2
i = i_save
j = j_save
return
end if

if ( k <= 1 ) then
exit
end if

k = k - 1
k1 = k

end do

if ( n1 == 1 ) then
i_save = 0
j_save = 0
indx = 0
i = i_save
j = j_save
else
i_save = n1
n1 = n1 - 1
j_save = 1
indx = 1
i = i_save
j = j_save
end if

return
end
subroutine tet_mesh_order4_to_order10_compute ( tetra_num, tetra_node1, &
node_num1, node_xyz1, edge_data, tetra_node2, node_num2, node_xyz2 )

!*****************************************************************************80
!
!! TET_MESH_ORDER4_TO_ORDER10_COMPUTE: quadratic tet mesh from a linear one.
!
!  Discussion:
!
!    A quadratic (10 node) tet mesh can be derived from a linear
!    (4 node) tet mesh by interpolating nodes at the midpoint of
!    every edge of the mesh.
!
!    The mesh is described indirectly, as the sum of individual
!    tetrahedrons.  A single physical edge may be a logical edge of
!    any number of tetrahedrons.  It is important, however, that a
!    new node be created exactly once for each edge, assigned an index,
!    and associated with every tetrahedron that shares this edge.
!
!    This routine handles that problem.
!
!    The primary amount of work occurs in sorting a list of 6 * TETRA_NUM
!    data items, one item for every edge of every tetrahedron.  Each
!    data item records, for a given tetrahedron edge, the global indices
!    of the two endpoints, the local indices of the two endpoints,
!    and the index of the tetrahedron.
!
!    Through careful sorting, it is possible to arrange this data in
!    a way that allows the proper generation of the interpolated nodes.
!
!    The node ordering for the quadratic tetrahedron is somewhat
!    arbitrary.  In the current scheme, the vertices are listed
!    first, followed by the 6 midside nodes.  Each midside node
!    may be identified by the two vertices that bracket it.  Thus,
!    the node ordering may be suggested by:
!
!      1  2  3  4 (1+2) (1+3) (1+4) (2+3) (2+4) (3+4)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TETRA_NUM, the number of tetrahedrons in the
!    linear mesh.
!
!    Input, integer ( kind = 4 ) TETRA_NODE1(4,TETRA_NUM), the indices of the
!    nodes in the linear mesh.
!
!    Input, integer ( kind = 4 ) NODE_NUM1, the number of nodes for the
!    linear mesh.
!
!    Input, real ( kind = 8 ) NODE_XYZ1(3,NODE_NUM1), the coordinates of
!    the nodes that make up the linear mesh.
!
!    Input, integer ( kind = 4 ) EDGE_DATA(5,6*TETRA_NUM), edge data.
!
!    Output, integer ( kind = 4 ) TETRA_NODE2(10,TETRA_NUM), the indices of
!    the nodes in the quadratic mesh.
!
!    Input, integer ( kind = 4 ) NODE_NUM2, the number of nodes for the
!    quadratic mesh.
!
!    Output, real ( kind = 8 ) NODE_XYZ2(3,NODE_NUM2), the coordinates of
!    the nodes that make up the quadratic mesh.
!
implicit none

integer ( kind = 4 ) node_num1
integer ( kind = 4 ) node_num2
integer ( kind = 4 ) tetra_num

integer ( kind = 4 ) edge
integer ( kind = 4 ) edge_data(5,6*tetra_num)
integer ( kind = 4 ) n1
integer ( kind = 4 ) n1_old
integer ( kind = 4 ) n2
integer ( kind = 4 ) n2_old
integer ( kind = 4 ) node
real ( kind = 8 ) node_xyz1(3,node_num1)
real ( kind = 8 ) node_xyz2(3,node_num2)
integer ( kind = 4 ) tetra
integer ( kind = 4 ) tetra_node1(4,tetra_num)
integer ( kind = 4 ) tetra_node2(10,tetra_num)
integer ( kind = 4 ) v
integer ( kind = 4 ) v1
integer ( kind = 4 ) v2
!
!  Generate the index and coordinates of the new midside nodes,
!  and update the tetradehron-node data.
!
node_xyz2(1:3,1:node_num1) = node_xyz1(1:3,1:node_num1)

tetra_node2(1:4,1:tetra_num) = tetra_node1(1:4,1:tetra_num)

node = node_num1

n1_old = -1
n2_old = -1

do edge = 1, 6 * tetra_num
!
!  Read the data defining the edge.
!
n1 = edge_data(1,edge)
n2 = edge_data(2,edge)
!
!  If this edge is new, create the coordinates and index.
!
if ( n1 /= n1_old .or. n2 /= n2_old ) then

node = node + 1

if ( node_num2 < node ) then
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'TET_MESH_ORDER4_TO_ORDER10_COMPUTE - Fatal error!'
write ( *, '(a)' ) '  Node index exceeds NODE_NUM2.'
stop
end if

node_xyz2(1:3,node) = &
( node_xyz2(1:3,n1) + node_xyz2(1:3,n2) ) / 2.0D+00

n1_old = n1
n2_old = n2

end if
!
!  Assign the node to the tetrahedron.
!
v1 = edge_data(3,edge)
v2 = edge_data(4,edge)
!
!  Here is where the local ordering of the nodes is effected:
!
if ( v1 == 1 .and. v2 == 2 ) then
v = 5
else if ( v1 == 1 .and. v2 == 3 ) then
v = 6
else if ( v1 == 1 .and. v2 == 4 ) then
v = 7
else if ( v1 == 2 .and. v2 == 3 ) then
v = 8
else if ( v1 == 2 .and. v2 == 4 ) then
v = 9
else if ( v1 == 3 .and. v2 == 4 ) then
v = 10
end if

tetra = edge_data(5,edge)

tetra_node2(v,tetra) = node

end do

return
end
subroutine tet_mesh_order4_to_order10_size ( tetra_num, tetra_node1, &
node_num1, edge_data, node_num2 )
!*****************************************************************************80
!
!! TET_MESH_ORDER4_TO_ORDER10_SIZE sizes a quadratic tet mesh from a linear one.
!
!  Discussion:
!
!    A quadratic (10 node) tet mesh can be derived from a linear
!    (4 node) tet mesh by interpolating nodes at the midpoint of
!    every edge of the mesh.
!
!    The mesh is described indirectly, as the sum of individual
!    tetrahedrons.  A single physical edge may be a logical edge of
!    any number of tetrahedrons.  It is important, however, that a
!    new node be created exactly once for each edge, assigned an index,
!    and associated with every tetrahedron that shares this edge.
!
!    This routine handles that problem.
!
!    The primary amount of work occurs in sorting a list of 6 * TETRA_NUM
!    data items, one item for every edge of every tetrahedron.  Each
!    data item records, for a given tetrahedron edge, the global indices
!    of the two endpoints, the local indices of the two endpoints,
!    and the index of the tetrahedron.
!
!    Through careful sorting, it is possible to arrange this data in
!    a way that allows the proper generation of the interpolated nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TETRA_NUM, the number of tetrahedrons in the
!    linear mesh.
!
!    Input, integer ( kind = 4 ) TETRA_NODE1(4,TETRA_NUM), the indices of the
!    nodes in the linear mesh.
!
!    Input, integer ( kind = 4 ) NODE_NUM1, the number of nodes for the
!    linear mesh.
!
!    Output, integer ( kind = 4 ) EDGE_DATA(5,6*TETRA_NUM), edge data.
!
!    Output, integer ( kind = 4 ) NODE_NUM2, the number of nodes for the
!    quadratic mesh.
!
implicit none

integer ( kind = 4 ) node_num1
integer ( kind = 4 ) tetra_num

integer ( kind = 4 ) a
integer ( kind = 4 ) b
integer ( kind = 4 ) edge
integer ( kind = 4 ) edge_data(5,6*tetra_num)
integer ( kind = 4 ) i
integer ( kind = 4 ) j
integer ( kind = 4 ) k
integer ( kind = 4 ) l
integer ( kind = 4 ) n1
integer ( kind = 4 ) n1_old
integer ( kind = 4 ) n2
integer ( kind = 4 ) n2_old
integer ( kind = 4 ) node_num2
integer ( kind = 4 ) tetra
integer ( kind = 4 ) tetra_node1(4,tetra_num)
!
!  Step 1.
!  From the list of nodes for tetrahedron T, of the form: (I,J,K,L)
!  construct the six edge relations:
!
!    (I,J,1,2,T)
!    (I,K,1,3,T)
!    (I,L,1,4,T)
!    (J,K,2,3,T)
!    (J,L,2,4,T)
!    (K,L,3,4,T)
!
!  In order to make matching easier, we reorder each pair of nodes
!  into ascending order.
!
do tetra = 1, tetra_num

i = tetra_node1(1,tetra)
j = tetra_node1(2,tetra)
k = tetra_node1(3,tetra)
l = tetra_node1(4,tetra)

call i4i4_sort_a ( i, j, a, b )

edge_data(1:5,6*(tetra-1)+1) = (/ a, b, 1, 2, tetra /)

call i4i4_sort_a ( i, k, a, b )

edge_data(1:5,6*(tetra-1)+2) = (/ a, b, 1, 3, tetra /)

call i4i4_sort_a ( i, l, a, b )

edge_data(1:5,6*(tetra-1)+3) = (/ a, b, 1, 4, tetra /)

call i4i4_sort_a ( j, k, a, b )

edge_data(1:5,6*(tetra-1)+4) = (/ a, b, 2, 3, tetra /)

call i4i4_sort_a ( j, l, a, b )

edge_data(1:5,6*(tetra-1)+5) = (/ a, b, 2, 4, tetra /)

call i4i4_sort_a ( k, l, a, b )

edge_data(1:5,6*(tetra-1)+6) = (/ a, b, 3, 4, tetra /)

end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!  We only intend to sort on rows 1:2; the routine we call here
!  sorts on the full column but that won't hurt us.
!
!  What we need is to find all cases where tetrahedrons share an edge.
!  By sorting the columns of the EDGE_DATA array, we will put shared edges
!  next to each other.
!
call i4col_sort_a ( 5, 6*tetra_num, edge_data )
!
!  Step 3. All the tetrahedrons which share an edge show up as consecutive
!  columns with identical first two entries.  Figure out how many new
!  nodes there are, and allocate space for their coordinates.
!
node_num2 = node_num1

n1_old = -1
n2_old = -1

do edge = 1, 6 * tetra_num
n1 = edge_data(1,edge)
n2 = edge_data(2,edge)
if ( n1 /= n1_old .or. n2 /= n2_old ) then
node_num2 = node_num2 + 1
n1_old = n1
n2_old = n2
end if
end do

return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
implicit none

character ( len = 8 ) ampm
integer ( kind = 4 ) d
integer ( kind = 4 ) h
integer ( kind = 4 ) m
integer ( kind = 4 ) mm
character ( len = 9  ), parameter, dimension(12) :: month = (/ &
'January  ', 'February ', 'March    ', 'April    ', &
'May      ', 'June     ', 'July     ', 'August   ', &
'September', 'October  ', 'November ', 'December ' /)
integer ( kind = 4 ) n
integer ( kind = 4 ) s
integer ( kind = 4 ) values(8)
integer ( kind = 4 ) y

call date_and_time ( values = values )

y = values(1)
m = values(2)
d = values(3)
h = values(5)
n = values(6)
s = values(7)
mm = values(8)

if ( h < 12 ) then
ampm = 'AM'
else if ( h == 12 ) then
if ( n == 0 .and. s == 0 ) then
ampm = 'Noon'
else
ampm = 'PM'
end if
else
h = h - 12
if ( h < 12 ) then
ampm = 'PM'
else if ( h == 12 ) then
if ( n == 0 .and. s == 0 ) then
ampm = 'Midnight'
else
ampm = 'AM'
end if
end if
end if

write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

return
end