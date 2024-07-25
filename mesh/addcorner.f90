!!## May 22, 2014, zlabel is introduced to identify zero-depth nodes.
module addcorner11m
use coastline_data
implicit none
contains
!########################################################### addcorner11
subroutine addcorner11(h_poly,i_poly,g_grd,g_bound,iflg,ilflg,xcorner,ycorner)
implicit none
type(grid_data),       intent(in)       :: g_grd
type(poly_data),       intent(in)       :: h_poly
type(poly_data),       intent(out)      :: i_poly
type(bound_data),      intent(inout)    :: g_bound
integer(4)                              :: ncmax, lpmax ! 2018.08.28
integer(4),allocatable,dimension(:)     :: zlabel,ibelong
integer(4),allocatable,dimension(:)     :: lpoly
integer(4),allocatable,dimension(:,:,:) :: ind2
real(8)   ,allocatable,dimension(:,:)   :: loc
integer(4)                              :: ln(2)
integer(4)                              :: ilflg(3)
! ilflg(i) : polygon # within nclose land polygons
!             which i-th corner nodes should be included
!loc(i,1)  : the starting area of i-th unclosed polygon,
!loc(i,2)  : end area of i-th unclosed polygon
! for the area, right boundary is area 0.5, bottom is 1.5, left is 2.5, top is 3.5
real(8),    allocatable,dimension(:,:)  :: xpoly,ypoly ! 4 for corner nodes
real(8),    allocatable,dimension(:)    :: xpoly2,ypoly2
! ncoast includes real # of coast points + 4 for corners
real(8),                dimension(3)    :: xcorner,ycorner
real(8),    allocatable,dimension(:)    :: x,y
integer(4)                              :: neast,node,nclose,lpoly0,iflg(3)
integer(4)                              :: i,j,ii,jj,lpoly2
!!!NOTE that this program assumes top right corner is in the ocean!!!
! ind2(i,j,1) is grid node # of original gebco file which i-th plygon's starting node(j=1), or end node (j=2)
! ind2(i,j,2) = 1 indicates ind(i,j,1) gebco node is on the left side of i-th coastline node
!  ind2(i,j,2)=2 -> gebco node ind(i,1) is on the topside
!          ind2(i,j,1)     <-   ind2(i,j,2)=1
!            ^
!         ind2(i,j,2)=2
!### zlabel(i) = 0 : depth at i-node is zero, 1 : depth at i-th node is not zero
!### iflg(i) is the polygon # after which the i-th corner node will be added #####
! 1st corner : bottom right, 2nd corner : bottom left, 3rd corner : top left

!#[0]## set input
 lpmax  = h_poly%lpmax    ! 2018.08.28
 node   = g_grd%node
 neast  = g_grd%neast
 ncmax  = g_bound%ncmax   ! 2018.08.28
 nclose = h_poly%nclose   ! # of unclosed polygons
 allocate(x(node),y(node))! 2018.08.28
 allocate(zlabel(ncmax),ibelong(ncmax)) ! 2018.08.28
 allocate(xpoly(lpmax,ncmax), ypoly(lpmax,ncmax)) ! 2018.08.28
 allocate(xpoly2(ncmax),      ypoly2(ncmax)  )    ! 2018.08.28
 allocate(lpoly(lpmax))   ! 2018.08.28
 allocate(ind2(nclose,2,2)) ! 2018.08.28
 allocate(loc(nclose,2)   ) ! 2018.08.28
 x(:)   = g_grd%xyz(2,:)  ! northward [km] 2018.08.28
 y(:)   = g_grd%xyz(1,:)  ! eastward [km]  2018.08.28
 !
 lpoly0 = h_poly%lpoly0
 lpoly  = h_poly%lpoly
 ypoly  = h_poly%xypoly(1,:,:) ! east
 xpoly  = h_poly%xypoly(2,:,:) ! north
 ind2   = h_poly%ind2
!
 if ( ncmax .ne. h_poly%ncmax ) then
  write(*,*) "GEGEGE! ncmax=",ncmax,"h_poly%ncmax=",h_poly%ncmax
  stop
 end if

iflg(1:3)=999 ! default
ii=0
!# [0] ### set loc(nclose,2)
loc(:,:)=0.d0
do i=1,nclose   ! nclose is new nclose after removing the small islands, comment on Oct. 24, 2015
  ii=ii+lpoly(i)    ! ii is accumulated number of unclosed polygon points, comment on Oct. 24, 2015
  ln(1)=ii-lpoly(i)+1; ln(2)=ii ! ln(1) is node # of first node of polygon, ln(2) is the end node #
  do j=1,2  ! start and end points of i-th polygons
    !write(*,*) "i,j=",i,j,"ln(j)=",ln(j),"ind(ln(j),1)=",ind(ln(j),1),"ind(ln(j),2)=",ind(ln(j),2)
    !write(*,*) "x,y=",x(ind(ln(j),1)),y(ind(ln(j),1))
    if ( mod(ind2(i,j,1),neast) .eq. 0 )                         loc(i,j)=0.5d0 ! right
    if ( mod(ind2(i,j,1),neast) .eq. 1 .and. ind2(i,j,2) .eq. 2) loc(i,j)=2.5d0 ! left
    if ( ind2(i,j,1) .ge. node-neast+1 )                         loc(i,j)=1.5d0 ! bottom
    if ( ind2(i,j,1) .le. neast-1 .and. ind2(i,j,2) .eq. 1 )     loc(i,j)=3.5d0 ! top
  end do
  if (loc(i,1) .eq. 0.d0 .or. loc(i,2) .eq. 0.d0) then
    write(*,*) "GEGEGE loc is not set!!"
    write(*,*) "i=",i,"loc(i,1)=",loc(i,1),"loc(i,2)=",loc(i,2)
    write(*,*) "ln(1)=",ln(1),"ln(2)=",ln(2)
    write(*,*) "ind2(i,1,1)=",ind2(i,1,1),"ind2(i,2,1)",ind2(i,2,1)
    write(*,*) "ind2(i,1,2)=",ind2(i,1,2),"ind2(i,2,2)",ind2(i,2,2)
    write(*,*) "lpoly(i)=",lpoly(i)
    write(*,'(a,2g15.7)')  ("x,y=",xpoly(i,j),ypoly(i,j),j=1,lpoly(i))
    stop
  end if
  write(*,*) "polygon #", i, "loc(i,1)=",loc(i,1)
  write(*,*) "polygon #", i, "loc(i,2)=",loc(i,2)
end do
! loc set end
!# [1-1] ####  before first boundary polygon
if ( ind2(1,1,1) .gt. node-neast .or. mod(ind2(1,1,1),neast) .eq. 1 .or. ind2(1,1,1) .lt. neast ) iflg(1)=0 ! 0 indicates befor the first polygon
if ( mod(ind2(1,1,1),neast) .eq. 1 .or. ind2(1,1,1) .lt. neast ) iflg(2)=0
if ( ind2(1,1,1) .lt. neast ) iflg(3)=0
!# [1-2] ####  between i-th and i+1th polygon
if ( nclose .ge. 2 ) then ! if nclose=1, [1-2] is not necessary.
  j=0
  do jj=1,nclose-1 !! investigate if corner is present between jj-th and jj+1-th boundary polygons
    j=j+lpoly(jj) ! j is the last node of jj-th unclosed polygon
!    if ( mod(ind2(jj,2,1),neast) .eq. 0 .and. ind2(jj+1,1,1) .gt. node-neast) iflg(1)=jj ! bottom right corner, deleted on Oct. 24, 2015
!    if ( ind2(jj,2,1) .gt. node-neast .and. mod(ind2(jj+1,1,1),neast) .eq. 1 )iflg(2)=jj ! bottom left corner
!    if ( mod(ind2(jj,2,1),neast) .eq. 1  .and. ind2(jj+1,1,1) .le. neast)     iflg(3)=jj ! top left corner
    if ( loc(jj,2) .lt. 1.d0 .and. loc(jj+1,1) .gt. 1.d0 ) iflg(1)=jj ! bottom right corner, inserted on Oct. 24, 2015
    if ( loc(jj,2) .lt. 2.d0 .and. loc(jj+1,1) .gt. 2.d0 ) iflg(2)=jj ! bottom left corner, inserted on Oct. 24, 2015
    if ( loc(jj,2) .lt. 3.d0 .and. loc(jj+1,1) .gt. 3.d0 ) iflg(3)=jj ! top left corner, inserted on Oct. 24, 2015
  end do
end if
!# [1-3] #### after the last polygon
j=j+lpoly(nclose)
if ( mod(ind2(nclose,2,1),neast) .eq. 0 ) iflg(1)=nclose ! bottom right boundary
if ( mod(ind2(nclose,2,1),neast) .eq. 0 .or. ind2(nclose,2,1) .gt. node-neast) iflg(2)=nclose ! bottom right boundary
if ( mod(ind2(nclose,2,1),neast) .eq. 0 .or. ind2(nclose,2,1) .gt. &
&    node-neast .or. mod(ind2(nclose,2,1),neast) .eq. 1) iflg(3)=nclose
!  iflg(1:3) is been set

!# [1-4] #### for i-th corner that iflg(i) =999, which should be included in the land polygon
do i=1,3
 if (iflg(i) .ne. 999) then
  ilflg(i)=0 ! ilflg(i)=0 indicates that i-th corner node is not included in any land polygons
 else
  do j=1,nclose
   if ( loc(j,1)  .lt. float(i) .and. float(i) .lt. loc(j,2) ) ilflg(i)=j ! select the land polygon which includes i-th corner
  end do
 end if
end do
write(*,*) "iflg(1:3)=",(iflg(i),i=1,3)
write(*,*) "ilflg(1:3)=",(ilflg(i),i=1,3)
!# [2-1] #### set corner nodes
xcorner(1)=x(node)         ; ycorner(1)=y(node)                       ! bottom right corner
xcorner(2)=x(node-neast+1) ; ycorner(2)=y(node-neast+1)! bottom left corner
xcorner(3)=x(1)            ; ycorner(3)=y(1)                               ! top left corner
!# [2-2] #### set top right corner node
xpoly2(1)=x(neast)
ypoly2(1)=y(neast)
lpoly2=1
zlabel(1)=1 !!! zlabel = 1 indicates depth here is not zero since this program assumes top right corner is in the ocean!!!
ibelong(1)=1 !!! ibelong(i) = 1 indicates boundary line starting from i-th node is not on coastlines

!# [3] #### integrate nclose boundary polygons including corner nodes ###
!# [3-1] #### before first boundary polygon
do j=1,3
 if (iflg(j) .eq. 0) then ! After iflg(j) polygon #, the j-th corner node is added. If iflg(j)=0, j-th corner is added before 1st polygon.
  xpoly2(lpoly2+1)=xcorner(j) ! lpoly2 is # of points for 1st ocean polygon
  ypoly2(lpoly2+1)=ycorner(j)
  lpoly2=lpoly2+1
  zlabel(lpoly2)=1 ! at added corner, the ocean depth is not zero
  ibelong(lpoly2)=1 ! On all lines starts from corner nodes, it is necessary to calculate ocean depths
 end if
end do

!# [3-2] #### after the first to the last polygon
do i=1,nclose !# nclose : # of polygon reaching calculation boundaries
 xpoly2(lpoly2+1:lpoly2+lpoly(i))=xpoly(i,1:lpoly(i))
 ypoly2(lpoly2+1:lpoly2+lpoly(i))=ypoly(i,1:lpoly(i))
 zlabel(lpoly2+1:lpoly2+lpoly(i))=0 ! All of the lpoly(i) nodes are on coastlines.
 ibelong(lpoly2+1:lpoly2+lpoly(i)-1)=0 ! the lines starting from the first and inner nodes are on coastlines
 ibelong(lpoly2+lpoly(i))=1 ! the last node is on coast line, but also on the boundaries.
 lpoly2=lpoly2+lpoly(i) ! Then on lines starting from last nodes, it is necessary to calculate ocean depths
do j=1,3
if (iflg(j) .eq. i) then
xpoly2(lpoly2+1)=xcorner(j)
ypoly2(lpoly2+1)=ycorner(j)
lpoly2=lpoly2+1
zlabel(lpoly2)=1 ! at added corner, the ocean depth is not zero
ibelong(lpoly2)=1 ! On lines starting from corner nodes, it is necessary to calculate depths
end if
end do
end do
do i=1,4
write(*,*) "xpoly2(),ypoly2()=",xpoly2(i),ypoly2(i)
end do

!# [4] #### deal with nclose +1 to lpoly0 th polygons
if ( lpoly0 - nclose + 1 .ge. 2 ) then
 xpoly(2:lpoly0-nclose+1,:)=xpoly(nclose+1:lpoly0,:)
 ypoly(2:lpoly0-nclose+1,:)=ypoly(nclose+1:lpoly0,:)
 lpoly(2:lpoly0-nclose+1)=lpoly(nclose+1:lpoly0)
 if ( lpoly0 - nclose +2 .le. h_poly%lpmax ) then
  xpoly(lpoly0-nclose+2:h_poly%lpmax,:)=0.d0
  ypoly(lpoly0-nclose+2:h_poly%lpmax,:)=0.d0
  lpoly(lpoly0-nclose+2:h_poly%lpmax)=0
 end if
end if
xpoly(1,:)=xpoly2(:) ! 1st ocean polygon
ypoly(1,:)=ypoly2(:)
lpoly(1)=lpoly2
lpoly0=lpoly0-nclose+1 ! nclose polygon has been combined to 1 polygon
zlabel(lpoly2+1:ncmax)=0 ! after the nclose polygon, all the nodes are on coastlines
ibelong(lpoly2+1:ncmax)=0 ! after the nclose polygon, all boundary lines are on coastlines
write(*,*) "lpoly0=",lpoly0,"lpoly(1)=",lpoly(1),"nclose=",nclose,"h_poly%lpoly0=",h_poly%lpoly0
!do i=1,lpoly(1)
!write(*,*) "xpoly(1,:),ypoly(1,:)=",xpoly(1,i),ypoly(1,i)
!end do

!#[5]## set output
call allocatepoly(i_poly,lpoly0,ncmax,0) ! nclose for i_poly = 0
i_poly%lpoly0          = lpoly0
i_poly%lpoly(1:lpoly0) = lpoly(1:lpoly0)
i_poly%xypoly(1,1:lpoly0,1:ncmax)   = ypoly(1:lpoly0,1:ncmax)
i_poly%xypoly(2,1:lpoly0,1:ncmax)   = xpoly(1:lpoly0,1:ncmax)
i_poly%nclose          = 0
g_bound%ibelong        = ibelong
g_bound%zlabel         = zlabel

write(*,*) "### addcorner11 end!###"
return
end subroutine addcorner11
end module addcorner11m
