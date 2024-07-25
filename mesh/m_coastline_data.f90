! Coded by T. Minami on 2016.09.03
module coastline_data
implicit none

type grid_data
 integer(4) :: node
 integer(4) :: neast,nsouth
 real(8),dimension(:,:),allocatable :: lonlatalt ! lon [deg], lat [deg], altitude [m]
 ! eastward[km], northward[km], upward[km]
  real(8),dimension(:,:),allocatable :: xyz
end type grid_data

! label    : polygon number, default = 0, meaning not belonging to any polygon
! lpoly0   : # of polygon, which is generated in makepolygon7
!           and modified in landpoligon10
! lpoly(i) : # of nodes consisting of i-th polygon
! nbound   : # of nodes on the calculation boundaries
! nclose   : # of unclosed polygon, which touch calculation boundaries

type coast_data
 integer(4) :: ncmax  ! max # of points on coastline
 integer(4) :: ncoast ! # of points on coastline
 integer(4) :: nbound ! # of boundary points
! cxy(1,i)  :  eastward, cxy(2,i) : northward
! ind(i,1)  :  grid node id of original gebco file on left or up side
! ind(i,2)  :  1 -> ind(1) indicates left node, 2 -> ind(2) is top side node
 real(8),   allocatable,dimension(:,:) :: cxy
 integer(4),allocatable,dimension(:,:) :: ind ! ind(ncmax,2)
 integer(4),allocatable,dimension(:,:) :: ind_r ! ind_r(ncmax,2) reordered ind
 integer(4),allocatable,dimension(:,:) :: label ! label(ncmax,2)
end type coast_data
! ind and ind_r are different in the order of node_num.
! former is in the same manner as topofile, but latter is in polygon manner

type poly_data
integer(4) :: lpmax  ! max # of polygons
integer(4) :: ncmax  ! max # of points for each polygon
integer(4) :: lpoly0 ! # of polygons
integer(4) :: nclose ! # of un-closed polygons, =nbound/2
integer(4),allocatable,dimension(:,:,:) :: ind2 ! ind2(nclose,2,2)
! ind2(i,1,:) : ind for start node of polygon
! ind2(i,2,:) : ind for end node of polygon
! ind2 is defined only for unclosed polygons because the start and end node is always
! assumed on the calculation boundary
integer(4),allocatable,dimension(:)     :: lpoly  ! # of points in each polygon
real(8),   allocatable,dimension(:,:,:) :: xypoly
! xypoly(1:2,1:ncmax,1:lpoly0) x,y component of
end type

type bound_data
 integer(4) :: ncmax
 integer(4),allocatable,dimension(:) :: zlabel
 integer(4),allocatable,dimension(:) :: ibelong
end type bound_data

contains

!######################################## ALLOCATEBOUND
subroutine allocatebound(g_bound, ncmax)
type(bound_data),intent(out) :: g_bound
integer(4),      intent(in)  :: ncmax

 g_bound%ncmax=ncmax
 allocate(g_bound%zlabel(ncmax))
 allocate(g_bound%ibelong(ncmax))

return
end subroutine allocatebound
!######################################## ALLOCATEPOLY
subroutine allocatepoly(g_poly,lpmax,ncmax,nclose)
implicit none
type(poly_data),intent(out) :: g_poly
integer(4),     intent(in)  :: ncmax, lpmax, nclose

g_poly%lpmax=lpmax
g_poly%ncmax=ncmax
g_poly%lpoly0=0
allocate(g_poly%lpoly(lpmax))
allocate(g_poly%xypoly(2,lpmax,ncmax))
allocate(g_poly%ind2(nclose,2,2))

return
end
!######################################## DEALLOCATEPOLY
subroutine deallocatepoly(g_poly)
implicit none
type(poly_data),intent(inout) :: g_poly

 deallocate(g_poly%lpoly)
 deallocate(g_poly%xypoly)

return
end subroutine deallocatepoly

!######################################## ALLOCATECOAST
subroutine allocatecoast(g_coast,ncmax)
implicit none
type(coast_data),intent(out) :: g_coast
integer(4),      intent(in)  :: ncmax

 g_coast%ncmax = ncmax
 allocate( g_coast%cxy(2,ncmax)  )
 allocate( g_coast%ind(ncmax,2)  )
 allocate( g_coast%ind_r(ncmax,2))
 allocate( g_coast%label(ncmax,2))

return
end subroutine allocatecoast

!######################################## DEALLOCATECOAST
subroutine deallocatecoast(g_coast)
implicit none
type(coast_data),intent(inout) :: g_coast

deallocate( g_coast%cxy)
deallocate( g_coast%ind)
deallocate( g_coast%ind_r)
deallocate( g_coast%label)

return
end subroutine deallocatecoast

!########################################  ALLOCATEGRD
subroutine ALLOCATEGRD(g_data, node)
implicit none
type(grid_data),intent(inout) :: g_data
integer(4),intent(in) :: node

g_data%node = node
allocate ( g_data%lonlatalt(3,node) )
allocate ( g_data%xyz(3,node) )

end subroutine ALLOCATEGRD

!########################################  DEALLOCATEGRD
subroutine DEALLOCATEGRD(g_data)
implicit none
type(grid_data),intent(inout) :: g_data

 deallocate (g_data%lonlatalt)
 deallocate (g_data%xyz)

return
end

end module coastline_data
