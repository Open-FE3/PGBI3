module geofem_util
implicit none
private
integer,parameter,public:: geofem_name_len = 64
integer,public:: kreal
integer,public:: kint
parameter(kreal = selected_real_kind(10))
parameter(kint = 4)

integer, allocatable:: input_data_initiated(:)

type,public:: local_mesh
! node info
integer n_node
real(kind=kreal),pointer:: node(:,:)
! elemental info
integer n_elem
integer,pointer:: elem_type(:)
integer,pointer:: elem(:,:)
! PE info
integer n_neighbor_pe
integer,pointer:: neighbor_pe(:)
integer           n_internal
integer,pointer:: import_index(:)
integer,pointer:: import_node(:)
integer,pointer:: export_index(:)
integer,pointer:: export_node(:)
integer,pointer:: global_node_id(:)
integer,pointer:: global_elem_id(:)
end type local_mesh

! interface structures for input of Analysis
type,public:: node_elem_grp
! enumerate
integer n_enum_grp
character(geofem_name_len),pointer:: enum_grp_name(:)
integer,pointer:: enum_grp_index(:)
integer,pointer:: enum_grp_node(:)

!     ! stepwise
!     integer n_step_grp
!     character(geofem_name_len),pointer:: step_grp_name(:)
!     integer,pointer:: step_grp_index(:)
!     integer,pointer:: step_grp_node(:)
end type node_elem_grp

type,public:: grp_data
! node groups
type(node_elem_grp) node_grp
! elemeng group
type(node_elem_grp) elem_grp

! surface group
integer n_surf_grp
character(geofem_name_len),pointer:: surf_grp_name(:)
integer,pointer:: surf_grp_index(:)
integer,pointer:: surf_grp_node(:,:)
end type grp_data

! interface structures for output of Analysis

! interface structures for output of analysis
type,public:: node_elem_data
integer n
integer n_component ! 成分数
integer,pointer:: n_free(:) ! 成分毎の自由度数
character(geofem_name_len),pointer:: label(:)  ! 成分毎のラベル(任意文字列)
real(kind=kreal),pointer:: data(:,:) ! 節点(また要素) x 自由度
integer,pointer:: id(:)              ! 各節点(要素)のローカルID
! つまり、描画対象になる
! 節点とは全てとは限らず、順も不同
integer,pointer:: mask(:)            ! マスク。詳細は不定
end type node_elem_data

type,public:: node_data_md
integer n
integer n_component ! 成分数
integer,pointer:: n_free(:) ! 成分毎の自由度数
character(geofem_name_len),pointer:: label(:)  ! 成分毎のラベル(任意文字列)
real(kind=kreal),pointer:: data(:,:) ! 節点(また要素) x 自由度
end type node_data_md
!  character(100) geofem_input_file_name
type(local_mesh), pointer:: local_mesh_buf_array(:)
type(grp_data), pointer:: grp_data_buf_array(:)
type(node_elem_data) tmpdat

! functions
integer, external  ::  geofem_strtoi
integer, external  ::  geofem_prepare_file
!####
! subroutines
external geofem_file_init

type tmp_grp_data
integer n
character(geofem_name_len),pointer:: name(:)
integer,pointer:: index(:)
integer,pointer:: item(:)
integer,pointer:: aux(:)
end type tmp_grp_data

integer,public:: geofem_app_rank
integer,public:: geofem_app_comm
integer,public:: geofem_app_size
integer,public:: geofem_is_server

integer geofem_coupler_comm
integer geofem_coupler_pesize
integer geofem_coupler_rank

character(124),public:: geofem_input_file_name
character(124),public:: geofem_output_file_name
character(124),public:: geofem_control_file_name
character(124),public:: geofem_visual_file_name

character(124),allocatable:: geofem_input_files(:)
character(124),allocatable:: geofem_output_files(:)

integer geofem_module_num

type geofem_direct_trans
integer count
real(kind=kreal) time
character(geofem_name_len) label
end type geofem_direct_trans

type(geofem_direct_trans), pointer:: geofem_direct_trans_map(:,:)

type geofem_send_map
integer:: initiated
integer,pointer:: index(:)
integer,pointer:: node(:)
end type geofem_send_map

type(geofem_send_map), pointer:: geofem_send_map_array(:,:)

type geofem_recv_map
integer,pointer:: tnodes(:)     ! local ids
integer,pointer:: selems(:)     ! element ids
integer,pointer:: selemtypes(:) ! element types
integer,pointer:: snodes(:,:)   ! nodes which construct element
! element id x element local node id
real(kind=kreal),pointer:: saxis(:,:)    ! axis node id x (x,y,z)
integer,pointer:: snid(:)      ! serial node ID in the PE
integer,pointer:: seid(:)      ! serial elemnt ID in the PE
end type geofem_recv_map

integer geofem_recv_map_array_init
type(geofem_recv_map), pointer:: geofem_recv_map_array(:,:,:) !my, from, pe

type geofem_trans
integer n
real(kind=kreal),pointer:: data(:,:)
end type geofem_trans

integer:: geofem_send_pending

type geofem_recv_conv
integer,pointer:: pe(:)
integer,pointer:: idx(:)
end type geofem_recv_conv

integer:: saved_visual_data_apid = -1

type(geofem_recv_conv), pointer:: geofem_recv_conv_remote_node(:,:)
type(geofem_recv_conv), pointer:: geofem_recv_conv_remote_elem(:,:)
type(geofem_recv_conv), pointer:: geofem_recv_conv_local(:,:)
!
public geofem_init, geofem_init_md
public geofem_get_input_datum, geofem_get_input_data_md
!public geofem_finalize

include 'mpif.h'

contains
!
!########################################  geofem_init
subroutine geofem_init(is_server, errno)
integer is_server
integer errno
call geofem_init_md(is_server, 1, errno)
end subroutine geofem_init
!########################################  geofem_init_md
subroutine geofem_init_md(is_server, n_module, errno)
integer is_server
integer n_module
integer errno
integer rank
integer i

integer debugging

!call geofem_special_init();
call mpi_init(errno)
call mpi_comm_rank(MPI_COMM_WORLD, rank, errno)

geofem_module_num = n_module

! allocation some regions
allocate(input_data_initiated(n_module), &
local_mesh_buf_array(n_module), &
grp_data_buf_array(n_module))

do i=1,n_module
input_data_initiated(i) = 0
enddo

!write(*,*) 'mpi_comm_rank'
call geofem_file_init(geofem_app_comm, geofem_app_rank, geofem_is_server)
!    call mpi_comm_size(geofem_app_comm, geofem_app_size)
call geofem_arg_init(n_module)

allocate(geofem_input_files(n_module))
allocate(geofem_output_files(n_module))

do i = 1,n_module
call geofem_get_file_name(0, i, len(geofem_input_files(i)), &
geofem_input_files(i))
call geofem_get_file_name(1, i, len(geofem_output_files(i)), &
geofem_output_files(i))
enddo

call geofem_get_file_name(0, 1, len(geofem_input_file_name), &
geofem_input_file_name)
call geofem_get_file_name(1, n_module, len(geofem_output_file_name), &
geofem_output_file_name);
call geofem_get_control_file_name(geofem_control_file_name)

is_server = geofem_is_server

if(n_module .ge. 2) then
call geofem_coupler_init(n_module, geofem_coupler_comm, &
geofem_coupler_pesize, geofem_coupler_rank)
endif

call geofem_debugging(debugging)

if(debugging .ne. 0) then
write(*,*) 'geofem_util: MPI initialization has been succeeded.'
endif

end subroutine geofem_init_md
!
!########################################  geofem_get_input_datum
subroutine geofem_get_input_datum(mesh, grp, errno)
type(local_mesh) mesh
type(grp_data) grp
integer errno

call geofem_get_input_data_md(1, mesh, grp, errno)

end subroutine geofem_get_input_datum
!########################################  geofem_get_input_data_md
subroutine geofem_get_input_data_md(apid, mesh, grp, errno)
integer apid
type(local_mesh) mesh
type(grp_data) grp
integer errno
integer n_module

n_module = geofem_module_num

if(apid .le. 0) then
apid = n_module
endif

if (input_data_initiated(apid) .ne. 0) then
call geofem_set_input_data(mesh, grp, errno)
return
endif

call geofem_get_input_data_inner(apid, errno)
mesh = local_mesh_buf_array(apid)
grp = grp_data_buf_array(apid)
end subroutine geofem_get_input_data_md
!###########################################  geofem_set_input_data
subroutine geofem_set_input_data(mesh, grp, errno)
type(local_mesh) mesh
type(grp_data) grp
integer errno

integer apid

if(saved_visual_data_apid .lt. 0) then
saved_visual_data_apid = 1
endif
apid = saved_visual_data_apid

mesh = local_mesh_buf_array(apid)
grp = grp_data_buf_array(apid)

errno = 0
end subroutine geofem_set_input_data
!########################################### geofem_get _input_data_inner
subroutine geofem_get_input_data_inner(apid, errno)
integer apid
integer errno
integer inpid
character(256) line
integer,pointer:: neighbor_pe(:)
integer,pointer:: global_node_id(:)
integer,pointer:: global_elem_id(:)
real(kind=kreal),pointer:: node(:,:)
integer,pointer:: import_index(:)
integer,pointer:: import_node(:)
integer,pointer:: export_index(:)
integer,pointer:: export_node(:)
integer,allocatable:: tmp_count(:)
integer,allocatable:: tmp_buf(:)
integer,pointer:: elem(:,:)
integer,pointer:: elem_type(:)
integer n
integer i,j
integer n_node
integer n_internal
integer n_neighbor_pe
integer itmp(geofem_name_len)
real(kind=kreal),allocatable:: rtmp(:)
integer n_elem

type(tmp_grp_data):: surf_data

integer,pointer:: surf_grp_node(:,:)

integer debugging
integer check

type(local_mesh),pointer:: local_mesh_buf
type(grp_data),pointer:: grp_data_buf

local_mesh_buf => local_mesh_buf_array(apid)
grp_data_buf => grp_data_buf_array(apid)

call geofem_debugging(debugging)

if(debugging .ne. 0) then
print'(i3, a, a)', geofem_app_rank, &
': opening input file: ', geofem_input_files(apid)
endif
call geofem_open(geofem_input_files(apid), 'r', 0, inpid)

if (inpid .lt. 0) then
call geofem_abort(216, trim(geofem_input_files(apid)) &
//': '//'cant open file')
endif

! dummy
call geofem_getline(inpid, line, 256, errno)
! write(*,*) 'dummy = ', line

call geofem_getint(inpid, 1, n_neighbor_pe, 0, errno)
local_mesh_buf%n_neighbor_pe = n_neighbor_pe

if(debugging .ne. 0) then
write(*,*) 'geofem_util:#', geofem_app_rank, &
'got # of neighbor PE: ', n_neighbor_pe
endif

! read neibe
allocate(neighbor_pe(n_neighbor_pe))
local_mesh_buf%neighbor_pe => neighbor_pe
call geofem_getint(inpid, n_neighbor_pe, neighbor_pe, 0, errno)

if(debugging .ne. 0) then
write(*,*) 'geofem_util:#', geofem_app_rank, &
'succeeded to read neighboring PEs'
write(*,*) (neighbor_pe(i),i=1,n_neighbor_pe)
endif

! mesh data

! node data
call geofem_getint(inpid, 1, n_node, 0, errno)
call geofem_getint(inpid, 1, n_internal, 1, errno)

local_mesh_buf%n_node = n_node
local_mesh_buf%n_internal = n_internal

if(debugging .ne. 0) then
write(*,*) 'geofem_util:#', geofem_app_rank, &
'succeeded to read # of node: ', n_node
write(*,*) 'geofem_util:#', geofem_app_rank, &
'succeeded to read # of internal node', n_internal
endif

! write(*,*) 'node = ', n_node, 'internal = ', n_internal

allocate(global_node_id(n_node), node(n_node,3), rtmp(3))

do i=1,n_node
call geofem_getint(inpid, 1, global_node_id(i), 0, errno)
call geofem_getdouble(inpid, 3, rtmp, 1, errno)
node(i,1:3) = rtmp(1:3)
! write(*,*) 'node = ', node(i, 1), node(i, 2), node(i, 3)
enddo

check = 0

if(debugging .ne. 0) then
write(*,*) 'geofem_util:#', geofem_app_rank, &
'succeeded to read ', n_node, ' node data'
do i=1, n_node
if(global_node_id(i) .le. 0) then
write(*,*) 'geofem_util:#', geofem_app_rank, &
' illegal node id: ', i, ' th: ', global_node_id(i)
check = 1
endif
enddo
if(check .ne. 0) then
write(*,*) 'geofem_util:#', geofem_app_rank, &
' all nodes id seem legal'
endif
endif

deallocate(rtmp)

local_mesh_buf%node => node
local_mesh_buf%global_node_id => global_node_id

! element data
call geofem_getint(inpid, 1, n_elem, 0, errno)
local_mesh_buf%n_elem = n_elem

allocate(global_elem_id(n_elem))
allocate(elem_type(n_elem), tmp_count(n_elem))

call geofem_getint(inpid, n_elem, elem_type, 0, errno)

if(debugging .ne. 0) then
write(*,*) 'geofem_util:#', geofem_app_rank, &
' succeeded to read element types'
endif

n = -1
do i=1,n_elem
call geofem_conv_n_node(elem_type(i), j)
tmp_count(i) = j
if (j .gt. n) then
n = j
endif
enddo

if(debugging .ne. 0) then
write(*,*) 'geofem_util:#', geofem_app_rank, &
' max element size is  ', n
endif

allocate(elem(n_elem, n))
allocate(tmp_buf(n))

! write(*,*) 'num elem = ', n

do i=1,n_elem
call geofem_getint(inpid, 1, itmp, 0, errno)
global_elem_id(i) = itmp(1)
call geofem_getint(inpid, tmp_count(i), tmp_buf, 1, errno)
elem(i, 1:n) = tmp_buf(1:n)
! write(*,*) 'elem = ', elem(i, 1), elem(i, 2), elem(i, 3)
enddo

if(debugging .ne. 0) then
check = 0
!       write(*,*) 'geofem_util:#', geofem_app_rank, &
!            ' succeeded to read element data'
do i=1,n_elem
do j=1,tmp_count(i)
if((elem(i,j) .le. 0) .and. (elem(i,j) .gt. &
size(local_mesh_buf%node))) then
!                write(*,*) 'geofem_util:#', geofem_app_rank, &
!                     ' suspicious element: elem(', i, ',', j, ') = ', &
!                     elem(i,j)
check = 1
endif
enddo
enddo
if(check .eq. 0) then
write(*,*) 'geofem_util:#', geofem_app_rank, &
' all elements seem legal'
endif
endif

deallocate(tmp_count, tmp_buf)

local_mesh_buf%global_elem_id => global_elem_id
local_mesh_buf%n_elem = n_elem
local_mesh_buf%elem_type => elem_type
local_mesh_buf%elem => elem

! read importing data
allocate(import_index(0:n_neighbor_pe))
call geofem_getint(inpid, n_neighbor_pe, import_index(1:), 0, errno)
import_index(0) = 0
n = import_index(n_neighbor_pe)
allocate(import_node(n))
do i = 1, n
call geofem_getint(inpid, 1, j, 0, errno)
import_node(i) = j
enddo

if(debugging .ne. 0) then
write(*,*) 'geofem_util:#', geofem_app_rank, &
' succeeded to read import data'
endif

local_mesh_buf%import_index => import_index
local_mesh_buf%import_node => import_node

! read exporting data
allocate(export_index(0:n_neighbor_pe))
call geofem_getint(inpid, n_neighbor_pe, export_index(1:), 0, errno)
export_index(0) = 0
n = export_index(n_neighbor_pe)
allocate(export_node(n))
do i = 1, n
call geofem_getint(inpid, 1, j, 0, errno)
export_node(i) = j
enddo

if(debugging .ne. 0) then
write(*,*) 'geofem_util:#', geofem_app_rank, &
' succeeded to read export data'
endif

local_mesh_buf%export_index => export_index
local_mesh_buf%export_node => export_node

!
! Group infomations
!

! Node group
call group_init(inpid, grp_data_buf%node_grp, errno)

if(debugging .ne. 0) then
call check_group_data(apid, grp_data_buf%node_grp%n_enum_grp, &
grp_data_buf%node_grp%enum_grp_index, &
grp_data_buf%node_grp%enum_grp_node, 1, 'node')
endif

! Elemntal group
call group_init(inpid, grp_data_buf%elem_grp, errno)

if(debugging .ne. 0) then
call check_group_data(apid, grp_data_buf%elem_grp%n_enum_grp, &
grp_data_buf%elem_grp%enum_grp_index, &
grp_data_buf%elem_grp%enum_grp_node, 2, 'element')
endif

! Surface
call group_init_aux(inpid, 2, surf_data, errno)
grp_data_buf%n_surf_grp = surf_data%n
n = surf_data%n
grp_data_buf%surf_grp_name => surf_data%name
grp_data_buf%surf_grp_index => surf_data%index
n = size(surf_data%item)
allocate(surf_grp_node(2,n))
surf_grp_node(1,:) = surf_data%item(:)
surf_grp_node(2,:) = surf_data%aux(:)
grp_data_buf%surf_grp_node => surf_grp_node

if(debugging .ne. 0) then
call check_group_data(apid, grp_data_buf%n_surf_grp, &
grp_data_buf%surf_grp_index, &
grp_data_buf%surf_grp_node(1,:), 3, 'surface')
endif

if(debugging .ne. 0) then
write(*,*) 'geofem_util:#', geofem_app_rank, &
' succeeded to read ', grp_data_buf%n_surf_grp, &
' element group data'
endif

input_data_initiated(apid) = 1
end subroutine geofem_get_input_data_inner
!########################################## group_init_aux
subroutine group_init_aux(inpid, flag, tmp_data, errno)
integer inpid
integer flag ! if flag .eq. 0 them ENUMERATION else STEPWISE
type(tmp_grp_data) tmp_data
integer errno

integer n_dat
character(geofem_name_len),pointer:: name(:)
integer,pointer:: index(:)
integer,pointer:: item(:)
integer,pointer:: aux(:)

integer i,j,n
integer,allocatable:: tmpi(:)

! enumerate type
call geofem_getint(inpid, 1, n_dat, 0, errno)
allocate(name(n_dat))
allocate(index(0:n_dat))

call geofem_getint(inpid, n_dat, index(1:), 0, errno)
index(0) = 0

n = index(n_dat)

if(flag .eq. 1) then
n = n * 3
endif

allocate(item(n))
! write(*,*) n, ' byte alocated'

if(flag .eq. 2) then
allocate(aux(index(n_dat)))
endif

n = 0
do i=1,n_dat
call geofem_getstr(inpid, name(i), geofem_name_len, 0, errno)
call geofem_fill_blank(name(i), geofem_name_len)
j = index(i) - n
if(flag .eq. 1) then
j = j * 3
endif
allocate(tmpi(j))
call geofem_getint(inpid, j, tmpi, 0, errno)
if(flag .ne. 1) then
item(n + 1:index(i)) = tmpi(1:j)
else
item(n * 3 + 1:index(i) * 3) = tmpi(1:j)
endif
if(flag .eq. 2) then
call geofem_getint(inpid, j, tmpi, 0, errno)
aux(n + 1:index(i)) = tmpi(1:j)
endif
deallocate(tmpi)
n = index(i)
enddo

tmp_data%n = n_dat
tmp_data%name => name
tmp_data%index => index
tmp_data%item => item
tmp_data%aux => aux

end subroutine group_init_aux
!#############################################  check_group_data
subroutine check_group_data(apid, n_data, index, contents, kind, label)
integer apid
integer n_data
integer index(:)
integer contents(:)
integer kind
integer i, j, c
integer maxnum
integer found
character(*) label
integer begin

type(local_mesh),pointer:: local_mesh_buf

local_mesh_buf => local_mesh_buf_array(apid)

if(kind .eq. 1) then
maxnum = local_mesh_buf%n_node
else if(kind .eq. 2) then
maxnum = local_mesh_buf%n_elem
else
maxnum = local_mesh_buf%n_elem
endif

found = 0
begin = 1

do i=1,n_data
c = 0
do j = begin, index(i)
if((contents(j) .le. 0) .or. (contents(j) .gt. maxnum)) then
!             write(*,*) 'geofem_util:#', geofem_app_rank, &
!                  'illegal data in ', label, c,  &
!                  'th data of group data(', i, ')  is illegal'
found = 1
endif
c = c + 1
enddo
begin = index(i)
enddo
end subroutine check_group_data


end module geofem_util

