!-----------------------------------------------------------------------
!   Copyright 1999, 2000, 2001, 2002, 2003 by the Research Organization
!   for Information Science & Technology (RIST)
!-----------------------------------------------------------------------

! module geofem_utilg

! $Id: geofem_util.f90,v 1.10 2001/01/26 06:56:41 sekita Exp $

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

  ! subroutines
  external geofem_file_init
  external geofem_open
  external geofem_close
  external geofem_printf
  external geofem_puts
  external geofem_write
  external geofem_gets
  external geofem_read
  external geofem_file_terminate
  external geofem_getenv
  external geofem_strtoni

  external geofem_getint
  external geofem_getdouble
  external geofem_getstr
  external geofem_getline

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

  public geofem_init, geofem_init_md
  public geofem_get_input_datum, geofem_get_input_data_md
  public geofem_finalize
  public geofem_put_data_md
!  public geofem_put_to_visual, geofem_put_to_visual_md
  public geofem_get_data_md

  public geofem_get_target_axis
  public geofem_get_source_node
  public geofem_get_source_axis
  public geofem_get_all_target_num
  public geofem_get_all_target_axis
  public geofem_get_all_source_elem_num
  public geofem_get_all_source_elem
  public geofem_get_all_source_axis_num
  public geofem_get_all_source_axis
  public geofem_get_all_source_elems
  public geofem_get_all_source_elemtypes
  public geofem_get_nmodule
  public geofem_set_input_data_c
  public geofem_set_input_datum_c

  include 'mpif.h'

contains

  subroutine geofem_init(is_server, errno)
    integer is_server
    integer errno
    call geofem_init_md(is_server, 1, errno)
  end subroutine geofem_init

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

  subroutine geofem_coupler_init(n_module, comm, pesize, rank)
    integer n_module
    integer comm, pesize, rank
    
    integer i, j, k
    integer errno

    call geofem_coupler_init_c(n_module, comm, pesize, rank, errno)

    allocate(geofem_direct_trans_map(n_module, n_module))

    allocate(geofem_send_map_array(n_module, n_module))
    do i=1,n_module
       do j=1,n_module
          geofem_send_map_array(i,j)%initiated = 0
          nullify(geofem_send_map_array(i, j)%index)
       enddo
    enddo

    allocate(geofem_recv_map_array(n_module, n_module, geofem_coupler_pesize))

    do i=1, n_module
       do j=1, n_module
          do k =1, geofem_coupler_pesize
             nullify(geofem_recv_map_array(i, j, k)%tnodes)
          enddo
       enddo
    enddo

    allocate(geofem_recv_conv_remote_node(n_module, n_module))
    allocate(geofem_recv_conv_remote_elem(n_module, n_module))
    allocate(geofem_recv_conv_local(n_module, n_module))
  end subroutine geofem_coupler_init


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

  subroutine geofem_set_grp_name(flag, enum_name)
    integer flag
    character(geofem_name_len), dimension(:):: enum_name
    external geofem_set_grp_name_aux
    external geofem_set_grp_name_init

    integer n
    integer i

    n = size(enum_name)

    call geofem_set_grp_name_init(flag, 0, n)

    do i=1,n
       call geofem_set_grp_name_aux(flag, 0, i, enum_name(i))
    enddo

  end subroutine geofem_set_grp_name

  subroutine geofem_set_surf_name(name)
    character(geofem_name_len),dimension(:):: name
    
    external geofem_set_grp_name_aux
    external geofem_set_grp_name_init

    integer n
    integer i

    n = size(name)

    call geofem_set_grp_name_init(3, 0, n)

    !write(*,*) 'surf_name'

    do i=1,n
       call geofem_set_grp_name_aux(3, 0, i, name(i))
    enddo
  end subroutine geofem_set_surf_name

  ! backward compatibility
  subroutine geofem_set_input_datum_c(errno)
    integer errno
    call geofem_set_input_data_c(errno)
  end subroutine geofem_set_input_datum_c

  subroutine geofem_set_input_data_c(errno)
    integer apid
    integer errno
    type(local_mesh), pointer:: local_mesh_buf
    type(grp_data), pointer:: grp_data_buf

    external geofem_set_mesh_ints
    external geofem_set_mesh_arrays
    external geofem_set_grp_info1
    external geofem_set_grp_info2

!    write(*,*) 'saved_visual_data_apid = ', saved_visual_data_apid
    if(saved_visual_data_apid .lt. 0) then
       saved_visual_data_apid = 1
    endif
 
    apid = saved_visual_data_apid
!    write(*,*) 'apid = ', apid

    local_mesh_buf => local_mesh_buf_array(apid)
    grp_data_buf => grp_data_buf_array(apid)

    ! write(*,*) local_mesh_buf%n_node, local_mesh_buf%n_elem
    if(input_data_initiated(apid) .eq. 0) then
       call geofem_get_input_data_inner(apid, errno)
    endif
    ! write(*,*) 'input_datum_inner'
    ! transfer mesh infos

    call geofem_set_mesh_ints(local_mesh_buf%n_node, &
         local_mesh_buf%n_elem, &
         local_mesh_buf%n_neighbor_pe, &
         local_mesh_buf%n_internal)
    !write(*,*) 'ints'

    call geofem_set_mesh_arrays( &
         size(local_mesh_buf%node),local_mesh_buf%node, &
         size(local_mesh_buf%elem_type), local_mesh_buf%elem_type, &
         size(local_mesh_buf%elem, 2), &
         size(local_mesh_buf%elem), local_mesh_buf%elem, &
         size(local_mesh_buf%neighbor_pe), local_mesh_buf%neighbor_pe, &
         size(local_mesh_buf%import_index), local_mesh_buf%import_index, &
         size(local_mesh_buf%import_node), local_mesh_buf%import_node, &
         size(local_mesh_buf%export_index), local_mesh_buf%export_index, &
         size(local_mesh_buf%export_node), local_mesh_buf%export_node, &
         size(local_mesh_buf%global_node_id), local_mesh_buf%global_node_id, &
         size(local_mesh_buf%global_elem_id), local_mesh_buf%global_elem_id)

    !write(*,*) 'arrays'

    ! transfer group infos
    call geofem_set_grp_info1(0, &
         grp_data_buf%node_grp%n_enum_grp, &
         grp_data_buf%node_grp%enum_grp_index, &
         grp_data_buf%node_grp%enum_grp_node)

    !write(*,*) 'group 1'

    call geofem_set_grp_info1(1, &
         grp_data_buf%elem_grp%n_enum_grp, &
         grp_data_buf%elem_grp%enum_grp_index, &
         grp_data_buf%elem_grp%enum_grp_node)

    !write(*,*) 'group 2'

    call geofem_set_grp_info2(grp_data_buf%n_surf_grp, &
         grp_data_buf%surf_grp_index, &
         grp_data_buf%surf_grp_node)

    !write(*,*) 'group 3'

    call geofem_set_grp_name(0, grp_data_buf%node_grp%enum_grp_name)

    !write(*,*) 'group 4'

    call geofem_set_grp_name(1, grp_data_buf%elem_grp%enum_grp_name)

    !write(*,*) 'group 5'

    call geofem_set_surf_name(grp_data_buf%surf_grp_name)

    !write(*,*) 'group 6'

  end subroutine geofem_set_input_data_c

  subroutine geofem_get_input_datum(mesh, grp, errno)
    type(local_mesh) mesh
    type(grp_data) grp
    integer errno

    call geofem_get_input_data_md(1, mesh, grp, errno)

  end subroutine geofem_get_input_datum

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

  subroutine group_init(inpid, g_buf, errno)
    type(node_elem_grp) g_buf
    integer inpid

    type(tmp_grp_data) tmp_data
    integer errno

    call group_init_aux(inpid, 0, tmp_data, errno)
    g_buf%n_enum_grp = tmp_data%n
    g_buf%enum_grp_name => tmp_data%name
    g_buf%enum_grp_index => tmp_data%index
    g_buf%enum_grp_node => tmp_data%item

  end subroutine group_init

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

  subroutine geofem_finalize
    integer ierr
    call geofem_file_terminate
    call MPI_Finalize(ierr)
  end subroutine geofem_finalize

!-----------
! put data
!-----------

  subroutine geofem_send_data_init(myid, toid, count, time, label)
    integer myid
    integer toid
    integer count
    real(kind=kreal) time
    character(*) label

    integer, pointer:: index(:)
    integer, pointer:: node(:)
    integer err
    type(geofem_send_map),pointer:: send_map

    call geofem_get_sender_data_init(myid, toid, err) ! init C region
    if(geofem_send_map_array(myid, toid)%initiated .eq. 0) then
       if(err .gt. 0) then
          allocate(index(0:geofem_coupler_pesize))
          send_map => geofem_send_map_array(myid, toid)
          send_map%index => index

          call geofem_get_sender_data_n(index)

          allocate(node(1:index(geofem_coupler_pesize)))
          send_map%node => node
          call geofem_get_sender_data_node(node)
       endif
       geofem_send_map_array(myid, toid)%initiated = 1
    endif

    call geofem_send_data_init_c(myid, toid, count, time, len(label), &
         label)

  end subroutine geofem_send_data_init

  subroutine geofem_send_data_label(lsize, labels)
    integer lsize
    character(*) labels(:)
    integer i
    call geofem_send_label_size(lsize)
    do i=1, lsize
       call geofem_send_label(i, len(labels(i)), labels(i))
    enddo
  end subroutine geofem_send_data_label

  subroutine geofem_put_data_md(myid, toid, count, time, label, nndat, errno)
    integer myid, toid, count
    real(kind=kreal) time
    character(*) label
    type(node_data_md) nndat
    integer errno
    integer pe
    integer isize

    call geofem_send_data_init(myid, toid, count, time, label)

    isize = size(nndat%label)
    call geofem_send_data_label(isize, nndat%label)
    call geofem_send_data_head(nndat%n, nndat%n_component, nndat%n_free)
    do pe = 1, geofem_coupler_pesize
       call geofem_send_data(pe - 1, &
            size(nndat%data, 1), size(nndat%data, 2), nndat%data)
       call geofem_really_send_data
    enddo
  end subroutine geofem_put_data_md

!----------
! get data
!----------

subroutine geofem_reader_data_init(myid, fromid)
  integer myid
  integer fromid

  integer n_tnodes, n_selems, n_maxnode, n_axis

  integer i, j, k

  integer n_node

  integer datasize
  type(geofem_recv_conv), pointer:: recv_conv
  type(geofem_recv_map),pointer:: recv_map

  integer size;
  
  do i = 1, geofem_coupler_pesize
     call geofem_reader_data_one_pe(myid, fromid, &
          i, n_tnodes, n_selems, n_maxnode, &
          n_axis)
     recv_map => geofem_recv_map_array(myid, fromid, i)
     if(associated(recv_map%tnodes)) then
        deallocate(recv_map%tnodes, recv_map%selems, &
             recv_map%selemtypes, recv_map%snodes, recv_map%saxis, &
             recv_map%snid, recv_map%seid)
     endif
     allocate(recv_map%tnodes(n_tnodes), recv_map%selems(n_tnodes), &
          recv_map%selemtypes(n_selems), &
          recv_map%snodes(n_selems, n_maxnode), &
          recv_map%saxis(n_axis, 3), &
          recv_map%snid(n_axis), recv_map%seid(n_selems))
     call geofem_reader_copy_data(recv_map%tnodes, &
          recv_map%selems, &
          recv_map%selemtypes, &
          recv_map%snodes, &
          recv_map%saxis)
!     do j=1,size(recv_map%selemtypes)
!        write(*,*) 'PE#', geofem_app_rank, i, ',', j, ': selemtypes = ', recv_map%selemtypes(j)
!     enddo
  enddo

!  write(*,*) 'geofem_reader_copy_data'

  ! serialisze node
  datasize = 0
  do i=1, geofem_coupler_pesize
     recv_map => geofem_recv_map_array(myid, fromid, i)
     datasize = datasize + size(recv_map%snid)
!     write(*,*) 'datasize = ', datasize
  enddo

  recv_conv => geofem_recv_conv_remote_node(myid, fromid)
  allocate(recv_conv%pe(datasize), recv_conv%idx(datasize))
!  write(*,*) 'recv_conv allocated'

  k = 1
  do i=1, geofem_coupler_pesize
     recv_map => geofem_recv_map_array(myid, fromid, i)
     datasize = size(recv_map%snid)
     do j = 1, datasize
!        write(*,*) i, ':', myid, ':', fromid, 'snid(', j, ') = ', k
        recv_map%snid(j) = k
        recv_conv%pe(k) = i
        recv_conv%idx(k) = j
        k = k + 1
     enddo
  enddo

!  write(*,*) 'serialize node'

  ! serialize element
  datasize = 0
  do i=1, geofem_coupler_pesize
     recv_map => geofem_recv_map_array(myid, fromid, i)
     datasize = datasize + size(recv_map%seid)
  enddo

  ! remote node transfer table
  recv_conv => geofem_recv_conv_remote_elem(myid, fromid)
  allocate(recv_conv%pe(datasize), recv_conv%idx(datasize))

  k = 1 
  do i=1, geofem_coupler_pesize ! PE loop
     recv_map => geofem_recv_map_array(myid, fromid, i)
     datasize = size(recv_map%seid)
     do j=1, datasize
        recv_map%seid(j) = k
        recv_conv%pe(k) = i
        recv_conv%idx(k) = j
        k = k + 1
     enddo
  enddo

  !write(*,*) 'serialize elem'

  ! local node transfer table
  recv_conv => geofem_recv_conv_local(myid, fromid)
  n_node = local_mesh_buf_array(myid)%n_node
  allocate(recv_conv%pe(n_node), recv_conv%idx(n_node))

!  write(*,*) 'n_node = ', n_node
  do i=1, n_node
     recv_conv%pe(i) = -1
  enddo

  ! 'tnodes' contains a target node ID in the local PE.
  do i=1, geofem_coupler_pesize
     recv_map => geofem_recv_map_array(myid, fromid, i)
     datasize = size(recv_map%tnodes)
!     write(*,*) 'PE#', i, ': tnode size = ', datasize
     do j=1, datasize
        k = recv_map%tnodes(j)
!        if(recv_conv%pe(k) .ne. -1) then
!           call GEOFEM_abort(216, 'Recieve map data is redundant')
!        endif
!        write(*,*) 'tnodes(', j, ') = ', k
        recv_conv%pe(k) = i
        recv_conv%idx(k) = j
     enddo
  enddo

  !write(*,*) 'local node transder'

end subroutine geofem_reader_data_init

subroutine geofem_get_data_md(myid, fromid, count, time, label, nndat, &
     interpf, errno)
  integer myid, fromid, count
  real(kind=kreal) time
  character(*)label
  type(node_data_md) nndat
  interface ! 以下は入力
     subroutine interpf(myid, fromid, stat, dim, tid, elemtype, elemid, &
          outstat, errno)
       parameter(kreal = selected_real_kind(10))
       integer myid ! 自己APID
       integer fromid ! 転送元APID
       integer tid  ! 内挿対象節点のID
       integer dim  ! 内挿対象点の次元(1/2 = 一次/二次)
       real(kind=kreal) stat(:,:) ! 状態量(全自由度分一度に与えられる)
       ! (node,自由度)
       integer elemtype  ! 内挿元要素の型
       integer elemid    ! 内挿元要素のID
       real(kind=kreal) outstat(:) ! 内挿された状態量。これだけが出力
       integer errno
     end subroutine interpf
  end interface
  integer errno

  type(geofem_trans),pointer:: trans(:)
  type(geofem_trans),pointer:: tmp_trans
  type(geofem_recv_conv),pointer:: local_conv
  type(geofem_recv_map),pointer:: recv_map
  real(kind=kreal),pointer:: stat(:,:)
  real(kind=kreal),pointer:: outstat(:)
  integer i, j, k, statsize, n_component
  integer pe, idx
  integer elem
  integer n_node, nodesize
  integer n_internal
  integer nlabel

  allocate(trans(geofem_coupler_pesize))

  call geofem_reader_data_init(myid, fromid)

  !write(*,*) 'geofem_reader_data_init'

  call geofem_recieve_data_init(myid, fromid, trans)
  !write(*,*) 'geofem_recieve_data_init'
  call geofem_recieve_data_init_c(myid, fromid, count, time, label)

  call geofem_recieve_data_head(nndat%n_component, nlabel)
  allocate(nndat%n_free(nndat%n_component))
  allocate(nndat%label(nlabel))
  do i = 1, nlabel
     call geofem_recieve_label(i, nndat%label(i))
!     print '(i3, a, a)', i, ': label = ', nndat%label(i)
  enddo
  call geofem_recieve_data_head_n_free(nndat%n_free)

  !call MPI_barrier(geofem_coupler_comm, i)

  statsize = 0
  n_component = nndat%n_component
  do i = 1, n_component
     statsize = statsize + nndat%n_free(i)
  enddo

  ! interpolation
  local_conv => geofem_recv_conv_local(myid, fromid)

  nodesize = local_mesh_buf_array(myid)%n_node
  nndat%n = nodesize
  n_internal = local_mesh_buf_array(myid)%n_internal

  allocate(nndat%data(nodesize, statsize))

  do i=1, n_internal
     pe = local_conv%pe(i)
     idx = local_conv%idx(i)
!     write(*,*) 'pe = ', pe, 'idx = ', idx
     if(pe .eq. -1) then
        call GEOFEM_abort(216, 'Recieve map lacks data')
     endif
     recv_map => geofem_recv_map_array(myid, fromid, pe)
     tmp_trans => trans(pe)

!     write(*,*) 'size of selems = ', size(recv_map%selems)

     elem = recv_map%selems(idx)

!     write(*,*) 'PE#', geofem_app_rank, ', selms(', idx, ') = ', elem
!     write(*,*) 'selms = ', recv_map%selems
!     write(*,*) 'selemtypes(', elem, ') = ', recv_map%selemtypes(elem)
!     write(*,*) 'selemtypes = ', recv_map%selemtypes

     call GEOFEM_conv_n_node(recv_map%selemtypes(elem), n_node)
!     write(*,*) 'size: tnodes = ', size(recv_map%tnodes), ' selemtypes = ',  &
!     size(recv_map%selemtypes)

!          do j=1,size(recv_map%selemtypes)
!             write(*,*) '#', j, recv_map%selemtypes(j)
!          enddo

          !write(*,*) '#', i, ' node: ', n_node, 'x', statsize
     allocate(stat(n_node, statsize), outstat(statsize))

     !call MPI_Barrier(geofem_coupler_comm, j)

!     write(*,*) 'PE#', geofem_app_rank
     do j = 1, n_node
        do k = 1, statsize
           stat(j, k) = tmp_trans%data(recv_map%snodes(elem, j), k)
!           write(*,*) 'stat(', j, ',', k, ') = ', stat(j, k)
        enddo
     enddo

!     write(*,*) '#', elem, ' elem: type = ', recv_map%selemtypes(elem)
     call interpf(myid, fromid, stat, 1, i, recv_map%selemtypes(elem), &
          recv_map%seid(elem), outstat, errno)

     nndat%data(i,:) = outstat(:)
     deallocate(stat, outstat)
  enddo

!  write(*,*) 'exchange_interp_data'

  call exchange_interp_data(myid, nndat)
!  write(*,*) 'PE#', geofem_coupler_rank, ': end of geofem_get_data_md'
  !call MPI_Barrier(geofem_coupler_comm, i)

!  write(*,*) 'PE#', geofem_coupler_rank, ': barriered'
  
end subroutine geofem_get_data_md

subroutine exchange_interp_data(myid, nndat)
  integer myid
  type(node_data_md) nndat

  type(local_mesh),pointer:: mesh
  real(kind=kreal),allocatable:: buf(:,:)
  integer, allocatable:: reqsend1(:), reqsend2(:)
  integer npe, dsize, pe, begin, end, i, j, num, idx, found
  integer stat(MPI_STATUS_SIZE)
  integer, allocatable:: stat1(:,:), stat2(:,:)
  integer ierr 

  mesh => local_mesh_buf_array(myid)

  npe = mesh%n_neighbor_pe
  dsize = size(nndat%data, 2)
  allocate(reqsend1(npe), reqsend2(npe))
  allocate(stat1(MPI_STATUS_SIZE, npe), stat2(MPI_STATUS_SIZE, npe))

  ! Export
  begin = mesh%export_index(0)
  do i=1, npe
     pe = mesh%neighbor_pe(i)
!     write(*,*) 'sending to PE#', pe, 'from', geofem_coupler_rank
     end = mesh%export_index(i)
     num = end - begin
     begin = begin + 1
     allocate(buf(num, dsize))
     do j=begin, end
        idx = mesh%export_node(j)
!        write(*,*) '#', geofem_coupler_rank, ': sending #', idx, ' node'
        buf(j - begin + 1, :) = nndat%data(idx, :)
     enddo
     !write(*,*) 'geofem_coupler_comm = ', geofem_coupler_comm
     call MPI_Isend(num, 1, MPI_INTEGER, pe, 130, geofem_coupler_comm, &
          reqsend1(i), ierr)
     call MPI_Isend(buf, num * dsize, MPI_DOUBLE_PRECISION, &
          pe, 130, geofem_coupler_comm, reqsend2(i), ierr)
!     write(*,*) '#', geofem_coupler_rank, ': sent ', num, 'data'
     begin = end
     deallocate(buf)
  enddo

  !call MPI_Barrier(geofem_coupler_comm, ierr)
!  write(*,*) 'PE#', geofem_coupler_rank, 'ALL SENT'

  ! Import
  begin = mesh%import_index(0)
  do i = 1, npe
     call MPI_Recv(num, 1, MPI_INTEGER, MPI_ANY_SOURCE, 130, &
          geofem_coupler_comm, stat, ierr)
     pe = stat(MPI_SOURCE)
!     write(*,*) 'PE#', geofem_coupler_rank, ': recieving from PE#', pe
     found = -1
     do j = 1, npe
        if(mesh%neighbor_pe(j) .eq. pe) then
           found = j
        endif
     enddo
     if(found .lt. 0) then
        call GEOFEM_abort(216, 'exchange_interp_data: message from invalid PE')
     endif
     begin = mesh%import_index(found-1)
     end = mesh%import_index(found)
!     write(*,*) 'PE#', geofem_coupler_rank, ': ', found, 'th PE'
!     write(*,*) 'PE#', geofem_coupler_rank, ': begin = ', begin, 'end = ',end
     if((end - begin) .ne. num) then
!        write(*,*) 'expecting: ', end - begin, ' recieved: ', num
        call GEOFEM_abort(216, 'exchange_interp_data: # of data inconsistent')
     endif
     allocate(buf(num, dsize))
     call MPI_Recv(buf, num * dsize, MPI_DOUBLE_PRECISION, &
          pe, 130, geofem_coupler_comm, stat, ierr)
!     write(*,*) 'PE#', geofem_coupler_rank, ': ', num, 'data recieved from pe#', pe
     begin = begin + 1
     do j = begin, end
        idx = mesh%import_node(j)
!        write(*,*) 'PE#', geofem_coupler_rank, ': idx = ', idx
        nndat%data(idx, :) = buf(j - begin + 1, :)
     enddo
     deallocate(buf)
  enddo


  !call MPI_Barrier(geofem_coupler_comm, ierr)
  !write(*,*) 'PE#', geofem_coupler_rank, ': ALL RECIEVED'

  call MPI_Waitall(npe, reqsend1, stat1, ierr)
  call MPI_Waitall(npe, reqsend2, stat2, ierr)
  !write(*,*) 'PE#', geofem_coupler_rank, ': ALL WAITED'
end subroutine exchange_interp_data

subroutine geofem_recieve_data_init(myid, fromid, trans_buf)
  integer myid, fromid
  type(geofem_trans):: trans_buf(:)

!  type(geofem_send_map), pointer:: send_map
!  type(geofem_recv_map), pointer:: recv_map
  integer i
  integer,allocatable:: n_num(:)
  integer,allocatable:: n_stat(:)

  call geofem_really_recieve_data(myid, fromid)

  allocate(n_num(geofem_coupler_pesize), n_stat(geofem_coupler_pesize))

  call geofem_recieve_data_num(myid, fromid, n_num, n_stat)

  do i=1, geofem_coupler_pesize
     trans_buf(i)%n = n_num(i)
!     write(*,*) 'n_num = ', n_num(i)
!     if(n_stat(i) .ne. statnum) then
!        call GEOFEM_abort(216, '# of status not agree with sender ')
!     endif
     allocate(trans_buf(i)%data(n_num(i), n_stat(i)))
     call geofem_recieve_data(myid, fromid, i - 1, trans_buf(i)%data)
  enddo
end subroutine geofem_recieve_data_init

!subroutine geofem_put_to_visual(count, time, label, nndat, nedat, errno)
!  integer count
!  real(kind=kreal) time
!  character(*) label
!  type(node_elem_data) nndat
!  type(node_elem_data) nedat
!  integer errno
!
!  call geofem_put_to_visual_md(geofem_module_num, count, time, label, &
!       nndat, nedat, errno)
!
!end subroutine geofem_put_to_visual

!!!!!!!!!!!!!!!!!!!!! Utilities !!!!!!!!!!!!!!!!!!!!!

subroutine geofem_get_nmodule(n)
  integer(kind=kint) n
  n = geofem_module_num
end subroutine geofem_get_nmodule

subroutine geofem_get_target_axis(myid, tid, taxis)
  integer myid
  integer tid
  real(kind=kreal) taxis(3)
  
  type(local_mesh),pointer:: mesh
  
  mesh => local_mesh_buf_array(myid)
  
  taxis(:) = mesh%node(tid, :)
end subroutine geofem_get_target_axis

  subroutine geofem_get_all_source_elems(myid, fromid, eid)
    integer myid
    integer fromid
    integer eid(:)
    integer n_internal
    type(geofem_recv_conv),pointer:: local_conv
    type(geofem_recv_map), pointer:: recv_map
    
    integer i, pe, idx

    n_internal = local_mesh_buf_array(myid)%n_internal
    local_conv => geofem_recv_conv_local(myid, fromid)

    do i=1, n_internal
       pe = local_conv%pe(i)
       idx = local_conv%idx(i)
       if(pe .eq. -1) then
          call GEOFEM_abort(216, 'Recieve map lacks data')
       endif
       recv_map => geofem_recv_map_array(myid, fromid, pe)
       eid(i) = recv_map%seid(recv_map%selems(idx))
    enddo
  end subroutine geofem_get_all_source_elems

  subroutine geofem_get_all_source_elemtypes(myid, fromid, etypes)
    integer myid
    integer fromid
    integer etypes(:)
    integer n_internal
    type(geofem_recv_conv),pointer:: local_conv
    type(geofem_recv_map),pointer:: recv_map
    integer i, pe, idx

    n_internal = local_mesh_buf_array(myid)%n_internal
    local_conv => geofem_recv_conv_local(myid, fromid)

    do i=1, n_internal
       pe = local_conv%pe(i)
       idx = local_conv%idx(i)
       if(pe .eq. -1) then
          call GEOFEM_abort(216, 'Recieve map lacks data')
       endif
       recv_map => geofem_recv_map_array(myid, fromid, pe)
       etypes(i) = recv_map%selemtypes(recv_map%selems(idx))
    enddo
  end subroutine geofem_get_all_source_elemtypes

  subroutine geofem_get_source_node(myid, fromid, eid, snode)
    integer myid, fromid, eid
    integer snode(:)
    integer pe, idx

    type(geofem_recv_conv),pointer:: conv
    type(geofem_recv_map),pointer:: recv
    integer,pointer:: nodes(:,:)
    integer i
    integer etype, nsize

    conv => geofem_recv_conv_remote_elem(myid, fromid)

    pe = conv%pe(eid)
    idx = conv%idx(eid)

    recv => geofem_recv_map_array(myid, fromid, pe)
    nodes => recv%snodes
    
    if(size(snode) .lt. size(nodes, 2)) then
       call geofem_abort(216, 'geofem_get_source_code: sizeof node array wrong')
    endif

    etype = recv%selemtypes(idx)
    call geofem_conv_n_node(etype, nsize)

    do i=1,nsize
!       write(*,*) 'recv%snid(', nodes(idx, i), ') = ', recv%snid(nodes(idx, i))
       snode(i) = recv%snid(nodes(idx, i))
    enddo
  end subroutine geofem_get_source_node

  subroutine geofem_get_source_axis(myid, fromid, snode, saxis)
    integer myid, fromid
    integer snode(:)
    real(kind=kreal) saxis(:,:)
    integer node_size
    integer pe, idx, i, n
    integer sn
    type(geofem_recv_conv),pointer:: conv
    type(geofem_recv_map),pointer:: recv_map

    conv => geofem_recv_conv_remote_node(myid, fromid)

    n = size(snode)
    node_size = size(conv%pe)
    
    do i=1, n
       sn = snode(i)
       if((sn .le. 0) .or. (sn .gt. node_size)) then 
          cycle
       endif
       pe = conv%pe(sn)
       idx = conv%idx(sn)
       recv_map => geofem_recv_map_array(myid, fromid, pe)

!       write(*,*) 'dest: ', size(saxis, 1), 'x', size(saxis, 2)  
!       if(associated(recv_map%saxis)) then
!          write(*,*) 'src:  ', size(recv_map%saxis, 1), 'x', size(recv_map%saxis, 2)
!       else
!          write(*,*) 'recv_map%saxis is not allocated'
!       endif

       saxis(i,:) = recv_map%saxis(idx, :)
    enddo
  end subroutine geofem_get_source_axis

  subroutine geofem_get_all_target_num(myid, num)
    integer myid
    integer num ! これは内部点と同じ数になるはず。

    type(local_mesh),pointer:: lm
    lm => local_mesh_buf_array(myid)
    num = lm%n_internal
  end subroutine geofem_get_all_target_num

  subroutine geofem_get_all_target_axis(myid, axis)
    integer myid
    real(kind=kreal) axis(:,:)
    integer i, n
    type(local_mesh),pointer:: lm

    lm => local_mesh_buf_array(myid)
    n = lm%n_internal
    
    do i=1,n
       axis(i,:) = lm%node(i,:)
    enddo
  end subroutine geofem_get_all_target_axis

  subroutine geofem_get_all_source_elem_num(myid, fromid, num)
    integer myid, fromid, num
    type(geofem_recv_conv),pointer:: recv
    recv => geofem_recv_conv_remote_elem(myid, fromid)
    num = size(recv%pe)
  end subroutine geofem_get_all_source_elem_num

  subroutine geofem_get_all_source_elem(myid, fromid, elemtype, node)
    integer myid, fromid
    integer elemtype(:)
    integer node(:,:)

    integer i, n, pe, idx

    type(geofem_recv_conv),pointer:: recv
    type(geofem_recv_map),pointer:: recv_map

    recv => geofem_recv_conv_remote_elem(myid, fromid)

    n = size(recv%pe)

    do i=1, n
       pe = recv%pe(i)
       idx = recv%idx(i)
       recv_map => geofem_recv_map_array(myid, fromid, pe)
       elemtype(i) = recv_map%selemtypes(idx)
       node(i,:) = recv_map%snodes(idx,:)
    enddo
  end subroutine geofem_get_all_source_elem

  subroutine geofem_get_all_source_axis_num(myid, fromid, num)
    integer myid, fromid, num
    integer i, n
    type(geofem_recv_map),pointer:: recv_map

    n = geofem_coupler_pesize
    num = 0

    do i=1, n
       recv_map => geofem_recv_map_array(myid, fromid, i)
       num = num + size(recv_map%saxis, 1)
    enddo
  end subroutine geofem_get_all_source_axis_num

  subroutine geofem_get_all_source_axis(myid, fromid, axis)
    integer myid, fromid
    real(kind=kreal) axis(:,:)
    integer i, j, n
    type(geofem_recv_conv),pointer:: recv_conv
    type(geofem_recv_map),pointer:: recv_map
    integer pe, idx

    recv_conv => geofem_recv_conv_remote_elem(myid, fromid)

    n = size(recv_conv%pe)

    do i=1, n
       pe = recv_conv%pe(i)
       idx = recv_conv%idx(i)
       recv_map => geofem_recv_map_array(myid, fromid, pe)
       axis(i, :) = recv_map%saxis(idx, :)
    enddo
  end subroutine geofem_get_all_source_axis

!include 'geofem_put_to_visual.inc'

!    subroutine dump_grp(grp)
!      type(node_elem_grp)grp

!      integer i
!      do i=1,grp%n_enum_grp
!         write(*,*) 'name: ', grp%enum_grp_name(i)
!         write(*,*) 'index: ', grp%enum_grp_index(i)
!      enddo

!      do i=1,grp%enum_grp_index(grp%n_enum_grp)
!         write(*,*) 'data: ', grp%enum_grp_node(i)
!      enddo
!    end subroutine dump_grp
end module geofem_util
