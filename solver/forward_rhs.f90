!#################################################################### GENMAT8
! Coded on May 16, 2016
! to calculate right hand side vector
!subroutine forward_rhs(A,em_mesh,nodes,l_line,b_vec,dmu,fsn1,Fxyz,vxyz,igroup,dt,it)
subroutine forward_rhs(A,em_mesh,h_ocean,l_line,b_vec,fsn0,fsn1,igroup,dt,it,g_param)
use  outerinnerproduct
use  iccg_var_takuto ! b_vec is not included, see m_iccg_var_takuto.f90
use  mesh_type       ! see m_mesh_type.f90
use  line_type       ! see m_line_type.f90
use  fem_util        ! for volume, intv, (see m_fem_utiil.f90 )
use  fem_edge_util   ! see fem_edge_util.f90
use  constants       ! for dmu (see m_constants.f90)
use  param           ! see m_param.f90
use  oceanFvxyz      ! for Fxyz, vxyz (see m_oceanFxyz.f90)
use spherical        ! see m_spherical.f90, added on 2016.11.19
implicit none
type(ocean_data),   intent(in) :: h_ocean
type(param_forward),intent(in) :: g_param
type(global_matrix),intent(in) :: A ! only for Avalue_bc, line_bc, bcorrect
type(mesh),intent(in)          :: em_mesh
integer(4),intent(in)          :: igroup(em_mesh%ntet),it
type(line_info),intent(in)     :: l_line
real(8),   intent(in)          :: dt
!real(8),   intent(in) :: vxyz(3,nodes) ! [mm/s]
!real(8),   intent(in) :: Fxyz(3,nodes) ! [nT]
real(8),   intent(in) :: fsn0(l_line%nline) ! [nT*km] added on July 19, 2016
real(8),   intent(in) :: fsn1(l_line%nline) ! [nT*km]
real(8),   intent(out):: b_vec(l_line%nline)
real(8)    :: elm_xyz(3,4),xx(3,6), gn(3,4), vF(3,4)
real(8),allocatable,dimension(:,:)  :: vxyz,Fxyz
real(8)    :: w(6,3),  v, sigma
real(8)    :: elm_k(6,6), S(6,6), S1(6,6),S2(6),F(3)
integer(4) :: n(4),node
integer(4) :: iele, i, j, k, l, ii, jj, idirection(6),id_group
!---------------  scales ----------------------------------------------
real(8), parameter  :: L0=1.d+3  ! [m]  scale length

!#[0]## set input
node = em_mesh%node
allocate(vxyz(3,node),Fxyz(3,node))
Fxyz = h_ocean%Fxyz ! [nT]
vxyz = h_ocean%vxyz ! [mm/sec]
b_vec(:)=0.d0

! tetrahedral mesh in sea_mesh is the same at the beginning of em_mesh
do iele=1, em_mesh%ntet
  sigma = em_mesh%cmodel(iele) ! 2018.08.31

!#
  !# [1] ## ! check the direction of edge, compared to the defined lines
  idirection(1:6)=1
  do j=1,6
    if ( l_line%n6line(iele, j) .lt. 0 ) idirection(j)=-1
  end do

  !# [2] ## Prepare the coordinates for 4 nodes of elements
!  write(*,*) "g_param%iflagspherical=",g_param%iflagspherical
!  stop
  do j=1,4
   n(j)=em_mesh%n4(iele,j) ! node index
   if ( g_param%iflagspherical .eq. 1) then   ! added on 2016.11.19
    elm_xyz(1:3,j)=em_mesh%xyzspherical(1:3,n(j)) ! [km] ! added on 2016.11.19
   else   ! added on 2016.11.19
    elm_xyz(1:3,j)=em_mesh%xyz(1:3,n(j)) ! [km]
   end if ! added on 2016.11.19
  end do

  !# [3] ## [x_mn]^T=L{x'34 x'14 x'42 x'23 x'31 x'12}=L[x'_lm]^T
  call calxmn(elm_xyz,xx)         ! see fem_util.f90
  call gradnodebasisfun(elm_xyz,gn,v) ! see fem_util.f90
  call LOCALMATRIX(elm_xyz,gn,xx,v,idirection,L0,dmu,sigma,S,S1) ! see forward_LHS.f90

  !# [4] ## Construct elemnt matrix, elm_k
   !elm_k(:,:)=S1(:,:)/dt ! for first-order Backward Euler
   elm_k(:,:)=S1(:,:)/2./dt ! second-order Backword Euler

  !# [5] ## source term
   S2(1:6)=0.d0
   if ( igroup(iele) .eq. 2 ) then ! when iele is included in ocean mesh
    do i=1,4     ! mean v*F in the tetrahedral element
!    F(1:3)= (/ 0.d0, 0.d0, Fxyz(3,n(i)) /)
     F(1:3)= Fxyz(1:3,n(i))
     vF(1:3,i) = outer(vxyz(1:3,n(i)),F(1:3)) ! [mm/s]*[nT]=[10^-12 V/m]
!     write(*,*) "vF(1:3)=",vF(1:3,i)
!     if (g_param%iflagspherical .eq. 1) then
!      call vectorspherical2xyz(em_mesh%lonlatalt(1:3,n(i)),vF(1:3,i),1)
!      write(*,*) "vF(1:3)=",vF(1:3,i)
!     end if
!     stop
    end do
    !#
    !#  mu*sigma*int[ w cdot N ]dv vF(1:4)
    !#  = mu*sigma*int[ {n_k*gn(:,l) - n_l*gn(:,k)}n_j]dv cdot vF(:,j)
    !#  = mu*sigma*{int(n_k*n_j)*[dv*gn(:,l) cdot vF(:)]
    !#             -int(n_l*n_j)*[dv*gn(:,k) cdot vF(:)]
    do i=1,6
     k=kl(i,1);l=kl(i,2)
     do j=1,4
     S2(i) = S2(i)+ dmu*sigma*( &
    &       intv(k,j,v)*inner(gn(:,l),vF(:,j)) &   ! for intv, see m_fem_util.f90
    &    -  intv(l,j,v)*inner(gn(:,k),vF(:,j)))*idirection(i) !
     end do
    end do
    end if

  !# [6] ## set right hand side vector, b_vec  ########################
   do i=1,6
    ii=l_line%n6line(iele,i)*idirection(i)
    b_vec(ii) = b_vec(ii) + S2(i)
    do j=1,6
     jj=l_line%n6line(iele,j)*idirection(j)
     !b_vec(ii)=b_vec(ii) + elm_k(i,j)*fsn1(jj)     ! for first-order Backward Euler
     b_vec(ii)=b_vec(ii) + elm_k(i,j)*(4.*fsn1(jj) - fsn0(jj)) ! for second order Backward Euler
    end do
   end do

   !# checkvalues
   !do i=1,6
   !ii=l_line%n6line(iele,i)*idirection(i)
   !if (elm_xyz(1,1) .ge. 0.d0 .and. igroup(iele) .eq. 1 .and. dabs(b_vec(ii)) .gt. 10 ) then
   ! write(*,*) "iele=",iele
   ! do j=1,4
   !  write(*,'(a,3g15.7)') "elm_xyz=",elm_xyz(1:3,j)
   ! end do
   ! do j=1,4
   !  write(*,'(a,3g15.7)') "vxyz=",vxyz(1:3,n(j))
   ! end do
   ! do j=1,4
   !  write(*,'(a,3g15.7)') "Fxyz=",Fxyz(1:3,n(j))
   ! end do
   ! do j=1,4
   !  write(*,'(a,3g15.7)') "v*F(1:3,j)=",vF(1:3,j)
   ! end do
   ! do j=1,6
   !  jj=l_line%n6line(iele,j)*idirection(j)
   !  write(*,'(a,g15.7,a,g15.7)') "b_vec(jj)",b_vec(jj),"fsn1(j)=",fsn1(jj)
   ! end do
   ! do j=1,6
   !  write(*,'(a,6g15.7)') "elm_k=",elm_k(j,1:6)
   ! end do
   !end if
   !goto 100
   !end do

   100 continue ! goto next element
end do ! element loop end

! apply boundry condition
do i=1, l_line%nline
 if ( A%line_bc(i) ) then
   b_vec(i) = A%Avalue_bc(i)
 else
   b_vec(i) = b_vec(i) - A%bcorrect(i)
 end if
end do

write(*,*) "### FORWARD_RHS END !! ###"
return
end subroutine forward_rhs
!
