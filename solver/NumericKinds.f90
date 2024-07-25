module NumericKinds
  implicit none

  ! KINDs for different types of REALs
  ! .............................
  integer, parameter :: &
       R_SP = selected_real_kind(  6,   37 ), &
       R_DP = selected_real_kind( 15,  307 ), &
       R_QP = selected_real_kind( 33, 4931 )
  ! Alias for precision that we use in the program (change this to any of the
  ! values 'R_SP', 'R_DP', or 'R_QP', to switch to another precision globally).
  integer, parameter :: RK = R_DP ! if changing this, also change RK_FMT

  ! KINDs for different types of INTEGERs
  ! .............................
  integer, parameter :: &
       I1B = selected_int_kind(2), & ! max = 127
       I2B = selected_int_kind(4), & ! max ~ 3.28x10^4
       I3B = selected_int_kind(9), & ! max ~ 2.15x10^9
       I4B = selected_int_kind(18)   ! max ~ 9.22x10^18
  ! Alias for integer-precision (analogue role to RK above).
  integer, parameter :: IK = I3B

  ! Edit-descriptors for real-values
  character(len=*), parameter :: R_SP_FMT = "f0.6", &
       R_DP_FMT = "f0.15", R_QP_FMT = "f0.33"
  ! Alias for output-precision to use in the program (keep this in sync with RK)
  character(len=*), parameter :: RK_FMT = R_DP_FMT
end module NumericKinds

