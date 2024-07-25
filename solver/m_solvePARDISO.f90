! Coded by T.M. on May 13, 2016
! based on pardiso_sym_f90.f90
!PROGRAM pardiso_sym_f90
INCLUDE 'mkl_pardiso.f90'
module solvePARDISO
USE mkl_pardiso
IMPLICIT NONE
!INTEGER, PARAMETER :: dp = KIND(1.0D0)
!.. Internal solver memory pointer
!################################################## type PARDISO_PARAM
type PARDISO_PARAM
TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: pt(:)
!.. All other variables
INTEGER maxfct, mnum, mtype, n, nrhs, error, msglvl, nnz
INTEGER error1
INTEGER, ALLOCATABLE :: iparm( : )
INTEGER, ALLOCATABLE :: ia( : )
INTEGER, ALLOCATABLE :: ja( : )
REAL(8), ALLOCATABLE :: amat( : )
REAL(8), ALLOCATABLE :: b( : )
REAL(8), ALLOCATABLE :: x( : )
!INTEGER i, j,idum(1)
INTEGER(4) :: idum(1)
REAL(8) ddum(1)
end type

contains
!#################################################################  PARDISOphase1
subroutine PARDISOphase1(A,B)
use iccg_var_takuto
type(global_matrix),intent(inout) :: A ! A is deallocated after pardiso
type(PARDISO_PARAM),intent(inout) :: B
integer(4) :: i,j,phase
!
write(*,*) "### PARDISOphase1 start!!###"
!.. Fill all arrays containing matrix data.
B%n = A%doftot               ! number of equations
B%nnz = A%iau_tot + A%doftot ! upper triangle + diagonal
B%nrhs = 1                 ! number of right hand side vector
B%maxfct = 1
B%mnum = 1
!################################################# for debug
!n = 8
!nnz = 18
!nrhs = 1
!maxfct = 1
!mnum = 1
!deallocate(A%AU,A%IAU,A%INU,A%D)
!allocate(A%AU(nnz-n),A%IAU(nnz),A%INU(0:n),A%D(n))
!A%D=(/7.d0,-4.d0,1.d0,7.d0,5.d0,-1.d0,11.d0,5.d0/)
!A%AU=(/1.d0,2.d0,7.d0,8.d0,2.d0,  5.d0,9.d0,1.d0,5.d0,5.d0/)
!A%INU=(/0,3,5,6,7,9,10,10,10/)
!A%IAU=(/3,6,7,3,5,8,7,6,7,8/)
!#################################################
ALLOCATE(B%ia(B%n + 1))
ALLOCATE(B%ja(B%nnz))
ALLOCATE(B%amat(B%nnz))
ALLOCATE(B%b(B%n))
ALLOCATE(B%x(B%n))
  B%IA(1)=1
  DO I = 2,B%n+1
   B%IA(I)=B%IA(I-1)+(A%INU(I-1)-A%INU(I-2))+1 !the first of Ith row
   B%JA(B%IA(I-1))=I-1 ! diagonal
   B%AMAT(B%IA(I-1))=A%D(I-1)
   DO J =1,(A%INU(I-1) - A%INU(I-2)) ! upper triangle of I-1 th row
   B%JA  (B%IA(I-1)+J) = A%IAU(A%INU(I-2)+J) !upper triangle
   B%AMAT(B%IA(I-1)+J) =  A%AU(A%INU(I-2)+J)
   END DO
  END DO
!do i=1,8
!write(*,'(a,3i7,2e15.7)') "I,INU,IAU,A",i,A%INU(i),A%IAU(i),A%AU(i)
!end do
deallocate (A%IAU, A%IAL, A%INU, A%INL, A%AU, A%AL, A%D)

!..
!.. Set up PARDISO control parameter
!..
ALLOCATE(B%iparm(64))

DO i = 1, 64
   B%iparm(i) = 0
END DO

B%iparm(1) = 1 ! no solver default
B%iparm(2) = 2 ! fill-in reordering from METIS
B%iparm(4) = 0 ! no iterative-direct algorithm
B%iparm(5) = 0 ! no user fill-in reducing permutation
B%iparm(6) = 0 ! =0: solution is stored in x, while b is not changed
B%iparm(8) = 2 ! numbers of iterative refinement steps
B%iparm(10) = 13 ! perturb the pivot elements with 1E-13
B%iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
B%iparm(13) = 1 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
! iparm(13) should be 1, if iparm(13)=0, the accuracy of solution dramatically fall down
! on May 13, 2016
B%iparm(14) = 0 ! Output: number of perturbed pivots
B%iparm(18) = -1 ! Output: number of nonzeros in the factor LU
B%iparm(19) = -1 ! Output: Mflops for LU factorization
B%iparm(20) = 0 ! Output: Numbers of CG Iterations

B%error  = 0 ! initialize error flag
B%msglvl = 1 ! print statistical information
!mtype  = 6 ! complex and symmetric
B%mtype = -2 ! real and symmetric indefinite

!.. Initialize the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.

ALLOCATE (B%pt(64))
DO i = 1, 64
   B%pt(i)%DUMMY =  0
END DO

!.. Reordering and Symbolic Factorization, This step also allocates
! all memory that is necessary for the factorization

phase = 11 ! only reordering and symbolic factorization

CALL pardiso (B%pt, B%maxfct, B%mnum, B%mtype, phase, B%n, B%AMAT, B%ia, B%ja, &
              B%idum, B%nrhs, B%iparm, B%msglvl, B%ddum, B%ddum, B%error)

WRITE(*,*) 'Reordering completed ... '
IF (B%error /= 0) THEN
   WRITE(*,*) 'The following ERROR was detected: ', B%error
!   GOTO 1000
call PARDISOphase4(B)
END IF
WRITE(*,*) 'Number of nonzeros in factors = ',B%iparm(18)
WRITE(*,*) 'Number of factorization MFLOPS = ',B%iparm(19)

return
end subroutine PARDISOphase1
!
!########################################################### phase2
subroutine PARDISOphase2(B)
type(PARDISO_PARAM),intent(inout) :: B
integer(4) :: phase
write(*,*) "### PARDISOphase2 start!! ###"
!.. Factorization.
phase = 22 ! only factorization
CALL pardiso (B%pt, B%maxfct, B%mnum, B%mtype, phase, B%n, B%AMAT, B%ia, B%ja, &
              B%idum, B%nrhs, B%iparm, B%msglvl, B%ddum, B%ddum, B%error)
WRITE(*,*) 'Factorization completed ... '
IF (B%error /= 0) THEN
   WRITE(*,*) 'The following ERROR was detected: ', B%error
!   GOTO 1000
call PARDISOphase4(B)
ENDIF
write(*,*) "### PARDISOphase2 END!! ###"
return
end subroutine PARDISOphase2

!########################################################### phase3
subroutine PARDISOphase3(B,b_vec,nline,xout)
type(PARDISO_PARAM),intent(inout) :: B
integer(4),intent(in) :: nline
real(8),intent(inout) :: b_vec(nline) !when iparm(6)=0, solution is in xout
real(8),intent(out) :: xout(nline)
integer(4) :: i,phase
!.. Back substitution and iterative refinement
B%iparm(8) = 2 ! max numbers of iterative refinement steps
phase = 33 ! only solving
B%b(:)=b_vec(:)
B%msglvl = 0 ! not print statistical info
!############################## for debug
!deallocate(b,x)
!allocate(b(n),x(n))
!DO i = 1, n
!   b(i) = 1.d0
!END DO
!##############################
CALL pardiso (B%pt, B%maxfct, B%mnum, B%mtype, phase, B%n, B%AMAT, B%ia, B%ja, &
              B%idum, B%nrhs, B%iparm, B%msglvl, B%b, B%x, B%error)
!WRITE(*,*) 'Solve completed ... '
IF (B%error /= 0) THEN
   WRITE(*,*) 'The following ERROR was detected: ', B%error
!   GOTO 1000
   CALL PARDISOphase4(B)
ENDIF
write(*,*) "### PARDISOphase3 END!! ###"
!########################################### for debug
!WRITE(*,*) 'The solution of the system is '
DO i = 1, B%n
!   WRITE(*,*) ' x(',i,') = ', x(i)
   xout(i) = B%x(i)
END DO
return
!######

!1000 CONTINUE
return
end subroutine PARDISOphase3
!
!########################################################### phase4
subroutine PARDISOphase4(B)
type(PARDISO_PARAM),intent(inout) :: B
integer(4) :: phase,error1
!.. Termination and release of memory
phase = -1 ! release internal memory
CALL pardiso (B%pt, B%maxfct, B%mnum, B%mtype, phase, B%n, B%ddum, B%idum, B%idum, &
              B%idum, B%nrhs, B%iparm, B%msglvl, B%ddum, B%ddum, error1)

IF (ALLOCATED(B%ia))      DEALLOCATE(B%ia)
IF (ALLOCATED(B%ja))      DEALLOCATE(B%ja)
IF (ALLOCATED(B%AMAT))    DEALLOCATE(B%AMAT)
IF (ALLOCATED(B%b))       DEALLOCATE(B%b)
IF (ALLOCATED(B%x))       DEALLOCATE(B%x)
IF (ALLOCATED(B%iparm))   DEALLOCATE(B%iparm)

IF (error1 /= 0) THEN
   WRITE(*,*) 'The following ERROR on release stage was detected: ', error1
   STOP 1
ENDIF

IF (B%error /= 0) STOP 1
!END PROGRAM pardiso_sym_f90

return
end subroutine PARDISOphase4

end module solvePARDISO
