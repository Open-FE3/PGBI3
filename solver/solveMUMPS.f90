! Modified by T.M. on April 18, 2016
! to solve the finite element linear problem
!
!############################################
!  This file is part of MUMPS 5.0.1, released
!  on Thu Jul 23 17:08:29 UTC 2015
!
      subroutine solveMUMPS(doftot,A,rhs,X,ip)
	use iccg_var_takuto
      IMPLICIT NONE
!==========================
type(global_matrix),intent(inout) :: A ! A is deallocated before mumps to increase memory
integer(4),intent(in)  :: doftot
complex(8),intent(in)  :: rhs(doftot)
complex(8),intent(out) :: X(doftot) ! solution vector
integer(4),intent(out) :: ip
!==========================
      INCLUDE 'mpif.h'
      INCLUDE 'zmumps_struc.h'
      TYPE (ZMUMPS_STRUC) mumps_par
      INTEGER IERR, I,J
	integer(4),allocatable,dimension(:) :: IRHS_SPARSE_work
	complex(8),allocatable,dimension(:) :: RHS_SPARSE_work
	integer(4) :: icount
      CALL MPI_INIT(IERR)
! Define a communicator for the package.
      mumps_par%COMM = MPI_COMM_WORLD
!  Initialize an instance of the package
!  for L U factorization (sym = 0, with working host)
      mumps_par%JOB = -1
      mumps_par%SYM = 2 ! symmetric 
      mumps_par%PAR = 1
      CALL ZMUMPS(mumps_par)
      IF (mumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
     &            "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), &
     &            "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
       GOTO 500
      END IF
!  Define problem on the host (processor 0)
      ip=mumps_par%MYID
      IF ( mumps_par%MYID .eq. 0 ) THEN
        mumps_par%N=doftot
        mumps_par%NZ=A%iau_tot + doftot ! upper triangle + diagonal
        ALLOCATE( mumps_par%IRN ( mumps_par%NZ ) )
        ALLOCATE( mumps_par%JCN ( mumps_par%NZ ) )
        ALLOCATE( mumps_par%A   ( mumps_par%NZ ) )
        ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )
	  ! should be allocated when the solution should be centralized
	  DO I = 1,mumps_par%N
          mumps_par%IRN(I)=I
	    mumps_par%JCN(I)=I
	    mumps_par%A(I)=A%D(I)
!	    mumps_par%RHS(I)=rhs(I)
	  END DO

	  ! gen sparse rhs added on May 11, 2016 ###########
	  allocate( IRHS_SPARSE_work(mumps_par%N))
	  allocate( RHS_SPARSE_work(mumps_par%N))
	  icount=0
	  DO I = 1,mumps_par%N
         if (cdabs(rhs(I)) .ge. 1.d-30) then
	    icount=icount+1
	    IRHS_SPARSE_work(icount)=I
	    RHS_SPARSE_work(icount)=rhs(I)
	   end if
	  END DO
	  allocate( mumps_par%IRHS_SPARSE(icount))
	  allocate(  mumps_par%RHS_SPARSE(icount))
        mumps_par%IRHS_SPARSE(1:icount)=IRHS_SPARSE_work(1:icount)
        mumps_par%RHS_SPARSE(1:icount) = RHS_SPARSE_work(1:icount)
	  deallocate(IRHS_SPARSE_work,RHS_SPARSE_work)
	  allocate( mumps_par%IRHS_PTR(2) )
	  mumps_par%NZ_RHS=icount
	  mumps_par%NRHS=1 ! number of vector
	  mumps_par%IRHS_PTR(1:2)=(/1,icount+1/)
	  mumps_par%ICNTL(20)=1 ! exploit the sparsity of right hand side
        ! sparse rhs id generated  ##########################

        DO I = 1, doftot
	    DO J=A%INU(I-1)+1,A%INU(I)
           mumps_par%IRN(J+doftot)=I
	     mumps_par%JCN(J+doftot)=A%IAU(J)
	     mumps_par%A  (J+doftot)=A%AU(J)
	    END DO
        END DO
      END IF
	deallocate (A%IAU, A%IAL, A%INU, A%INL, A%AU, A%AL, A%D)
!  Call package for solution
      mumps_par%JOB = 6
      CALL ZMUMPS(mumps_par)
      IF (mumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
     &            "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1),&
     &            "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
       GOTO 500
      END IF
!  Solution has been assembled on the host
!      IF ( mumps_par%MYID .eq. 0 ) THEN
!        WRITE( 6, * ) ' Solution is ',(mumps_par%RHS(I),I=1,mumps_par%N)
!      END IF
!  Deallocate user data
      IF ( mumps_par%MYID .eq. 0 )THEN
        DO I=1,doftot
	   X(I)=mumps_par%RHS(I)
 	  END DO
        DEALLOCATE( mumps_par%IRN )
        DEALLOCATE( mumps_par%JCN )
        DEALLOCATE( mumps_par%A   )
        DEALLOCATE( mumps_par%RHS )
      END IF
      call MPI_BCAST (X(1) , doftot, MPI_DOUBLE_COMPLEX, 0, mumps_par%COMM, ierr)
!  Destroy the instance (deallocate internal data structures)
      mumps_par%JOB = -2
      CALL ZMUMPS(mumps_par)
      IF (mumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",&
     &            "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1),&
     &            "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2)
       GOTO 500
      END IF
 500  CALL MPI_FINALIZE(IERR)
      RETURN
      STOP
      END

