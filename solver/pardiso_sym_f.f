********************************************************************************
*                              INTEL CONFIDENTIAL
*   Copyright(C) 2004-2008 Intel Corporation. All Rights Reserved.
*   The source code contained  or  described herein and all documents related to
*   the source code ("Material") are owned by Intel Corporation or its suppliers
*   or licensors.  Title to the  Material remains with  Intel Corporation or its
*   suppliers and licensors. The Material contains trade secrets and proprietary
*   and  confidential  information of  Intel or its suppliers and licensors. The
*   Material  is  protected  by  worldwide  copyright  and trade secret laws and
*   treaty  provisions. No part of the Material may be used, copied, reproduced,
*   modified, published, uploaded, posted, transmitted, distributed or disclosed
*   in any way without Intels prior express written permission.
*   No license  under any  patent, copyright, trade secret or other intellectual
*   property right is granted to or conferred upon you by disclosure or delivery
*   of the Materials,  either expressly, by implication, inducement, estoppel or
*   otherwise.  Any  license  under  such  intellectual property  rights must be
*   express and approved by Intel in writing.
*
********************************************************************************
*   Content : MKL PARDISO Fortran-77 example
*
********************************************************************************
C----------------------------------------------------------------------
C Example program to show the use of the "PARDISO" routine
C for symmetric linear systems
C---------------------------------------------------------------------
C This program can be downloaded from the following site:
C www.pardiso-project.org
C
C (C) Olaf Schenk, Department of Computer Science,
C University of Basel, Switzerland.
C Email: olaf.schenk@unibas.ch
C
C---------------------------------------------------------------------
C      PROGRAM pardiso_sym

      IMPLICIT NONE
      include 'mkl_pardiso.f77'
C.. Internal solver memory pointer for 64-bit architectures
C.. INTEGER*8 pt(64)
C.. Internal solver memory pointer for 32-bit architectures
C.. INTEGER*4 pt(64)
C.. This is OK in both cases
      INTEGER*8 pt(64)
C.. All other variables
      INTEGER maxfct, mnum, mtype, phase, n, nrhs, error, msglvl
      INTEGER iparm(64)
      INTEGER ia(9)
      INTEGER ja(18)
      REAL*8 a(18)
      REAL*8 b(8)
      REAL*8 x(8)
      INTEGER i, idum
      REAL*8 waltime1, waltime2, ddum
C.. Fill all arrays containing matrix data.
      DATA n /8/, nrhs /1/, maxfct /1/, mnum /1/
      DATA ia /1,5,8,10,12,15,17,18,19/
      DATA ja
     1 /1,  3,    6,7,
     2    2,3,  5,
     3      3,        8,
     4        4,    7,
     5          5,6,7,
     6            6,  8,
     7              7,
     8                8/
      DATA a
     1 /7.d0,     1.d0,          2.d0,7.d0,
     2       -4.d0,8.d0,     2.d0,
     3            1.d0,                    5.d0,
     4                 7.d0,     9.d0,
     5                      5.d0,1.d0,5.d0,
     6                           -1.d0,     5.d0,
     7                                11.d0,
     8                                     5.d0/
C..
C.. Set up PARDISO control parameter
C..
      do i = 1, 64
         iparm(i) = 0
      end do
      iparm(1) = 1 ! no solver default
      iparm(2) = 2 ! fill-in reordering from METIS
      iparm(3) = 1 ! numbers of processors
      iparm(4) = 0 ! no iterative-direct algorithm
      iparm(5) = 0 ! no user fill-in reducing permutation
      iparm(6) = 0 ! =0 solution on the first n compoments of x
      iparm(7) = 0 ! not in use
      iparm(8) = 9 ! numbers of iterative refinement steps
      iparm(9) = 0 ! not in use
      iparm(10) = 13 ! perturbe the pivot elements with 1E-13
      iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
      iparm(12) = 0 ! not in use
      iparm(13) = 0 ! not in use
      iparm(14) = 0 ! Output: number of perturbed pivots
      iparm(15) = 0 ! not in use
      iparm(16) = 0 ! not in use
      iparm(17) = 0 ! not in use
      iparm(18) = -1 ! Output: number of nonzeros in the factor LU
      iparm(19) = -1 ! Output: Mflops for LU factorization
      iparm(20) = 0 ! Output: Numbers of CG Iterations
      error = 0 ! initialize error flag
      msglvl = 1 ! print statistical information
      mtype = -2 ! symmetric, indefinite
C.. Initiliaze the internal solver memory pointer. This is only
C necessary for the FIRST call of the PARDISO solver.
      do i = 1, 64
         pt(i) = 0
      end do
C.. Reordering and Symbolic Factorization, This step also allocates
C all memory that is necessary for the factorization
      phase = 11 ! only reordering and symbolic factorization
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,
     1 idum, nrhs, iparm, msglvl, ddum, ddum, error)
      WRITE(*,*) 'Reordering completed ... '
      IF (error .NE. 0) THEN
         WRITE(*,*) 'The following ERROR was detected: ', error
         STOP
      END IF
      WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
      WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)
C.. Factorization.
      phase = 22 ! only factorization
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,
     1 idum, nrhs, iparm, msglvl, ddum, ddum, error)
      WRITE(*,*) 'Factorization completed ... '
      IF (error .NE. 0) THEN
         WRITE(*,*) 'The following ERROR was detected: ', error
         STOP
      ENDIF
C.. Back substitution and iterative refinement
      iparm(8) = 2 ! max numbers of iterative refinement steps
      phase = 33 ! only factorization
      do i = 1, n
         b(i) = 1.d0
      end do
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,
     1 idum, nrhs, iparm, msglvl, b, x, error)
      WRITE(*,*) 'Solve completed ... '
      WRITE(*,*) 'The solution of the system is '
      DO i = 1, n
         WRITE(*,*) ' x(',i,') = ', x(i)
      END DO
C.. Termination and release of memory
      phase = -1 ! release internal memory
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum,
     1 idum, nrhs, iparm, msglvl, ddum, ddum, error)
      END
