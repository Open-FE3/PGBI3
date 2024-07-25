!coded on Jan. 20, 2016 by T.M.
!based on sort.for in Numerical recipes in fortran
!program main
!implicit none
!integer(4) :: n4(1:4)=(/  5,   10,  1,  20 /)
!real(8)   ::  r4(1:4)=(/3.0,  2.0,1.0, 4.0 /)
!call sort_index(4,n4,r4)
!write(*,*) n4(1:4)
!write(*,*) r4(1:4)
!end program main
!########################################################
! modified from SORT.for
SUBROUTINE SORT_INDEX(N,NA,RA)
      implicit real(selected_real_kind(8))(a-h,o-z)
	INTEGER(4) :: N
	INTEGER(4),DIMENSION(N) :: NA
      REAL(8)   ,DIMENSION(N) :: RA
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
	    IIA=NA(L) ! inserted
        ELSE
          RRA=RA(IR)
	    IIA=NA(IR) ! inserted
          RA(IR)=RA(1)
	    NA(IR)=NA(1) ! inserted
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
	      NA(1)=IIA ! inserted
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
		NA(I)=NA(J) ! inserted
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        NA(I)=IIA ! inserted
      GO TO 10
      END

