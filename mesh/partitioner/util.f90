!-----------------------------------------------------------------------
!   Copyright 1999, 2000, 2001, 2002, 2003 by the Research Organization
!   for Information Science & Technology (RIST)
!-----------------------------------------------------------------------

!C
!C***
!C*** mSORT
!C***
!C
      subroutine mSORT (STEM,INUM,NN)
      integer*4  STEM(NN), INUM(NN)

      do ii = 1,NN
        INUM(ii)= ii
      enddo

      do ii= 1,NN-1
      do jj= 1,NN-ii
        if (STEM(INUM(jj)) .lt. STEM(INUM(jj+1))) then
          ITEM      = INUM(jj+1)
          INUM(jj+1)= INUM(jj)
          INUM(jj)  = ITEM
        endif
      enddo
      enddo

      return
      end

!C
!C***
!C*** DATA_COMPRESS
!C***
!C
      subroutine DATA_COMPRESS (I1, I2, I_INN, I_OUT)
      integer  I_INN(0:I1), I_OUT(0:I1)

      ITOT= I1
      do i= 0, I1
        I_OUT(i)= I_INN(i)
      enddo

      I2= ITOT

      return
      end




