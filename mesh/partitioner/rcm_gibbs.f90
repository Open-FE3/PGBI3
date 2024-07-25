!-----------------------------------------------------------------------
!   Copyright 1999, 2000, 2001, 2002, 2003 by the Research Organization
!   for Information Science & Technology (RIST)
!-----------------------------------------------------------------------

!C
!C***
!C*** RCM_GIBBS
!C***
!C
      subroutine RCM_GIBBS(IA, L, MA, M, NP, IFLAG, C, IER)
!C
!C    Reverse Cuthil-McKee Method for Symmetric Sparse Matrix
!C                                                                     
!C    on ENTRY:                                                        
!C      IA     ARRAY which contains ROW INDEX of NON-ZERO ELEMENTs 
!C      L      LEADIMG DIMENSION of IA.                  
!C      MA     acc. sum of # non-zero ELEMs in each COL. of the MATRIX
!C      M      COLUMN number of the MATRIX                     
!C      IFLAG  =0:pure RCM, =1:RCM with GIBBS
!C
!C    on RETURN:                                                       
!C      C      ARRAY which contains NEW NODE NUMBERs
!C      R      the ROW ACCOUNT in EACH ROW.                            
!C      IW     LEVEL NUMBER in EACH NODE.                              
!C      IER    ERROR CODE. if IER=0, NORMAL RETURN.
!C
!C    others: WORKING PARAMETERS.
!C
       implicit REAL*8 (A-H,O-Z)

       integer(kind=4), dimension(L)   :: IA
       integer(kind=4), dimension(0:M) :: MA
       integer(kind=4), dimension(NP)  :: C
       integer(kind=4)                 :: TOTAL

       integer(kind=4), dimension(:), allocatable ::                    &
     &                  R, IW, IR, IC, JC, KC
!C
!C +-------+
!C | INIT. |
!C +-------+
!C===
      allocate (R(M), IW(M), IR(M), IC(M), JC(M), KC(0:M))     

      R = 0
      IW= 0
      IR= 0
      IC= 0
      JC= 0
      KC= 0
!C===

       IER = 0
      TOTAL = 0
      LEVEL = 0
      
      IW= -1
       R=  0

!C
!C-- TRIANGULAR ELEMENT MUST BE NON-ZWRO
      do j= 1, M
        R(j) = MA(j) - MA(j-1) - 1
      enddo

      MIN = 10000
      IROW= 0
      do i= 1, M
        if (R(i) .lt. MIN) then
          MIN = R(i)
          IROW=   i
        endif
      enddo

!C
!C +---------------+
!C | without GIBBS |
!C +---------------+
!C===
      if (IFLAG .eq. 0) then
        IW(IROW) = LEVEL
           TOTAL = TOTAL + 1

   45   continue
        LEVEL = LEVEL + 1
        do j= 1, M
          if (IW(j) .eq. (LEVEL-1)) then
            do n= MA(j-1)+1,MA(j)
              if (IA(n) .ne. j) then
                if (IW(IA(n)) .lt. 0) then
                  IW(IA(n)) = LEVEL
                      TOTAL = TOTAL + 1
                endif
              endif
            enddo
          endif
        enddo

        if (TOTAL .lt. M) goto 45
!C
!C-- END. LEVELing
           LL  = LEVEL
        LEVEL  = 0
        C(IROW)= 1
        TOTAL  = 1
!C
   90   continue
        LEVEL= LEVEL + 1
           LN= 0

        do j= 1, M
          if (IW(j) .eq. LEVEL) then
               LN = LN + 1
            IR(LN)= j
            IC(LN)= R(j)
          endif
        enddo

        if (LN .eq. 0) then
          write (*,*) ' (SUBR. RCM) STOP AT ', LEVEL
          IER = 1
          deallocate (R,IW,IR,IC,JC,KC)     
          return
        endif

!C
!C-- ORDERING of MODES in THE SAME LEVEL
        do k= 1, LN
          TOTAL   = TOTAL + 1
          C(IR(k))= TOTAL
        enddo

        if (LEVEL .lt. LL) goto 90
!C===

!C
!C +------------+
!C | with GIBBS |
!C +------------+
!C===
       else
        JC= 0
        IC= 0
         C= 0

        call GIBBS (IA,L,MA,M,NP,R,IW,IROW,LEVEL,IR,IC)

         C=  0
        IW= -1

        IW(IROW)= LEVEL
          TOTAL = TOTAL + 1
             LL = 1  

        IC(1) = IROW
        JC(1) = IROW

  220   continue
             LN= LL
             LL= 0
          LEVEL= LEVEL + 1

          do k= 1, LN
            j= IC(k)
            do n= MA(j-1)+1, MA(j)
              if (IA(n) .ne. j) then
                if (IW(IA(n)) .lt. 0) then
                  IW(IA(n)) = LEVEL
                      TOTAL = TOTAL + 1
                         LL = LL + 1
                      IR(LL)= IA(n)
                   JC(TOTAL)= IA(n)
                endif
              endif
            enddo
          enddo

          if (LL .ne. 0) then
            do k= 1, LL
              IC(k) = IR(k)
            enddo
           else
            MIN = 100000
            do k= 1, M 
              if (IW(i) .lt. 0) then
                if (R(i) .lt. MIN) then
                   MIN = R(i)
                  IROW = i
                endif
              endif
            enddo
            
            IW(IROW) = LEVEL
               TOTAL = TOTAL + 1
                  LL = 1
               IC(1) = IROW
            JC(TOTAL)= IROW
          endif

          KC(LEVEL) = TOTAL
          if (TOTAL .lt. M) goto 220

!C
!C-- END. LEVELing
          TOTAL = 0
          do k= 1        , LEVEL
          do n= KC(k-1)+1, KC(k)
            TOTAL   = TOTAL + 1
            C(JC(n))= TOTAL  
          enddo
          enddo
        endif
!C===
      do i= 1, M
        C(i) = M + 1 - C(i)
      enddo
!C
      deallocate (R,IW,IR,IC,JC,KC)     

      return
      end

!C
!C***
!C*** GIBBS
!C***
      subroutine GIBBS (IA, L, MA, M, NP, R, IW, IROW, LEVEL, IR, IC)
!C
!C     selection of ROOT NODE with GIBBS METHOD for A SPARSE MATRIX.
!C
!C     on RETURN: 
!C       IROW   THE INITIAL NODE NUMBER SELECTED.
!C
      
      implicit REAL*8(A-H,O-Z)
      integer*4 R, ROOT, TOTAL, ROOTX
      dimension IA(L), R(M), MA(0:M), IW(M), IC(M), IR(M)

      IS   = 0
      ROOT = IROW
      ROOTX= IROW
      
   20 continue
        TOTAL= 0
        LEVEL= 0
        
        do i= 1, M
          IW(i) = - 1
        enddo

        IW(ROOT)= LEVEL
           TOTAL= TOTAL + 1
              LL= 1
           IC(1)= ROOT
!C
   45   continue
             LN= LL
             LL= 0
          LEVEL= LEVEL + 1
                    
          do k= 1, LN
            j= IC(k)
            do n= MA(j-1)+1, MA(j)
              if (IA(n) .ne. j) then
                if (IW(IA(n)) .lt. 0) then
                  IW(IA(n)) = LEVEL
                      TOTAL = TOTAL + 1
                         LL = LL + 1
                      IR(LL)= IA(n)
                endif
              endif
            enddo
          enddo
          
          if (LL .ne. 0) then
            do k= 1, LL
              IC(k)= IR(k)
            enddo
           else
            MIN = 100000
            do i= 1, M
              if (IW(i) .lt. 0) then
                if (R(i) .lt. MIN) then
                  MIN = R(i)
                  IROW=   i
                endif
              endif
            enddo
           
            IW(IROW)= LEVEL
               TOTAL= TOTAL + 1
               IC(1)= IROW
                  LL= 1
          endif
                             
          if (TOTAL .lt. M) goto 45
!C
!C-- END. LEVELing
          MIN = 100000
          do i= 1, M
            if (IW(i) .eq. LEVEL) then
              if (R(i) .lt. MIN) then
                MIN= R(i)
                 IX=   i
              endif
            endif
          enddo

          if (IS .ne. 1) then
            IS    = 1
            ROOT  = IX
            ILEVEL= LEVEL
            goto 20
          endif

         IROW = ROOTX
         IER = 0
      
      return     
      end






