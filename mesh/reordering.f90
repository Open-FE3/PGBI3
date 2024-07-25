! Coppied from partitioner/local_data.f90 on 2015.06.05
!  This file includes
! subroutine REORDERING, RCM_GIBBS, GIBBS, and mSORT
	  subroutine REORDERING(N, IELMTOT)
        implicit real(selected_real_kind(8))(a-h,o-z)
        integer(kind=kint), dimension(:),  allocatable :: INL  , INU
        integer(kind=kint), dimension(:,:),allocatable :: IAL, IAU
        integer(kind=kint), dimension(:),  allocatable :: OLDtoNEW_gibbs
        integer(kind=kint), dimension(:),  allocatable ::  ITEM_gibbs
        integer(kind=kint), dimension(:),  allocatable :: STACK_gibbs
        integer(kind=kint), dimension(100) :: NCOL1, NCOL2
        integer(kind=kint), dimension(:),  allocatable :: REGION
        integer(kind=kint) :: REGIONtot
        integer(kind=kint), dimension(:),  allocatable :: REGIONstack
        integer(kind=kint), dimension(:),  allocatable :: REGIONstackWK
        integer(kind=kint), dimension(:),  allocatable :: REGIONitem
!C +------------------+
!C | connected region |
!C +------------------+
      allocate (REGION(N))
      REGION= 0
      icouT= 0
OUTER:do iterO= 1, N
        iflagO= 0
INNER:  do iter= 1, N
            icou= 0
          do icel= 1, IELMTOT        
            iflag= 0
            do k= 1, NODELM(IELMLOCAL(icel))
              in= ICELNOD(icel,k)
              if (iflagO.eq.0.and.REGION(in).eq.0) then
                REGION(in)= iterO
                iflagO    = 1
                icouT     = icouT + 1
                icou      = 1
              endif
              if (REGION(in).eq.iterO) then
                iflag= 1
                exit
              endif
            enddo
            if (iflag.eq.1) then
              do k= 1, NODELM(IELMLOCAL(icel))
                in = ICELNOD(icel,k)
                if (REGION(in).eq.0) then
                  icou = icou  + 1
                  icouT= icouT + 1
                  REGION(in)= iterO
                  if (icouT.eq.N) exit OUTER
                endif
              enddo
            endif
          enddo
          if (icou.eq.0) exit INNER
        enddo INNER
      enddo   OUTER

      REGIONtot= iterO
      allocate (REGIONstack  (0:REGIONtot), REGIONitem(N))
      allocate (REGIONstackWK(0:REGIONtot))
                REGIONstack  = 0
                REGIONstackWK= 0
                REGIONitem   = 0
      do i= 1, N
        ic= REGION(i)
        REGIONstack(ic)= REGIONstack(ic) + 1
      enddo
      do ic= 1, REGIONtot
        REGIONstack(ic)= REGIONstack(ic-1) + REGIONstack(ic)
      enddo
      do i= 1, N
        ic= REGION(i)
        REGIONstackWK(ic)= REGIONstackWK(ic) + 1
        k = REGIONstack(ic-1) + REGIONstackWK(ic)
        REGIONitem(k)= i
      enddo
!C===

!C
!C +---------------------+
!C | matrix connectivity |
!C +---------------------+
!C===
        NU= 50
        NL= 50

  125   continue
        allocate (INL(N2), IAL(N2,NL))
        allocate (INU(N2), IAU(N2,NU))

        INL= 0
        IAL= 0
        INU= 0
        IAU= 0

        do icel= 1, IELMTOT        
        do k1  = 1, NODELM(IELMLOCAL(icel))
        do k2  = 1, NODELM(IELMLOCAL(icel))
          if (k1.ne.k2) then
            in1= ICELNOD(icel,k1)
            in2= ICELNOD(icel,k2)
            if (in1.le.N2.and.in2.le.N2) then
            if (in1.gt.in2) then
              do kk= 1, INL(in1)
                if (in2.eq.IAL(in1,kk)) goto 130
              enddo
              icou= INL(in1) + 1
              if (icou.gt.NL) then
                deallocate (INL, INU, IAL, IAU)
                NL= NL + 2
                goto 125
              endif
              IAL(in1,icou)= in2
              INL(in1     )= icou
             else
              do kk= 1, INU(in1)
                if (in2.eq.IAU(in1,kk)) goto 130
              enddo
              icou= INU(in1) + 1
              if (icou.gt.NU) then
                deallocate (INL, INU, IAL, IAU)
                NU= NU + 2
                goto 125
              endif
              IAU(in1,icou)= in2
              INU(in1     )= icou
            endif
            endif
          endif          
  130     continue
        enddo 
        enddo
        enddo

        do in= 1, N2
          NN= INL(in)
          do k= 1, NN
            NCOL1(k)= IAL(in,k)
          enddo
            call mSORT (NCOL1, NCOL2, NN)
          do k= NN, 1, -1
            IAL(in,NN-k+1)= NCOL1(NCOL2(k))
          enddo        

          NN= INU(in)
          do k= 1, NN
            NCOL1(k)= IAU(in,k)
          enddo
            call mSORT (NCOL1, NCOL2, NN)
          do k= NN, 1, -1
            IAU(in,NN-k+1)= NCOL1(NCOL2(k))
          enddo        
        enddo
!C===
        NUmax= 0
        NLmax= 0
        do i= 1, N2
          NUmax= max (NUmax, INU(i))
          NLmax= max (NLmax, INL(i))
        enddo

!C
!C +----------------------+
!C | RCM/Gibbs reordering |
!C +----------------------+
!C===
        allocate (STACK_gibbs(0:N))
        NLU        = 0
        STACK_gibbs= 0
        do i= 1, N2
          NLU           = NLU              + INL(i) + 1 + INU(i)
          STACK_gibbs(i)= STACK_gibbs(i-1) + INL(i) + 1 + INU(i)
        enddo

        allocate (ITEM_gibbs(NLU))
        do i= 1, N2
          do k= 1, INL(i)
                       in = STACK_gibbs(i-1) + k
            ITEM_gibbs(in)= IAL(i,k)
          enddo
          ITEM_gibbs(STACK_gibbs(i-1)+INL(i)+1)= i
          do k= 1, INU(i)
                       in = STACK_gibbs(i-1) + INL(i) + 1 + k
            ITEM_gibbs(in)= IAU(i,k)
          enddo
        enddo
!C===
        deallocate (INL, INU, IAL, IAU)
          allocate (OLDtoNEW_gibbs(INODTOT))

        call RCM_GIBBS (ITEM_gibbs, NLU, STACK_gibbs, N2, N2, 1,        &
     &                  OLDtoNEW_gibbs, IER)

        deallocate (STACK_gibbs, ITEM_gibbs)

        do i= 1, N2
          in= OLDtoNEW_gibbs(i)
          OLDtoNEW(i )= in
          NEWtoOLD(in)= i
        enddo

        do i= N2+1, INODTOT
          OLDtoNEW(i)= i
          NEWtoOLD(i)= i
        enddo

        do icel= 1, IELMTOT
          do k= 1, NODELM(IELMLOCAL(icel))
            in= OLDtoNEW(ICELNOD(icel,k))
            ICELNOD(icel,k)= in
          enddo
        enddo

        end subroutine REORDERING
!##########################################################
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
!#######################################################
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








