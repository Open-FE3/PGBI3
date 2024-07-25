!-----------------------------------------------------------------------
!   Copyright 1999, 2000, 2001, 2002, 2003 by the Research Organization
!   for Information Science & Technology (RIST)
!-----------------------------------------------------------------------

      subroutine LOCAL_DATA
      use analyzer

      character*80 LINE
      integer(kind=kint), dimension(:)  , allocatable :: IWKM
      integer(kind=kint), dimension(:)  , allocatable :: OLDtoNEW
      integer(kind=kint), dimension(:)  , allocatable :: NEWtoOLD
      integer(kind=kint), dimension(:)  , allocatable :: IMASK1, IMASK2
      integer(kind=kint), dimension(:)  , allocatable :: IELMLOCAL

!C
!C-- init.
      allocate (IWKM  (NP))

      N2= max (IELMTOT, N)
      allocate (IMASK1(N2))
      allocate (IMASK2(N2))

      allocate (IELMLOCAL(IELMTOT))

      INODTOTG= N       
      IELMTOTG= IELMTOT

      write (*,'(//,a)') ' Do you need Re-Ordering (RCM/Gibbs) ?'
      write (*,'( /,a)') ' 0.NO, 1.YES'
      write (*,'(a)')    ' Please TYPE 0 or 1 !!'
      write (*,'(/,">>>")')
      read  (*,*) iopt

      do ip= 1, NP

!C
!C +--------------------------+
!C | read INITIAL LOCAL files |
!C +--------------------------+
!C===
      open (11,file=WORKFIL(ip), status='unknown',form='unformatted')
        rewind (11)
          read (11) ID
          read (11) INODTOT
          read (11) N2
          read (11) N3         
          read (11) (IWKM(i),i=1,N3)

          do i= 1, INODTOT
            read (11) INODLOCAL(i)
          enddo

          STACK_IMPORT(0,ip)= 0
          read (11) (STACK_IMPORT(k,ip), k=1, N3)
          do is= 1, STACK_IMPORT(N3,ip)
            read (11) NOD_IMPORT(is)          
          enddo

          STACK_EXPORT(0,ip)= 0
          allocate (NOD_EXPORT(STACK_EXPORT(N3,ip)))
          read (11) (STACK_EXPORT(k,ip), k=1, N3)
          do is= 1, STACK_EXPORT(N3,ip)
            read (11) NOD_EXPORT(is)          
          enddo

          read (11) IELMTOT
          do i= 1, IELMTOT
            read (11) IELMLOCAL(i),                                     &
     &               (ICELNOD(i,k), k=1,NODELM(IELMLOCAL(i)))
          enddo
        close (11)
!C===

!C
!C +-----------------+
!C | LOCAL NUMBERING |
!C +-----------------+
!C===
          do i= 1, INODTOTG
            IMASK1(i)= 0
          enddo

          do i= 1, IELMTOTG
            IMASK2(i)= 0
          enddo

          do i= 1, INODTOT
            in= INODLOCAL(i)
            IMASK1(in)= i 
          enddo

          do i= 1, IELMTOT
            in= IELMLOCAL(i)
            IMASK2(in)= i 
          enddo

          ELMGRPSTACK(0)= 0
          do ig= 1, ELMGRPTOT
            icou= 0
            do is= ELMGRPSTACKG(ig-1)+1,ELMGRPSTACKG(ig)
              i= ELMGRPITEMG(is)
              if (IMASK2(i).ne.0) then
                icou= icou + 1
                  in= ELMGRPSTACK(ig-1) + icou
                ELMGRPITEM(in)= IMASK2(i)
              endif
            enddo
            ELMGRPSTACK(ig)= ELMGRPSTACK(ig-1) + icou
          enddo

          NODGRPSTACK(0)= 0
          do ig= 1, NODGRPTOT
            icou= 0
            do is= NODGRPSTACKG(ig-1)+1,NODGRPSTACKG(ig)
              i= NODGRPITEMG(is)
              if (IMASK1(i).ne.0) then
                icou= icou + 1
                  in= NODGRPSTACK(ig-1) + icou
                NODGRPITEM(in)= IMASK1(i)
              endif
            enddo
            NODGRPSTACK(ig)= NODGRPSTACK(ig-1) + icou
          enddo

          SUFGRPSTACK(0)= 0
          do ig= 1, SUFGRPTOT
            icou= 0
            do is= SUFGRPSTACKG(ig-1)+1,SUFGRPSTACKG(ig)
              i= SUFGRPITEMG(1,is)
              if (IMASK2(i).ne.0) then
                icou= icou + 1
                  in= SUFGRPSTACK(ig-1) + icou
                SUFGRPITEM(in,1)= IMASK2(i)
                SUFGRPITEM(in,2)= SUFGRPITEMG(2,is)
              endif
             enddo
            SUFGRPSTACK(ig)= SUFGRPSTACK(ig-1) + icou
          enddo
!C===

!C
!C +-------------------------+
!C | write FINAL LOCAL files |
!C +-------------------------+
!C===
        allocate (OLDtoNEW(INODTOT), NEWtoOLD(INODTOT))
        do i= 1, INODTOT
          OLDtoNEW(i)= i
          NEWtoOLD(i)= i
        enddo

        if (iopt.eq.1) call REORDERING

        open (12,file=FILNAME(ip), status='unknown')

        rewind (12)
          if (iopt.eq.0) then
            write (* ,'("PE:", 3i10)') ip, INODTOT, N2
           else
            write (* ,'("PE:", 4i10)') ip, INODTOT, N2, NCOLORtot
          endif

          write(12,'( a )'  ) '!'
          write(12,'( a )'  )                                           &
     &      '! 1.parallel information'
          write(12,'(10i10)')  ip-1
          write(12,'(10i10)')  N3
          write(12,'(10i10)') (IWKM(inei)-1,inei=1,N3)
          write(12,'( a )'  )  '!'
          write(12,'( a )'  )                                           &
     &      '! 2.mesh information (nodes and elements in partition)'
          write(12,'( a )  ')                                           &
     &      '! 2.1 node'
          write(12,'(10i10)') INODTOT, N2

          do i= 1, INODTOT
            in = NEWtoOLD(i)
            inn= INODLOCAL(in)
            write (12,'(i10,3( 1pe16.6 ))') G_NODE(INODLOCAL(in)),      & 
     &                                  (XYZ(inn,k),k=1,3)
          enddo

          write(12,'( a )'  )  '! 2.2 element (connection)'
          write(12,'(10i10)') NPC(ip)

          write(12,'(10i10)') (IELMTYP(IELMLOCAL(i)), i=1,NPC(ip))
          do i= 1, IELMTOT
            write (12,'(30i10)') G_ELEM(IELMLOCAL(i)),                  &
     &                          (ICELNOD(i,k),k=1,NODELM(IELMLOCAL(i)))
          enddo

          write(12,'( a )'  )  '!'
          write(12,'( a )'  )  '! 3.import / export information'
          write(12,'( a )'  )  '! 3.1 import'
          STACK_IMPORT(0,ip)= 0
            write (12,'(10i10)') (STACK_IMPORT(k,ip), k=1, N3)
            write (12,'(10i10)') (NOD_IMPORT(is),                       &
     &                            is= 1, STACK_IMPORT(N3,ip))

          write(12,'( a )'  )  '!'
          write(12,'( a )'  )  '! 3.2 export'
          write(12,'(10i10)') (STACK_EXPORT(inei,ip),                   &
     &                                     inei= 1, N3)
          write(12,'(10i10)') (OLDtoNEW(NOD_EXPORT(is)),                &
     &                         is= 1, STACK_EXPORT(N3,ip))
          deallocate (NOD_EXPORT)

          write(12,'( a )') '!'
          write(12,'( a )') '! boundary condition'
          write(12,'( a )') '! 4. group information'
          write(12,'( a )') '! 4.1 node group'

          call DATA_COMPRESS (NODGRPTOT, IGTOT, NODGRPSTACK, IWORK)

          if (IWORK(IGTOT).eq.0) IGTOT= 0
          write (12, '(  i10)')                  NODGRPTOT
          write (12, '(10i10)') (IWORK(ig), ig=1,NODGRPTOT)

          do ig= 1, NODGRPTOT
            write (12, '(a64)')    NODGRPNAME(ig)
            write (12, '(10i10)') (OLDtoNEW(NODGRPITEM(is)),            &
     &                            is= NODGRPSTACK(ig-1)+1,              &
     &                                NODGRPSTACK(ig) )
          enddo

          write(12,'( a )') '!'
          write(12,'( a )') '! 4.2 element group'

          call DATA_COMPRESS (ELMGRPTOT, IGTOT, ELMGRPSTACK, IWORK)

          if (IWORK(IGTOT).eq.0) IGTOT= 0
          write (12, '(  i10)')                  ELMGRPTOT
          write (12, '(10i10)') (IWORK(ig), ig=1,ELMGRPTOT)

          do ig= 1, ELMGRPTOT
            write (12, '(a64)')    ELMGRPNAME(ig)
            write (12, '(10i10)') (ELMGRPITEM(is),                      &
     &                            is= ELMGRPSTACK(ig-1)+1,              &
     &                                ELMGRPSTACK(ig) )
          enddo

          write(12,'( a )') '!'
          write(12,'( a )') '! 4.3 element-surf. group'

          call DATA_COMPRESS (SUFGRPTOT, IGTOT, SUFGRPSTACK, IWORK)

          if (IWORK(IGTOT).eq.0) IGTOT= 0
          write (12, '(  i10)')                  SUFGRPTOT
          write (12, '(10i10)') (IWORK(ig), ig=1,SUFGRPTOT)

          do ig= 1, SUFGRPTOT
              write (12, '(a64)')   SUFGRPNAME (ig)
              write (12, '(10i10)')(SUFGRPITEM(is,1),                   &
     &                              is= SUFGRPSTACK(ig-1)+1,            &
     &                                  SUFGRPSTACK(ig) )
              write (12, '(10i10)')(SUFGRPITEM(is,2),                   &
     &                              is= SUFGRPSTACK(ig-1)+1,            &
     &                                  SUFGRPSTACK(ig) )
          enddo

        close (12)
        deallocate (OLDtoNEW, NEWtoOLD)
!C===
      enddo

      return

      contains
     
        subroutine REORDERING

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

!C
!C +------------------+
!C | connected region |
!C +------------------+
!C===
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
    
      end subroutine LOCAL_DATA








