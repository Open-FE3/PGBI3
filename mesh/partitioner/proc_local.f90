!-----------------------------------------------------------------------
!   Copyright 1999, 2000, 2001, 2002, 2003 by the Research Organization
!   for Information Science & Technology (RIST)
!-----------------------------------------------------------------------

      subroutine PROC_LOCAL
      use analyzer
      character*80 UCDFIL
      character* 6 ETYPE

      open (21,file='partition.log',status='unknown')

      if (NTYP.eq.1) then
        write ( *,'(/,"Recursive Coordinate Bisection")')
        write (21,'(/,"Recursive Coordinate Bisection")')
      endif

      if (NTYP.eq.2) then
        write ( *,'(/,"Greedy")')
        write (21,'(/,"Greedy")')
      endif

      if (NTYP.eq.3) then
        write ( *,'(/,"MeTiS")')
        write (21,'(/,"MeTiS")')
      endif

      write  ( *,'(/,"*** GRID  file   ", a80)')  GRIDFIL
      write  (21,'(/,"*** GRID  file   ", a80)')  GRIDFIL

      if (NTYP.eq.3) then
        write  ( *,'("*** MeTiS file   ", a80)')  METISFIL
        write  (21,'("*** MeTiS file   ", a80)')  METISFIL
      endif

      write ( *,'(/,i5, " PEs")') NP
      write (21,'(/,i5, " PEs")') NP
!C
!C-- UCD file

      write (*,'(//,a)') ' Do you want to create file in UCD format?'
      write (*,'( /,a)') ' 0.NO, 1.YES'
      write (*,'(a)')    ' Please TYPE 0 or 1 !!'
      write (*,'(/,">>>")')
      read  (*,*) iopt

      if (iopt.eq.1) then
        write (*,*      ) 'name of UCD file?'
        write (*,'(/,">>>")')
        read  (*,'(a80)')  UCDFIL

        open (22, file= UCDFIL, status='unknown')
          N1= 1
          N0= 0

          icou= 0
          do ie= 1, IELMTOT
            ityp= IELMTYP(ie)
            call UCDelemtype (ityp, ETYPE)
            if (ETYPE.ne.'NO') icou= icou + 1
          enddo

          write (22,'(5i12)') N, icou, N0, N1, N0
          do i= 1, N
            write (22,'(i12,3(1pe16.6))') i, (XYZ(i,k),k=1,3)
          enddo

          icou= 0
          do ie= 1, IELMTOT
            ityp= IELMTYP(ie)
            call UCDelemtype (ityp, ETYPE)
            if (ETYPE.ne.'NO') then
              icou= icou + 1
              write (22,'(i12,i3,1x,a6,1x,8i12)') icou, N1, ETYPE,      &
     &              (ICELNOD(ie,k), k= 1, NODELM(ie))
            endif
          enddo
          write (22,'(2i3)')  N1, N1
          write (22,'(a  )') 'COLOR, color'
          do ie= 1, IELMTOT
            in1= ICELNOD(ie,1)
            ig = IGROUP(in1)
            do k= 2, NODELM(ie)
              inK= ICELNOD(ie,k)
              igK= IGROUP (inK)
              if (igK.ne.ig) exit
            enddo
            if (ig.ne.igK) then
              ig= NP + 5
             else
              ig= ig-1
            endif
            write (22,'(2i12)') ie, ig
          enddo
        close (22)

  100   continue
        write (*,'(//,a)') ' Do you want to finish NOW ?'
        write (*,'( /,a)') ' 0.NO, 1.YES'
        write (*,'(a)')    ' Please TYPE 0 or 1 !!'
        read  (*,*) iopt2
        if (iopt2.eq.1) then
          write (*,'(//,a)') ' ARE YOU SURE ?'
          write (*,'( /,a)') ' 0.NO, 1.YES'
          write (*,'(a)')    ' Please TYPE 0 or 1 !!'
          read  (*,*) iopt3
          if (iopt3.eq.0) goto 100
          call  MPI_FINALIZE (errno)
          stop ' * normal termination'
        endif
      endif
!C
!C-- create LOCAL DATA
      call CALC_EDGCUT
      call CRE_LOCAL_DATA

!C
!C-- OVERLAPPED ELEMENTs

      do icel= 1, IELMTOT
        ISTACK(icel)= 0
      enddo

      do icel= 1, IELMTOT
        do k1= 1, NODELM(icel)
        do k2= 1, NODELM(icel)
          ig1= IGROUP(ICELNOD(icel,k1))
          ig2= IGROUP(ICELNOD(icel,k2))
          if (ig1.ne.ig2) ISTACK(icel)= 1
        enddo
        enddo
      enddo

      icou= 0
      do icel= 1, IELMTOT
        if (ISTACK(icel).eq.1) icou= icou + 1
      enddo
      write ( *,'(/,"OVERLAPPED ELEMENTS", i12)')  icou
      write (21,'(/,"OVERLAPPED ELEMENTS", i)')  icou

!C
!C-- NEIGHBORING PEs
      call NEIB_PE

!C
!C-- INTERFACE info.
      call INTERFACE_NODES
      close (21)

!C
!C-- distributed Local DATA
      call LOCAL_DATA

      return
      end

      subroutine UCDelemtype (ityp, ETYPE)
      character* 6 ETYPE

      if (ityp.eq.111 .or. ityp.eq.611) then
        ETYPE='line  '
        return
      endif

      if (ityp.eq.211 .or. ityp.eq.711) then
        ETYPE='tri   '
        return
      endif

      if (ityp.eq.221 .or. ityp.eq.721) then
        ETYPE='quad  '
        return
      endif

      if (ityp.eq.311) then
        ETYPE='tet   '
        return
      endif

      if (ityp.eq.321 .or. ityp.eq.511) then
        ETYPE='prism '
        return
      endif

      if (ityp.eq.331 .or. ityp.eq.512) then
        ETYPE='hex   '
        return
      endif

      if (ityp.eq.112) then
        ETYPE='line2 '
        return
      endif

      if (ityp.eq.212 .or. ityp.eq.712) then
        ETYPE='tri2  '
        return
      endif

      if (ityp.eq.222 .or. ityp.eq.722) then
        ETYPE='quad2 '
        return
      endif

      if (ityp.eq.312) then
        ETYPE='tet2  '
        return
      endif

      if (ityp.eq.322 .or. ityp.eq.512) then
        ETYPE='prism2'
        return
      endif

      if (ityp.eq.332 .or. ityp.eq.522) then
        ETYPE='hex2  '
        return
      endif

      ETYPE= 'NO'
      return

      end



