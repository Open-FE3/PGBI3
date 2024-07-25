!-----------------------------------------------------------------------
!   Copyright 1999, 2000, 2001, 2002, 2003 by the Research Organization
!   for Information Science & Technology (RIST)
!-----------------------------------------------------------------------

      subroutine PARASET
      use analyzer
      character*80 HEADER, HEADW

      N2= max (IELMTOT, N)
      allocate (IGROUP(N2))

      write (*,'(/," *********************************")')
      write (*,'(  " ||                             ||")')
      write (*,'(  " || GeoFEM Partitioner Ver.5.00 ||")')
      write (*,'(  " ||                             ||")')
      write (*,'(  " *********************************")')
      write (*,*) ' '

!C
!C-- RCB/GREEDY

  100   continue
        write (*,'(/,"# select PARTITIONING METHOD")')
        write (*,'(  "  RCB                      (1)")')
        write (*,'(  "  no longer available      (2)")')
        write (*,'(  "  MeTiS/RSB                (3)")')
        write (*,'(  "  generate MeTiS/RSB INPUT (4)")')
        write (*,*) ' '
        write (*,*) 'Please TYPE 1,2,3 or 4 !!'
        write (*,'(/,">>>")')
        read  (*,*)   NTYP
        write (*,*)  NTYP
        if (NTYP.lt.1 .or. NTYP.gt.4 .or. NTYP.eq.2) goto 100

      if (NTYP.eq.1) then
  110 continue
          write (*,'(/,"*** RECURSIVE COORDINATE BISECTION (RCB)")')
          write (*,*) 'How many partitions (2**n)?'
          write (*,'(/,">>>")')
          read  (*,*)  NPOWER
          write (*,*)  NPOWER
          if (NPOWER.lt.0) goto 110
        NP= 2**NPOWER
      endif

      if (NTYP.eq.2) then
          write (*,'(/,"*** GREEDY")')
          write (*,*) 'How many partitions ?'
          write (*,'(/,">>>")')
          read  (*,*)  NP
      endif

      if (NTYP.eq.3) then
        write (*,'(/,"*** MeTiS/RSB")')
        write (*,*      ) 'MeTiS/RSB partition-table FILE ?'
        write (*,'(/,">>>")')
        read  (*,'(a80)')  METISFIL

        open (85,file=METISFIL,status='unknown')
          NPARTMAX=-100
          do i= 1, N
            read  (85,*,err=998, end=999) ii
            IGROUP(i)= ii + 1
            NPARTMAX = max(NPARTMAX,ii+1)
          enddo

          NP= NPARTMAX
          if (NP.lt.1) call ERROR_EXIT(32,0)

          write (*,*) 'How many partitions ?'
          write (*,'(/,">>>")')
          write (*,*)  NP
      endif

      if (NTYP.eq.4) then
        call CREATE_METIS_INPUT
      endif

      write (*,'(//,"*** ",i3," REGIONS",//)') NP

      if (NP.gt.N) call ERROR_EXIT(32,N)
!C
!C-- allocation
      allocate (STACK_EXPORT(0:NP,NP))
      allocate (STACK_IMPORT(0:NP,NP))
      allocate (NPN(NP))
      allocate (NPC(NP))
      allocate (ISTACKN(0:NP))
      allocate (ISTACKC(0:NP))
      allocate (NEIBPE   (NP,NP))
      allocate (NEIBPETOT(NP))
      allocate (   NODTOT(NP))
      allocate (INTNODTOT(NP))

      allocate (IWORK(0:N2))
      allocate (IMASK(-N2:+N2))
      allocate (IDEAD (N2))
      allocate (ISTACK(N2))

      allocate (ICOND1(N2))
      allocate (ICOND2(N2))

      allocate (INODLOCAL(G_NODE_MAX))

      allocate (NOD_IMPORT(  N2))
!C
!C-- FILE NAME
      do 
        write (*,'(a)') '# HEADER of the OUTPUT file ?'
        write (*,'(a)') '  HEADER should not be <work>'
        write (*,'(/,">>>")')
        read  (*,'(a80)') HEADER
        if (HEADER.ne.'work') exit
      enddo

      HEADW= 'work'
      allocate (FILNAME(NP), WORKFIL(NP))

      do my_rank= 0, NP-1
        call DEFINE_FILE_NAME (HEADER, my_rank, FILNAME(my_rank+1))
        write (*,'(i5,5x,a80)') my_rank, FILNAME(my_rank+1)
      enddo

      do my_rank= 0, NP-1
        call DEFINE_FILE_NAME (HEADW, my_rank, WORKFIL(my_rank+1))
      enddo

!C
!C-- PARAMETER SET
 
      inum   = N / NP
      idev   = N - inum * NP

      do ip= 1, NP
        NPN(ip)= inum
      enddo

      do ip= 1, idev
        NPN(ip)= NPN(ip) + 1
      enddo

      do i= 1, N
        if (NTYP.ne.3) IGROUP(i)= 0
        RHO   (i)= 0
        IMASK (i)= 0
        IDEAD (i)= 0
      enddo

      return

 998  continue
        call ERROR_EXIT (21,0)
 999  continue
        call ERROR_EXIT (22,0)

      end

!C
!C***
!C*** DEFINE_FILE_NAME
!C***
!C
      subroutine DEFINE_FILE_NAME (HEADER, my_rank, filname)

      character (len=80) ::  filname
      character (len=80) ::  HEADER
      character (len= 1) ::  SUBindex1
      character (len= 2) ::  SUBindex2
      character (len= 3) ::  SUBindex3
      character (len= 4) ::  SUBindex4
      character (len= 5) ::  SUBindex5
      character (len= 6) ::  SUBindex6

      HEADER= adjustL (HEADER)
      LENGTH= len_trim(HEADER)

      if (my_rank.le.9) then
        ID= 1
        write(SUBindex1 ,'(i1.1)') my_rank
       else if (my_rank.le.99) then
        ID= 2
        write(SUBindex2 ,'(i2.2)') my_rank
       else if (my_rank.le.999) then
        ID= 3
        write(SUBindex3 ,'(i3.3)') my_rank
       else if (my_rank.le.9999) then
        ID= 4
        write(SUBindex4 ,'(i4.4)') my_rank
       else if (my_rank.le.99999) then
        ID= 5
        write(SUBindex5 ,'(i5.5)') my_rank
       else if (my_rank.le.999999) then
        ID= 6
        write(SUBindex6 ,'(i6.6)') my_rank
      endif

      if (ID.eq.1) filname= HEADER(1:LENGTH)//'.'//SUBindex1
      if (ID.eq.2) filname= HEADER(1:LENGTH)//'.'//SUBindex2
      if (ID.eq.3) filname= HEADER(1:LENGTH)//'.'//SUBindex3
      if (ID.eq.4) filname= HEADER(1:LENGTH)//'.'//SUBindex4
      if (ID.eq.5) filname= HEADER(1:LENGTH)//'.'//SUBindex5
      if (ID.eq.6) filname= HEADER(1:LENGTH)//'.'//SUBindex6

      return
      end



