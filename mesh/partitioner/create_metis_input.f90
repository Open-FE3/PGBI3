!-----------------------------------------------------------------------
!   Copyright 1999, 2000, 2001, 2002, 2003 by the Research Organization
!   for Information Science & Technology (RIST)
!-----------------------------------------------------------------------

      subroutine CREATE_METIS_INPUT
      use analyzer
!C
!C-- init.
      allocate (NEIBNODTOT(N))
      do i= 1, N
        NEIBNODTOT(i)= 0
      enddo
      do ie= 1, IEDGTOT
        in1= IEDGNOD   (ie,1)
        in2= IEDGNOD   (ie,2)
        ik1= NEIBNODTOT(in1) + 1
        ik2= NEIBNODTOT(in2) + 1
        NEIBNODTOT(in1)    = ik1
        NEIBNODTOT(in2)    = ik2
      enddo

      NEIBMAXmetis= -100
      do i= 1, N
        NEIBMAXmetis= max(NEIBMAXmetis,NEIBNODTOT(i))
      enddo

      allocate (NEIBNOD   (N,NEIBMAXmetis))
      do i= 1, N
        NEIBNODTOT(i)= 0
      do k= 1, NEIBMAXmetis
        NEIBNOD   (i,k)= 0
      enddo
      enddo
!C
!C-- neighboring NODEs
      do ie= 1, IEDGTOT
        in1= IEDGNOD   (ie,1)
        in2= IEDGNOD   (ie,2)
        ik1= NEIBNODTOT(in1) + 1
        ik2= NEIBNODTOT(in2) + 1
        NEIBNOD   (in1,ik1)= in2
        NEIBNOD   (in2,ik2)= in1
        NEIBNODTOT(in1)    = ik1
        NEIBNODTOT(in2)    = ik2
      enddo
!C
!C-- output
      open (21,file='partition.log',status='unknown')
        write ( *,'(/,"generate MeTiS/RSB input")')
        write (21,'(/,"generate MeTiS/RSB input")')
      
        write  ( *,'(/,"*** GRID  file   ", a80)')  GRIDFIL
        write  (21,'(/,"*** GRID  file   ", a80)')  GRIDFIL

        write  ( *,'(/,"# MeTiS  input file ?  ")')
        write  ( *,'(/,">>>")')
        read   ( *,'(a80)') METISFIL

        write  ( *,'(/,"*** MeTiS/RSB input  ", a80)')  METISFIL
        write  (21,'(/,"*** MeTiS/RSB input  ", a80)')  METISFIL
      close (21)

      open  (22,file=METISFIL,status='unknown')
        write (22,'(2i12)') N, IEDGTOT
        do in= 1, N
          write (22,'(5120i12)') (NEIBNOD(in,k),k=1,NEIBNODTOT(in)) 
        enddo
      close (22)

      call  MPI_FINALIZE (ierr)
      stop
      end
