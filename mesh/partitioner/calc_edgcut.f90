!-----------------------------------------------------------------------
!   Copyright 1999, 2000, 2001, 2002, 2003 by the Research Organization
!   for Information Science & Technology (RIST)
!-----------------------------------------------------------------------

      subroutine CALC_EDGCUT
      use analyzer
!C
!C-- calc. EDGECUT

      IEDGCUT= 0
      do ie= 1, IEDGTOT
        in1= IEDGNOD(ie,1)
        in2= IEDGNOD(ie,2)
        ig1= IGROUP(in1)
        ig2= IGROUP(in2)
        if (ig1.ne.ig2) IEDGCUT= IEDGCUT + 1
      enddo

      write ( *,'(/,"TOTAL EDGE     #   ", i12)') IEDGTOT
      write ( *,'(  "TOTAL EDGE CUT #   ", i12)') IEDGCUT
      write (21,'(/,"TOTAL EDGE     #   ", i12)') IEDGTOT
      write (21,'(  "TOTAL EDGE CUT #   ", i12)') IEDGCUT

      return
      end



