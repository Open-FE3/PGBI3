!-----------------------------------------------------------------------
!   Copyright 1999, 2000, 2001, 2002, 2003 by the Research Organization
!   for Information Science & Technology (RIST)
!-----------------------------------------------------------------------

      subroutine CRE_LOCAL_DATA
      use analyzer
      implicit REAL*8 (A-H,O-Z)

      N2n= 2 * max (IELMTOT, N)
      N2c= 2 * max (IELMTOT, N)
      allocate (NPNID(N2n))
      allocate (NPCID(N2c))

      ISTACKN(0)= 0
      ISTACKC(0)= 0

  100 continue
      do ip= 1, NP
        do icel= 1, IELMTOT
          ISTACK(icel)= 0
        enddo

        do icel= 1, IELMTOT
        do    k= 1, NODELM(icel)
          in= ICELNOD(icel,k)
          ig= IGROUP (in)

          if (ig.eq.ip) ISTACK(icel)= 1
        enddo
        enddo

        icou= 0
        do icel= 1, IELMTOT
          if (ISTACK(icel).eq.1) then
            icou= icou + 1
              is= ISTACKC(ip-1) + icou
              if (is.gt.N2c) then
                deallocate (NPCID)
                N2c= N2c * 11/10 + 1
                  allocate (NPCID(N2c))
                goto 100
              endif
            ISTACKC(ip)= is
              NPC  (ip)= icou
              NPCID(is)= icel
          endif
        enddo
      enddo

      do ip= 1, NP
        NPN(ip)= 0
      enddo

      do i= 1, N
        ig= IGROUP(i)
        NPN(ig)= NPN(ig) + 1
      enddo

      do ip= 1, NP
        ISTACKN(ip)= ISTACKN(ip-1) + NPN(ip)
        NPN(ip)= 0
      enddo

      do i= 1, N
        ip = IGROUP(i)
        icou= NPN(ip) + 1
        is = ISTACKN(ip-1) + icou
              if (is.gt.N2n) then
                deallocate (NPNID)
                N2n= N2n * 11/10 + 1
                  allocate (NPNID(N2n))
                goto 100
              endif
          NPN(ip)= icou
        NPNID(is)= i
      enddo

      MAXN= NPN(1)
      MINN= NPN(1)
      MAXC= NPC(1)
      MINC= NPC(1)
      
      do ip= 2, NP
        MAXN= max (MAXN,NPN(ip))
        MINN= min (MINN,NPN(ip))
        MAXC= max (MAXC,NPC(ip))
        MINC= min (MINC,NPC(ip))
      enddo

      write ( *,'(/,"TOTAL NODE     #   ", i12)') N
      write ( *,'(  "TOTAL CELL     #   ", i12)') IELMTOT
      write ( *,'(/," PE    NODE#   CELL#")')
      write (21,'(/,"TOTAL NODE     #   ", i12)') N
      write (21,'(  "TOTAL CELL     #   ", i12)') IELMTOT
      write (21,'(/," PE    NODE#   CELL#")')

      do ip= 1, NP
        write ( *,'(i3,5i8)') ip, NPN(ip), NPC(ip)
        write (21,'(i3,5i8)') ip, NPN(ip), NPC(ip)
      enddo

      write ( *,'(/,"MAX.node/PE        ", i12)') MAXN
      write ( *,'(  "MIN.node/PE        ", i12)') MINN
      write ( *,'(  "MAX.cell/PE        ", i12)') MAXC
      write ( *,'(  "MIN.cell/PE        ", i12)') MINC
      write (21,'(/,"MAX.node/PE        ", i12)') MAXN
      write (21,'(  "MIN.node/PE        ", i12)') MINN
      write (21,'(  "MAX.cell/PE        ", i12)') MAXC
      write (21,'(  "MIN.cell/PE        ", i12)') MINC

      return
      end







