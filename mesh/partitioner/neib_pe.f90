!-----------------------------------------------------------------------
!   Copyright 1999, 2000, 2001, 2002, 2003 by the Research Organization
!   for Information Science & Technology (RIST)
!-----------------------------------------------------------------------

      subroutine NEIB_PE
      use analyzer

      do ip= 1, NP
        NEIBPETOT(ip)= 0
      do is= ISTACKC(ip-1)+1, ISTACKC(ip)
        icel= NPCID(is)
      do  k= 1, NODELM(icel)
        in= ICELNOD(icel,k)
        ig= IGROUP (in)

        if (ig.ne.ip) call FIND_NEIBPE (ig, ip)
      enddo
      enddo
      enddo
      
      write ( *,'(/," PE/NEIB-PE#    NEIB-PEs")')
      write (21,'(/," PE/NEIB-PE#    NEIB-PEs")')
      do ip= 1, NP
        write ( *,'(i3,i4,5x, 31i4)') ip, NEIBPETOT(ip),                &
     &           (NEIBPE(ip,k),k=1,NEIBPETOT(ip))
        write (21,'(i3,i4,5x, 31i4)') ip, NEIBPETOT(ip),                &
     &           (NEIBPE(ip,k),k=1,NEIBPETOT(ip))
      enddo

      return
      end




      subroutine FIND_NEIBPE (ig, ip)
      use analyzer

      do inei= 1, NEIBPETOT(ip)
        if (ig.eq.NEIBPE(ip,inei)) return
      enddo

                NEIBPETOT(ip) = NEIBPETOT(ip) + 1
      NEIBPE(ip,NEIBPETOT(ip))= ig

      return
      end
