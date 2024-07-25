!########################################## LOCAL_DATA
! Coppied from partitioner/local_data.f90 on 2015.06.04
! This file includes subroutine LOCAL_DATA and DATA_COMPRESS
      subroutine LOCAL_DATA( &
&	                                    INODTOTG, IELMTOTG, IELMGRPTOT, NP, xyz,&
&	                                    WORKFIL, FILNAME, ELMGRPNAME, &
&                                       ELMGRPSTACKG, ELMGRPITEMG, fileheader)
	implicit real(selected_real_kind(8))(a-h,o-z)
	integer(4),intent(in)       :: INODTOTG, IELMTOTG, IELMGRPTOT, NP
	real(8),intent(in)             :: xyz(INODTOTG,3)
	character(50),intent(in) :: WORKFIL(NP), FILNAME(NP), ELMGRPNAME(IELMGRPTOT)
	integer(4),intent(in)       :: ELMGRPSTACKG(0:IELMGRPTOT), ELMGRPITEMG(IELMTOTG)
	integer(4) :: G_NODE(INODTOTG), G_ELEM(IELMTOTG)
      integer(4), dimension(:)  , allocatable :: IWKM, OLDtoNEW, NEWtoOLD
      integer(4), dimension(:)  , allocatable :: IMASK1, IMASK2, IELMLOCAL, ELMGRPITEM
	integer(4), dimension(:,:), allocatable :: ICELNOD
	integer(4), dimension(:)  , allocatable :: INODLOCAL, NOD_IMPORT, NOD_EXPORT
	integer(4) :: STACK_IMPORT(0:NP, NP), STACK_EXPORT(0:NP,NP)
	integer(4) :: ELMGRPSTACK(0:IELMGRPTOT)
	integer(4) :: IWORK(0:IELMGRPTOT)
	integer(4), dimension(:,:), allocatable :: n4a !# for msh file
	real(8), dimension(:,:), allocatable :: xyza !# for msh file
	character(50) :: em3dmshfile, fileheader
	character(2) :: num
      parameter( iopt = 1 )
!      INODTOTG        ! total # of global nodes
!      IELMTOTG         ! total # of global elements
      G_NODE(1:INODTOTG)=(/(i,i=1,INODTOTG)/)
	G_ELEM(1:IELMTOTG)=(/(i,i=1,IELMTOTG)/)
!    Do you need Re-Ordering (RCM/Gibbs) ?'   iopt= 0.NO, 1.YES'
      N2= max (IELMTOTG, INODTOTG)
      allocate ( IWKM(NP), IMASK1(N2), IMASK2(N2))
	do ip= 1, NP !=================================  ip loop ========
!#[1]##  read  work."ip"
      open (11,file=WORKFIL(ip), status='unknown',form='unformatted') ! made in INTERFACE_NODES
          rewind (11)
          read (11) ID              ! PE ID
          read (11) INODTOT  ! # of nodes of internal + external nodes of ip-th PE
          read (11) N2             ! # of internal nodes of ip-th PE
          read (11) N3             ! # of neiboring PEs of ip-th PE
          read (11) (IWKM(i),i=1,N3)   ! PE id of ID's neighboring PEs
	            allocate ( INODLOCAL(INODTOT) )
          do i= 1, INODTOT
            read (11) INODLOCAL(i)  ! Global node id of internal and external nodes
          enddo
          STACK_IMPORT(0,ip)= 0
          read (11) (STACK_IMPORT(k,ip), k=1, N3) ! stacking count of import nodes from neighboring PEs
	            allocate( NOD_IMPORT(STACK_IMPORT(N3, ip)) )
          do is= 1, STACK_IMPORT(N3,ip)
            read (11) NOD_IMPORT(is)          ! Local node id of import nodes for ip-th PE
          enddo
          STACK_EXPORT(0,ip)= 0
          read (11) (STACK_EXPORT(k,ip), k=1, N3)  ! stacking count of export nodes to neighboring PEs
                  allocate ( NOD_EXPORT(STACK_EXPORT(N3,ip)) ) ! N3 is # of neighboring PEs
          do is= 1, STACK_EXPORT(N3,ip)
            read (11) NOD_EXPORT(is)          ! Local node id of export nodes from ip-th PE
          enddo
          read (11) IELMTOT ! # of elements at least one of whose nodes is related to ip-th PE
	           allocate( IELMLOCAL(IELMTOT), ICELNOD(IELMTOT,4) )
          do i= 1, IELMTOT
            read (11) IELMLOCAL(i),   (ICELNOD(i,k), k=1,4) ! Global element id
          enddo      ! IELMLOCAL(i) : global element id of i-th element of ip-th PE
        close (11) !  ICELNOD(i,k)  : local node id of k-th node of i-th element in ip-th PE
!
!#[2]##  LOCAL NUMBERING
	    IMASK1(1:INODTOTG)= 0
	    IMASK2(1:IELMTOTG)= 0
          do i= 1, INODTOT     !  # of internal + esternal nodes of ip-th PE
            IMASK1(INODLOCAL(i))= i   ! INODLOCAL: Global node id of internal + external nodes of ip-th PE
          enddo
!#[2-1]##  Generate ELMGRPSTACK and ELMGRPITEM [element group]
          do i= 1, IELMTOT    !  # of elements related to ip-th PE
            IMASK2( IELMLOCAL(i) )= i ! IELMLOCAL:  Local elemetn id of elements related to ip-th PE
          enddo
	       allocate( ELMGRPITEM(IELMTOT) ) ! ELMGRPITEM is the local element id for element group
          ELMGRPSTACK(0)= 0
          do ig= 1, IELMGRPTOT ! # of element groups
            icou= 0
            do is= ELMGRPSTACKG(ig-1)+1,ELMGRPSTACKG(ig)
              i= ELMGRPITEMG(is)     ! ELMGRPSTACKG is count for ELMGRPITEMG,
              if (IMASK2(i).ne.0) then
                icou= icou + 1
		    in= ELMGRPSTACK(ig-1) + icou ! ELMGRPSTACK for Local element in ip-th PE
                ELMGRPITEM(in)= IMASK2(i)    ! Enter local element id to ELMGRPITEM
              endif
            enddo
            ELMGRPSTACK(ig)= ELMGRPSTACK(ig-1) + icou
          enddo
!#[2-2]##  Generate NODGRPSTACK and NODGRPITEM [node group]
!          NODGRPSTACK(0)= 0
!          do ig= 1, NODGRPTOT
!            icou= 0
!            do is= NODGRPSTACKG(ig-1)+1,NODGRPSTACKG(ig)
!              i= NODGRPITEMG(is)
!              if (IMASK1(i).ne.0) then
!                icou= icou + 1
!                  in= NODGRPSTACK(ig-1) + icou
!                NODGRPITEM(in)= IMASK1(i)
!              endif
!            enddo
!            NODGRPSTACK(ig)= NODGRPSTACK(ig-1) + icou
!          enddo
!#[2-3]##  Generate SUFGRPSTACK and SUFGRPITEM [surface group]
!          SUFGRPSTACK(0)= 0 ! count of SURFGRPITEM
!          do ig= 1, SUFGRPTOT
!            icou= 0
!            do is= SUFGRPSTACKG(ig-1)+1,SUFGRPSTACKG(ig)
!              i= SUFGRPITEMG(1,is)
!              if (IMASK2(i).ne.0) then
!                icou= icou + 1
!                  in= SUFGRPSTACK(ig-1) + icou
!                SUFGRPITEM(in,1)= IMASK2(i)
!                SUFGRPITEM(in,2)= SUFGRPITEMG(2,is)
!              endif
!             enddo
!            SUFGRPSTACK(ig)= SUFGRPSTACK(ig-1) + icou
!          enddo
!
!#[3]## write FINAL LOCAL files |
        allocate (OLDtoNEW(INODTOT), NEWtoOLD(INODTOT)) ! INODTOT is # of internal and external nodes
        do i= 1, INODTOT
          OLDtoNEW(i)= i ! At first, Old and New order is the same
          NEWtoOLD(i)= i ! these may be changed by subroutine REORDERING
        enddo
!        if (iopt.eq.1) call REORDERING(INODTOTG, IELMTOT) !############# Reverse Cut Hil McKee
        open (12,file=FILNAME(ip), status='unknown')
        rewind (12)
          if (iopt.eq.0) then
            write (* ,'("PE:", 3i10)') ip, INODTOT, N2
           else
            write (* ,'("PE:", 4i10)') ip, INODTOT, N2 !, NCOLORtot
          endif
          write(12,'( a )'  ) '!'
          write(12,'( a )'  ) '! 1.parallel information'
          write(12,'(10i10)')  ip-1                             ! PE id should start with 0
          write(12,'(10i10)')  N3
          write(12,'(10i10)') (IWKM(inei)-1,inei=1,N3)
          write(12,'( a )'  )  '!'
          write(12,'( a )'  )   '! 2.mesh information (nodes and elements in partition)'
          write(12,'( a )  ')   '! 2.1 node'
          write(12,'(10i10)') INODTOT, N2
	    allocate(xyza(INODTOT,3), n4a(IELMTOT, 5) )   !# for msh file
          do i= 1, INODTOT
            in = NEWtoOLD(i)
            inn= INODLOCAL(in)
!            write (12,'(i10,3( 1pe16.6 ))') G_NODE(INODLOCAL(in)), (XYZ(inn,k),k=1,3)
            write (12,*) G_NODE(INODLOCAL(in)), (XYZ(inn,k),k=1,3)
		xyza(i,1:3)=xyz(inn,1:3)     !# for msh file
          enddo
          write(12,'( a )'  )  '! 2.2 element (connection)'
!          write(12,'(10i10)') NPC(ip)
!          write(12,'(10i10)') (IELMTYP(IELMLOCAL(i)), i=1,NPC(ip))
          IELMTYP=311 ! 1st order tetrahedron
          write(12,'(10i10)') IELMTOT
          write(12,'(10i10)') (IELMTYP, i=1,IELMTOT)
          do i= 1, IELMTOT
            write (12,'(30i10)') G_ELEM(IELMLOCAL(i)),  (ICELNOD(i,k),k=1,4)
		n4a(i,2:5)=ICELNOD(i,1:4) !# for msh file
          enddo
          write(12,'( a )'  )  '!'
          write(12,'( a )'  )  '! 3.import / export information'
          write(12,'( a )'  )  '! 3.1 import'
          STACK_IMPORT(0,ip)= 0
            write (12,'(10i10)') (STACK_IMPORT(k,ip), k=1, N3)
            write (12,'(10i10)') (NOD_IMPORT(is),  is= 1, STACK_IMPORT(N3,ip))
          write(12,'( a )'  )  '!'
          write(12,'( a )'  )  '! 3.2 export'
          write(12,'(10i10)') (STACK_EXPORT(inei,ip),  inei= 1, N3)
          write(12,'(10i10)') (OLDtoNEW(NOD_EXPORT(is)), is= 1, STACK_EXPORT(N3,ip))
          write(12,'( a )') '!'
          write(12,'( a )') '! boundary condition'
          write(12,'( a )') '! 4. group information'
          write(12,'( a )') '! 4.1 node group'
!          call DATA_COMPRESS (NODGRPTOT, IGTOT, NODGRPSTACK, IWORK)
!          if (IWORK(IGTOT).eq.0) IGTOT= 0
!          write (12, '(  i10)')                  NODGRPTOT
!          write (12, '(10i10)') (IWORK(ig), ig=1,NODGRPTOT)
!          do ig= 1, NODGRPTOT
!	     write (12, '(a64)')    NODGRPNAME(ig)
!	     write (12, '(10i10)') (OLDtoNEW(NODGRPITEM(is)), is= NODGRPSTACK(ig-1)+1, NODGRPSTACK(ig) )
!          enddo
          write(12,'( a )') '!'
          write(12,'( a )') '! 4.2 element group'
          call DATA_COMPRESS (IELMGRPTOT, IGTOT, ELMGRPSTACK, IWORK)
          if (IWORK(IGTOT).eq.0) IGTOT= 0
          write (12, '(  i10)')                  IELMGRPTOT
          write (12, '(10i10)') (IWORK(ig), ig=1,IELMGRPTOT)
          do ig= 1, IELMGRPTOT
            write (12, '(a64)')    ELMGRPNAME(ig)
            write (12, '(10i10)') (ELMGRPITEM(is), is= ELMGRPSTACK(ig-1)+1,  ELMGRPSTACK(ig) )
		do is= ELMGRPSTACK(ig-1)+1,  ELMGRPSTACK(ig) !# for msh file
		n4a(ELMGRPITEM(is),1)=ig                                       !# for msh file
		end do                                                                         !# for msh file
          enddo
          write(12,'( a )') '!'
          write(12,'( a )') '! 4.3 element-surf. group'
!         call DATA_COMPRESS (SUFGRPTOT, IGTOT, SUFGRPSTACK, IWORK)
!         if (IWORK(IGTOT).eq.0) IGTOT= 0
!          write (12, '(  i10)')                  SUFGRPTOT
!          write (12, '(10i10)') (IWORK(ig), ig=1,SUFGRPTOT)
!          do ig= 1, SUFGRPTOT
!              write (12, '(a64)')   SUFGRPNAME (ig)
!              write (12, '(10i10)')(SUFGRPITEM(is,1), is= SUFGRPSTACK(ig-1)+1, SUFGRPSTACK(ig) )
!              write (12, '(10i10)')(SUFGRPITEM(is,2), is= SUFGRPSTACK(ig-1)+1, SUFGRPSTACK(ig) )
!          enddo
        close (12)
! Output msh file for ip-th PE
        write(num,'(i2.2)') ip
	  em3dmshfile=fileheader(1:len_trim(fileheader))//num//".msh"
	  call OUTEM3D(em3dmshfile,xyza,INODTOT,n4a,IELMTOT)
        deallocate (OLDtoNEW, NEWtoOLD, NOD_EXPORT, NOD_IMPORT)
	  deallocate ( INODLOCAL, IELMLOCAL, ICELNOD, ELMGRPITEM)
	  deallocate (xyza, n4a )
      enddo !========================================  ip loop ========
	deallocate ( IWKM, IMASK1, IMASK2)
      return
      end subroutine LOCAL_DATA

!###################################################################
!C*** DATA_COMPRESS
!C***
      subroutine DATA_COMPRESS (I1, I2, I_INN, I_OUT)
      integer  I_INN(0:I1), I_OUT(0:I1)
      ITOT= I1
      do i= 0, I1
        I_OUT(i)= I_INN(i)
      enddo
      I2= ITOT
      return
      end
!
!############################################  OUTEM3d
! Copied from combine3d.f90 on 2015.06.09
subroutine OUTEM3D(em3dfile,xyza,nodea,n4a,nteta)
implicit none
integer(4),intent(in) :: nodea,nteta
!integer(4),intent(in) :: n4s(ntets,5), n1s(npois), npois,ntets
integer(4),intent(inout) :: n4a(nteta,5)
real(8),dimension(nodea,3),intent(in) :: xyza
character(50),intent(in) :: em3dfile
integer(4) :: i, j, l
integer(4) :: ishift
open(1,file=em3dfile)
write(1,'(a11)') "$MeshFormat"
write(1,'(a7)') "2.2 0 8"
write(1,'(a14)') "$EndMeshFormat"
write(1,'(a6)') "$Nodes"
write(1,*) nodea
do i=1,nodea
write(1,*) i, (xyza(i,j),j=1,3)
end do
write(1,'(a9)') "$EndNodes"
! out oceanfile ### tetrahedron elements
write(1,'(a9)') "$Elements"
!write(1,*) npois+ntets+nteta  !ntets+nteta
write(1,*) nteta  !ntets+nteta
!do i=1,npois
!write(1,*) i,"  15   2   0   1", n1s(i) ! 1 for ocean
!end do
!ishift=npois
ishift=0
!do i=1,ntets
!write(1,*) ishift+i,"  4   2   0   1", (n4s(i,l),l=2,5) ! 1 for ocean
!end do
!n4a(:,1)=n4a(:,1)+1 ! originaly 1 for air, 2 for land
!ishift=npois+ntets
do i=1,nteta
write(1,*) i+ishift,"  4   2   0",(n4a(i,l),l=1,5) ! n4a(1)=2 for air, 3 for land
end do
write(1,'(a12)') "$EndElements"
close(1)
write(*,*) "### OUTEM3D END ###"
return
end subroutine outem3d

