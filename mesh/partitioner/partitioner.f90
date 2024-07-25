!-----------------------------------------------------------------------
!   Copyright 1999, 2000, 2001, 2002, 2003 by the Research Organization
!   for Information Science & Technology (RIST)
!-----------------------------------------------------------------------

      module partitioner
        use geofem_util
        implicit REAL*8 (A-H,O-Z)
        integer(kind=kint), dimension (:)  , allocatable :: RHO
        integer(kind=kint), dimension (:,:), allocatable :: STACK_EXPORT
        integer(kind=kint), dimension (:,:), allocatable :: STACK_IMPORT

        integer(kind=kint), dimension (:)  , allocatable ::             &
     &           ELMGRPITEM, NODGRPITEM,                                &
     &           NEIBNODTOT, IWORK, IACTEDG, IEDGFLAG, IMASK,           &
     &           IDEAD, ISTACK, IGROUP, ICOND1, ICOND2

        integer(kind=kint), dimension (:,:), allocatable ::             &
     &           SUFGRPITEM, IEDGNOD, NEIBNOD, NEIBPE

        integer(kind=kint), pointer::                                   &
     &                              NODGRPITEMG(:), ELMGRPITEMG(:),     &
     &                              SUFGRPITEMG(:,:),                   &
     &                              ELMGRPSTACKG(:), NODGRPSTACKG(:),   &
     &                              SUFGRPSTACKG(:)

        integer(kind=kint ), pointer::  ICELNOD(:,:), IELMTYP(:)
        integer(kind=kint ), pointer::  G_NODE(:), G_ELEM(:)
        real   (kind=kreal), pointer::  XYZ(:,:)

        integer(kind=kint ) :: G_NODE_MAX, G_ELEM_MAX

        integer(kind=kint ), dimension (:) , allocatable ::             &
     &           NPN, NPNID, ISTACKN, NPC, NPCID, ISTACKC, NEIBPETOT,   &
     &           NODTOT, INTNODTOT, INODLOCAL, NOD_EXPORT, NOD_IMPORT,  &
     &           ELMGRPSTACK, NODGRPSTACK , SUFGRPSTACK

        integer(kind=kint ), dimension (:) , allocatable ::  NODELM
        integer(kind=kint ), dimension (:) , allocatable ::  NUMSUF

        integer(kind=kint )                   ::                        &
     &           RHOMAX, RHOMIN, ELMGRPTOT, NODGRPTOT, SUFGRPTOT,       &
     &           N, NP, NPOWER, NTYP

        character(geofem_name_len),pointer::                            &
     &                NODGRPNAME(:), ELMGRPNAME(:), SUFGRPNAME(:)
      

        character*80                          :: GRIDFIL, METISFIL
        character (len=80), dimension(:), allocatable, save :: FILNAME
        character (len=80), dimension(:), allocatable, save :: WORKFIL

        integer(kind=kint )                   ::                        &
     &         IOPTFLAG, ITERMAX, IEDGCUT, IEDGTOT, IELMTOT, IACTEDGTOT

      end module partitioner
