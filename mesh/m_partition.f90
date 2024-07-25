! Coded by T. MINAMI on 2015.06.02
! to share the allocatable variable between msh2GeoFEM.f90 and the subroutine CRE_LOCAL_DATA in it.
module partition
implicit real(selected_real_kind(8))(a-h,o-z)
integer(4) :: N2n, N2c ! dimension for NPNID and NPCID
integer(4),allocatable,dimension(:) :: NPNID, NPCID
integer(4),allocatable,dimension(:) :: NOD_EXPORT
end module
