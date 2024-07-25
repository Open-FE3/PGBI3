!# modified on 2018.08.28
module bgelem
implicit none
!#
type bgele
 integer(4)             :: nfield ! 2018.08.28
 real(8),allocatable,dimension(:,:) :: xbg
 real(8),allocatable,dimension(:,:) :: ybg
 real(8),allocatable,dimension(:,:) :: val
end type bgele

contains
!################################################  initbgele
!# coded on 2018.08.28
subroutine initbgele(nfield,pos)
implicit none
type(bgele),intent(out) :: pos
integer(4), intent(in)  :: nfield

 pos%nfield = nfield
 allocate(pos%xbg(nfield,2))
 allocate(pos%ybg(nfield,2))
 allocate(pos%val(nfield,4))

return
end
!################################################  setbgele
!# modified on 2018.08.28
subroutine setbgele(pos,i,y1,y2,x1,x2,v)
! in this subroutine, y is the horizontal axis, x is the vertical axis
implicit none
integer(4)                :: i,j
real(8)                   :: y1,y2,x1,x2,v(4)
type(bgele),intent(inout) :: pos

 pos%ybg(i,1:2)=(/y1,y2/)
 pos%xbg(i,1:2)=(/x1,x2/)
 pos%val(i,1:4)=v(1:4)

!write(*,*) "bgelement # =",i
!write(*,*) "pos%val(i,1-4)=",(pos%val(i,j),j=1,4)
return
end

!############################################### sizebgele
!# modified on 2018.08.28
subroutine sizebgele(x,y,pos,vout)
! x is the vertical coordinate, y is the horizontal
type(bgele) :: pos
integer(4)  :: i,j
real(8)     :: x,y,x1,x2,y1,y2,v(4),a(4),vout,area
integer(4)  :: nfield

!#[1]## set nfield 2018.08.28
 nfield = pos%nfield

!#[2]## search which element the coord is within
do i=1,nfield
 x1=pos%xbg(i,1)
 x2=pos%xbg(i,2)
 y1=pos%ybg(i,1)
 y2=pos%ybg(i,2)
 if (x1 .le. x .and. x .le. x2 .and. y1 .le. y .and. y .le. y2) then
  v(1:4)=pos%val(i,1:4)
  area=(x2-x1)*(y2-y1)
  a(1)=(x-x1)*(y2-y)/area
  a(2)=(x2-x)*(y2-y)/area
  a(3)=(x2-x)*(y-y1)/area
  a(4)=(x-x1)*(y-y1)/area
  vout=v(1)*a(1)+v(2)*a(2)+v(3)*a(3)+v(4)*a(4)
  !write(*,*) "x,y=",x,y,"bgele,ent #=",i
  !write(*,*) "a(1-4)=",(a(j),j=1,4)
  !write(*,*) "v(1-4)=",(v(j),j=1,4)
  return
 end if
end do

write(*,*) "GEGEGE no bgmesh element is selected for the coord"
write(*,*) "(x,y)=",x,y
write(*,*) "x(1)=",pos%xbg(1,2),"x(node)=",pos%xbg(9,1)
write(*,*) "y(1)=",pos%ybg(1,1),"y(node)=",pos%ybg(9,2)
stop
end subroutine sizebgele
!
subroutine outposgeo(i,pos,outfile)
implicit none
type(bgele) :: pos
integer(4) :: i
character(50) :: outfile
open(1,file=outfile)
write(1,*) "> -Z"
write(1,*)pos%ybg(i,1),pos%xbg(i,2)
write(1,*)pos%ybg(i,1),pos%xbg(i,1)
write(1,*)pos%ybg(i,2),pos%xbg(i,1)
write(1,*)pos%ybg(i,2),pos%xbg(i,2)
write(1,*)pos%ybg(i,1),pos%xbg(i,2)
close(1)
end subroutine outposgeo
end module bgelem


