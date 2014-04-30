subroutine coordLocalToGlobal(O,U,V,W,orderaxes,nx,ny,nz,&
  & dx,dy,dz,xsou,ysou,zsou)

implicit none

double precision, intent(in) :: O(3),U(3),V(3),W(3)
character (len=3), intent(in) :: orderaxes  ! fast-med-slow axes ordering: xyz, zxy, xyz
double precision, intent(inout) :: xsou,ysou,zsou
integer, intent(in) :: nx,ny,nz
real, intent(in) :: dx,dy,dz

double precision :: f(3)
double precision :: x,y,z

x = xsou / ((nx-1)*dx)
y = ysou / ((ny-1)*dy)
z = zsou / ((nz-1)*dz)

if (orderaxes == 'xyz') then
   f = x*U + y*V + z*W + O
elseif (orderaxes == 'zxy') then
   f = x*V + y*W + z*U + O
elseif (orderaxes == 'zyx') then
   f = x*W + y*V + z*U + O
endif

write(*, '(F12.2,A,F12.2,A,F12.2,A,F12.2,A,F12.2,A,F12.2)'), xsou,', ',ysou,', ',zsou,' ==> ',f(1),', ',f(2),', ',f(3)
xsou = f(1)
ysou = f(2)
zsou = f(3)

return
end subroutine
