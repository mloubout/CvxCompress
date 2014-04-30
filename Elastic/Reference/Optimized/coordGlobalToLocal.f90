subroutine coordGlobalToLocal(O,U,V,W,orderaxes,nx,ny,nz,dx,dy,dz,&
  & xsou,ysou,zsou)

implicit none

double precision, intent(in) :: O(3),U(3),V(3),W(3)
character (len=3), intent(in) :: orderaxes  ! fast-med-slow axes ordering: xyz, zxy, xyz
double precision, intent(inout) :: xsou,ysou,zsou
integer, intent(in) :: nx,ny,nz
real, intent(in) :: dx,dy,dz

double precision :: f(3), g(3), h(3)

f(1) = xsou
f(2) = ysou
f(3) = zsou

g = f - O

h(1) = dot_product(g,U) / dot_product(U,U)
h(2) = dot_product(g,V) / dot_product(V,V)
h(3) = dot_product(g,W) / dot_product(W,W)

if (orderaxes == 'xyz') then
   xsou = h(1)
   ysou = h(2)
   zsou = h(3) 
elseif (orderaxes == 'zxy') then
   xsou = h(2)
   ysou = h(3)
   zsou = h(1)
elseif (orderaxes == 'zyx') then
   xsou = h(3)
   ysou = h(2)
   zsou = h(1)
endif
xsou = xsou * (nx-1) * dx
ysou = ysou * (ny-1) * dy
zsou = zsou * (nz-1) * dz

write(*, '(F12.2,A,F12.2,A,F12.2,A,F12.2,A,F12.2,A,F12.2)'), f(1),', ',f(2),', ',f(3),' ==> ',xsou,', ',ysou,', ',zsou

return
end subroutine
