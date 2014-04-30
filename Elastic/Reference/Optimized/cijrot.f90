 subroutine cijrot(theta,c11,c22,c33,c44,c55,c66,c12,c13,c14,c15,c16,c23,c24,c25,c26, &
                 & c34,c35,c36,c45,c46,c56)

! Program computes the elastic moduli after a counterclockwise
! rotation of theta about the 2-axis.
!
!
!  o Last Modified: 
!    09-10-09  Changed sign on sine to be consistent with a positive angle 
!              --> clockwise rotation about x2.
!    03-31-09  Changed argument list so that the full diagonal + upper 
!              triangle of cIJ's is passed to this routine.
!
!  o Problem Geometry (note: 3-axis is pointing down):
!
!                3         
!                ^        3'
!                |      + 
!                |    +   
!                |  +  
!                |+   
!              2 o----------->1
!                  + ) theta = pos.angle for clockwise rot. from 1 -> 1'
!                    +
!                      +
!                        1'
!                            
!
!
! o Program Notes:
!
!  Performs Bond transformation for a clockwise rotation about the 2-axis:
!
!         [c11  c12  c13  c14  c15  c16 ]
!         [c12  c22  c23  c24  c25  c26 ]
!         [c13  c23  c33  c34  c35  c36 ]
!  cIJ' = [c14  c24  c34  c44  c45  c46 ]
!         [c15  c25  c35  c45  c55  c56 ]
!         [c16  c26  c36  c46  c56  c66 ]
!
!                      T
!  cIJ = [M] [cIJ'] [M]
!
!  where [M] is the Bond transformation matrix for the 6x6 cIJ for a
!  rotation about the y-axis (see Mavko et al., The Rock Physics Handbook,
!  pp.14, 1998).
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
! Dimension Arrays and Variables:
parameter(tol=-1.e30)
real, dimension(6,6) :: cm, c, cdum, m, mt
real, dimension(3,3) :: b

! Zero stiffness matrix:
c = 0.

! 21 independent cIJ's:
c(1,1) = c11
c(2,2) = c22
c(3,3) = c33
c(4,4) = c44
c(5,5) = c55
c(6,6) = c66

c(1,2) = c12
c(1,3) = c13
c(1,4) = c14
c(1,5) = c15
c(1,6) = c16

c(2,3) = c23
c(2,4) = c24
c(2,5) = c25
c(2,6) = c26

c(3,4) = c34
c(3,5) = c35
c(3,6) = c36

c(4,5) = c45
c(4,6) = c46

c(5,6) = c56

! symmetric (sub-diagonal) cIJ's:
c(2,1) = c12
c(3,1) = c13
c(4,1) = c14
c(5,1) = c15
c(6,1) = c16

c(3,2) = c23
c(4,2) = c24
c(5,2) = c25
c(6,2) = c26

c(4,3) = c34
c(5,3) = c35
c(6,3) = c36

c(5,4) = c45
c(6,4) = c46

c(6,5) = c56

! .Compute 3x3 Rotation matrix [b] for Counter-Clockwise Rotation About 2-Axis:
b(1,1) =  cos(theta)
b(1,2) =  0.
b(1,3) = -sin(theta)  ! minus sign is consistent with a clockwise rotation about x2 being a positive angle
b(2,1) =  0.
b(2,2) =  1.
b(2,3) =  0.
b(3,1) =  sin(theta)  ! plus sign is consistent with a clockwise rotation about x2 being a positive angle
b(3,2) =  0.
b(3,3) =  cos(theta)

! .Compute Bond Transformation Matrix [M] using [b]:
m(1,1)=b(1,1)**2
m(1,2)=b(1,2)**2
m(1,3)=b(1,3)**2
m(1,4)=2.*b(1,2)*b(1,3)
m(1,5)=2.*b(1,3)*b(1,1)
m(1,6)=2.*b(1,1)*b(1,2)

m(2,1)=b(2,1)**2
m(2,2)=b(2,2)**2
m(2,3)=b(2,3)**2
m(2,4)=2.*b(2,2)*b(2,3)
m(2,5)=2.*b(2,3)*b(2,1)
m(2,6)=2.*b(2,1)*b(2,2)

m(3,1)=b(3,1)**2
m(3,2)=b(3,2)**2
m(3,3)=b(3,3)**2
m(3,4)=2.*b(3,2)*b(3,3)
m(3,5)=2.*b(3,3)*b(3,1)
m(3,6)=2.*b(3,1)*b(3,2)

m(4,1)=b(2,1)*b(3,1)
m(4,2)=b(2,2)*b(3,2)
m(4,3)=b(2,3)*b(3,3)
m(4,4)=b(2,2)*b(3,3)+b(2,3)*b(3,2)
m(4,5)=b(2,1)*b(3,3)+b(2,3)*b(3,1)
m(4,6)=b(2,2)*b(3,1)+b(2,1)*b(3,2)
  
m(5,1)=b(3,1)*b(1,1)
m(5,2)=b(3,2)*b(1,2)
m(5,3)=b(3,3)*b(1,3)
m(5,4)=b(1,2)*b(3,3)+b(1,3)*b(3,2)
m(5,5)=b(1,1)*b(3,3)+b(1,3)*b(3,1)
m(5,6)=b(1,1)*b(3,2)+b(1,2)*b(3,1)

m(6,1)=b(1,1)*b(2,1)
m(6,2)=b(1,2)*b(2,2)
m(6,3)=b(1,3)*b(2,3)
m(6,4)=b(2,2)*b(1,3)+b(1,2)*b(2,3)
m(6,5)=b(1,1)*b(2,3)+b(1,3)*b(2,1)
m(6,6)=b(2,2)*b(1,1)+b(1,2)*b(2,1)

! .Perform Matrix Multiplication:  cIJ'=[M][cIJ][M]T

! ..Zero arrays:
cdum = 0
cm = 0

! ..First compute [M][cIJ] multiplication:
do i=1,6

  do j=1,6
    cdum(i,j)=0.

    do k=1,6
      cdum(i,j)=cdum(i,j)+m(i,k)*c(k,j)
    enddo

  enddo

enddo
   
! ..Compute [M]T:
do i=1,6

  do j=1,6
    mt(i,j)=m(j,i)
  enddo

enddo

! ..Compute cIJ'=([M][cIJ])[M]T:
do i=1,6

  do j=1,6
    cm(i,j)=0.

    do k=1,6
      cm(i,j)=cm(i,j)+cdum(i,k)*mt(k,j)
    enddo

  enddo

enddo

! Re-assign for return to main routine:
c11=cm(1,1)
c22=cm(2,2)
c33=cm(3,3)
c44=cm(4,4)
c55=cm(5,5)
c66=cm(6,6)
c12=cm(1,2)
c13=cm(1,3)
c14=cm(1,4)
c15=cm(1,5)
c16=cm(1,6)
c23=cm(2,3)
c24=cm(2,4)
c25=cm(2,5)
c26=cm(2,6)
c34=cm(3,4)
c35=cm(3,5)
c36=cm(3,6)
c45=cm(4,5)
c46=cm(4,6)
c56=cm(5,6)

return
end subroutine
