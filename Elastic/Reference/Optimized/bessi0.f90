real function bessi0(x)

! Modified Bessel Function of zero order.
! From Numerical Recipes, Press et al. (1986), pp. 177

implicit none

real :: x, AX
real :: Y,P1,P2,P3,P4,P5,P6,P7 
real :: Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9

data P1,P2,P3,P4,P5,P6,P7/1.0,3.5156229,3.0899424,1.2067492, &
     & 0.2659732,0.360768e-1,0.45813e-2/
data Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228,0.1328592e-1, &
     & 0.225319e-2,-0.157565e-2,0.916281e-2,-0.2057706e-1, &
     & 0.2635537e-1,-0.1647633e-1,0.392377e-2/

if (abs(x) < 3.75) then
  Y = (x/3.75)**2
  bessi0 = P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
else
  AX = abs(x)
  Y = 3.75/AX
  bessi0 = (exp(AX)/sqrt(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4 &
          & +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
endif

return
end function
