 SUBROUTINE FIT(X,Y,NDATA,SIG,MWT,A,B,SIGA,SIGB,CHI2,Q)

! Linear regression routine from Numerical Recipes, pg. 508-509.
!
! Given a set of NDATA points X(I) and Y(I) with standard deviations SIG(I),
! fit them to a straight line y = a + bx by minimizing chi-squared CHI2. 
! Returned are A, B and their respective probable uncertainties SIGA and 
! SIGB, the CHI2, and the goodness-of-fit probability Q (that the fit would
! have CHI2 this large or larger). If MWT = 0 on input, then the standard
! deviations are assumed to be unavailable: Q is returned as 1.0 and the
! normalization of CHI2 is to unit standard deviation on all points.
! 

DIMENSION X(NDATA),Y(NDATA),SIG(NDATA)

SX=0.
SY=0.
ST2=0.
B=0.
IF(MWT.NE.0) THEN
  SS=0.
  DO I=1,NDATA
    WT=1./(SIG(I)**2)
    SS=SS+WT
    SX=SX+X(I)*WT
    SY=SY+Y(I)*WT
  ENDDO
ELSE
  DO I=1,NDATA
    SX=SX+X(I)
    SY=SY+Y(I)
  ENDDO
  SS=FLOAT(NDATA)
ENDIF

SXOSS=SX/SS
IF(MWT.NE.0) THEN
  DO I=1,NDATA
    T=(X(I)-SXOSS)/SIG(I)
    ST2=ST2+T*T
    B=B+T*Y(I)/SIG(I)
  ENDDO
ELSE
  DO I=1,NDATA
    T=X(I)-SXOSS
    ST2=ST2+T*T
    B=B+T*Y(I)
  ENDDO
ENDIF

B=B/ST2
A=(SY-SX*B)/SS
SIGA=SQRT((1.+SX*SX/(SS*ST2))/SS)
SIGB=SQRT(1./ST2)
CHI2=0.
IF(MWT.EQ.0) THEN
  DO I=1,NDATA
    CHI2=CHI2+(Y(I)-A-B*X(I))**2
  ENDDO
  Q=1.
  SIGDAT=SQRT(CHI2/(NDATA-2))
  SIGA=SIGA*SIGDAT
  SIGB=SIGB*SIGDAT
ELSE
  DO I=1,NDATA
    CHI2=CHI2+((Y(I)-A-B*X(I))/SIG(I))**2
  ENDDO
  !!!!Q=GAMMQ(0.5*(NDATA-2),0.5*CHI2)
ENDIF

RETURN
END
