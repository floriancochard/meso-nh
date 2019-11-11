!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
      SUBROUTINE gauher(x,w,n)
      INTEGER n,MAXIT
      REAL w(n),x(n)
      REAL EPS,PIM4
      PARAMETER (EPS=3.D-14,PIM4=.7511255444649425D0,MAXIT=10)
      INTEGER i,its,j,m
      REAL p1,p2,p3,pp,z,z1
C
      REAL SUMW
C
      m=(n+1)/2
      do 13 i=1,m
        if(i.eq.1)then
          z=sqrt(float(2*n+1))-1.85575*(2*n+1)**(-.16667)
        else if(i.eq.2)then
          z=z-1.14*n**.426/z
        else if (i.eq.3)then
          z=1.86*z-.86*x(1)
        else if (i.eq.4)then
          z=1.91*z-.91*x(2)
        else
          z=2.*z-x(i-2)
        endif
        do 12 its=1,MAXIT
          p1=PIM4
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=z*sqrt(2.d0/j)*p2-sqrt(dble(j-1)/dble(j))*p3
11        continue
          pp=sqrt(2.d0*n)*p2
          z1=z
          z=z1-p1/pp
          if(abs(z-z1).le.EPS)goto 1
12      continue
1       x(i)=z
        x(n+1-i)=-z
        pp=pp/PIM4 ! NORMALIZATION
        w(i)=2.d0/(pp*pp)
        w(n+1-i)=w(i)
13    continue
C
C NORMALISATION
C
      SUMW = 0.0
      DO 14 I=1,N
      SUMW = SUMW + W(I)
14    CONTINUE
      DO 15 I=1,N
      W(I) = W(I)/SUMW
15    CONTINUE
C
      return
      END
