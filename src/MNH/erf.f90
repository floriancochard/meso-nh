!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
! MASDEV4_7 operators 2006/05/18 13:07:25
!-----------------------------------------------------------------
!##############
MODULE MODI_ERF
!##############
!
INTERFACE
!
FUNCTION ERF(PX)  RESULT(PERF)
REAL, INTENT(IN)                                  :: PX
REAL                                              :: PERF
END FUNCTION ERF
!
END INTERFACE
!
END MODULE MODI_ERF
!     #############################
      FUNCTION ERF(PX) RESULT(PERF)
!     #############################
!
!*       0. DECLARATIONS
!           ------------
!
IMPLICIT NONE
!
!*       0.1 declarations of arguments and result
!
REAL, INTENT(IN)                     ::   PX
REAL                                 ::   PERF
!
!*       0.2 declarations of local variables
!
REAL, DIMENSION(4), SAVE             ::   ZP1,ZQ1
REAL, DIMENSION(8), SAVE             ::   ZP2,ZQ2
REAL, DIMENSION(5), SAVE             ::   ZP3,ZQ3
REAL,               SAVE             ::   ZCONST2
REAL                                 ::   ZA,ZR,ZS,ZT,ZU,ZY
INTEGER                              ::   JENTRY,JBASIC
LOGICAL, SAVE                        ::   FIRST_TIME = .TRUE.
!
IF(FIRST_TIME)THEN
  FIRST_TIME=.FALSE.
  ZP1(1) =  2.42667955230532E+2
  ZP1(2) =  2.19792616182942E+1
  ZP1(3) =  6.99638348861914E+0
  ZP1(4) = -3.56098437018154E-2
!
  ZQ1(1) =  2.15058875869861E+2
  ZQ1(2) =  9.11649054045149E+1
  ZQ1(3) =  1.50827976304078E+1
  ZQ1(4) =  1.00000000000000E+0
!
  ZP2(1) =  3.00459261020162E+2
  ZP2(2) =  4.51918953711873E+2
  ZP2(3) =  3.39320816734344E+2
  ZP2(4) =  1.52989285046940E+2
  ZP2(5) =  4.31622272220567E+1
  ZP2(6) =  7.21175825088309E+0
  ZP2(7) =  5.64195517478974E-1
  ZP2(8) = -1.36864857382717E-7
!
  ZQ2(1) =  3.00459260956983E+2
  ZQ2(2) =  7.90950925327898E+2
  ZQ2(3) =  9.31354094850610E+2
  ZQ2(4) =  6.38980264465631E+2
  ZQ2(5) =  2.77585444743988E+2
  ZQ2(6) =  7.70001529352295E+1
  ZQ2(7) =  1.27827273196294E+1
  ZQ2(8) =  1.00000000000000E+0
!
  ZP3(1) = -2.99610707703542E-3
  ZP3(2) = -4.94730910623251E-2
  ZP3(3) = -2.26956593539687E-1
  ZP3(4) = -2.78661308609648E-1
  ZP3(5) = -2.23192459734185E-2
!
  ZQ3(1) =  1.06209230528468E-2
  ZQ3(2) =  1.91308926107830E-1
  ZQ3(3) =  1.05167510706793E+0
  ZQ3(4) =  1.98733201817135E+0
  ZQ3(5) =  1.00000000000000E+0
!
!  ZCONST1 = 0.707106781186548            ! SQRT(1/2)
  ZCONST2 = 0.564189583547756            ! SQRT(1/PI)
END IF
!
JENTRY = 1
ZT = PX
ZA = ABS(ZT)
IF(ZA.LE.6.0)GO TO 11
PERF = SIGN(1.0,PX)
RETURN
!
11 ZS = ZT*ZT
IF(ZA.GT.0.47)GO TO 20
!
JBASIC = 1
ZY = ZT*(ZP1(1)+ZS*(ZP1(2)+ZS*(ZP1(3)+ZS*ZP1(4))))/                 &
        (ZQ1(1)+ZS*(ZQ1(2)+ZS*(ZQ1(3)+ZS*ZQ1(4))))
GO TO 30
!
20 JBASIC = 2
IF(ZA.GT.4.0)GO TO 21
ZY = EXP(-ZS)*(ZP2(1)+ZA*(ZP2(2)+ZA*(ZP2(3)+ZA*(ZP2(4)+             &
               ZA*(ZP2(5)+ZA*(ZP2(6)+ZA*(ZP2(7)+ZA*ZP2(8))))))))/   &
              (ZQ2(1)+ZA*(ZQ2(2)+ZA*(ZQ2(3)+ZA*(ZQ2(4)+             &
               ZA*(ZQ2(5)+ZA*(ZQ2(6)+ZA*(ZQ2(7)+ZA*ZQ2(8))))))))
GO TO 30
!
21 ZY = 0.0
IF(ZA.GT.26.0)GO TO 30
ZR = 1.0/ZA
ZU = ZR*ZR
ZY = ZR*EXP(-ZS)*(ZCONST2+                                          &
        ZU*(ZP3(1)+ZU*(ZP3(2)+ZU*(ZP3(3)+ZU*(ZP3(4)+ZU*ZP3(5)))))/  &
           (ZQ3(1)+ZU*(ZQ3(2)+ZU*(ZQ3(3)+ZU*(ZQ3(4)+ZU*ZQ3(5))))))
!
30 CONTINUE
IF(JBASIC.EQ.2)GO TO 31
PERF = ZY
RETURN
!
31 PERF = 1.0-ZY
IF(PX.LT.0.0)PERF = -PERF
RETURN
!
END FUNCTION ERF
