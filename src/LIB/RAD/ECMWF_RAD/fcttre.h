!*
!     ------------------------------------------------------------------
#ifdef DOC

!     This COMDECK includes the Thermodynamical functions for the cy39
!       ECMWF Physics package.
!       Consistent with YOMCST Basic physics constants, assuming the
!       partial pressure of water vapour is given by a first order
!       Taylor expansion of Qs(T) w.r.t. to Temperature, using constants
!       in YOETHF
!       Two sets of functions are available. In the first set only the
!       cases water or ice are distinguished by temperature.  This set 
!       consists of the functions FOEDELTA,FOEEW,FOEDE and FOELH.
!       The second set considers, besides the two cases water and ice 
!       also a mix of both for the temperature range RTICE < T < RTWAT.
!       This set contains FOEALFA,FOEEWM,FOEDEM,FOELDCPM and FOELHM.

!       Depending on the consideration of mixed phases either the first 
!       set (e.g. surface, post-processing) or the second set 
!       (e.g. clouds, condensation, convection) should be used.

!     ------------------------------------------------------------------
!     *****************************************************************

!                NO CONSIDERATION OF MIXED PHASES

!     *****************************************************************
#endif
REAL_B :: FOEDELTA
REAL_B :: PTARE
FOEDELTA (PTARE) = MAX (_ZERO_,SIGN(_ONE_,PTARE-RTT))
#ifdef DOC

!                  FOEDELTA = 1    water
!                  FOEDELTA = 0    ice

!     THERMODYNAMICAL FUNCTIONS .

!     Pressure of water vapour at saturation
!        INPUT : PTARE = TEMPERATURE
#endif
REAL_B :: FOEEW,FOEDE,FOEDESU,FOELH,FOELDCP
FOEEW ( PTARE ) = R2ES*EXP (&
  &(R3LES*FOEDELTA(PTARE)+R3IES*(_ONE_-FOEDELTA(PTARE)))*(PTARE-RTT)&
&/ (PTARE-(R4LES*FOEDELTA(PTARE)+R4IES*(_ONE_-FOEDELTA(PTARE)))))

FOEDE ( PTARE ) = &
  &(FOEDELTA(PTARE)*R5ALVCP+(_ONE_-FOEDELTA(PTARE))*R5ALSCP)&
&/ (PTARE-(R4LES*FOEDELTA(PTARE)+R4IES*(_ONE_-FOEDELTA(PTARE))))**2

FOEDESU ( PTARE ) = &
  &(FOEDELTA(PTARE)*R5LES+(_ONE_-FOEDELTA(PTARE))*R5IES)&
&/ (PTARE-(R4LES*FOEDELTA(PTARE)+R4IES*(_ONE_-FOEDELTA(PTARE))))**2

FOELH ( PTARE ) =&
         &FOEDELTA(PTARE)*RLVTT + (_ONE_-FOEDELTA(PTARE))*RLSTT

FOELDCP ( PTARE ) = &
         &FOEDELTA(PTARE)*RALVDCP + (_ONE_-FOEDELTA(PTARE))*RALSDCP
#ifdef DOC

!     *****************************************************************

!           CONSIDERATION OF MIXED PHASES

!     *****************************************************************

!     FOEALFA is calculated to distinguish the three cases:

!                       FOEALFA=1            water phase
!                       FOEALFA=0            ice phase
!                       0 < FOEALFA < 1      mixed phase

!               INPUT : PTARE = TEMPERATURE
#endif
REAL_B :: FOEALFA
FOEALFA (PTARE) = MIN(_ONE_,((MAX(RTICE,MIN(RTWAT,PTARE))-RTICE)&
 &/(RTWAT-RTICE))**2) 

#ifdef DOC

!     Pressure of water vapour at saturation
!        INPUT : PTARE = TEMPERATURE
#endif
REAL_B :: FOEEWM,FOEDEM,FOELDCPM,FOELHM
FOEEWM ( PTARE ) = R2ES *&
     &(FOEALFA(PTARE)*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))+&
  &(_ONE_-FOEALFA(PTARE))*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES)))

FOEDEM ( PTARE ) = FOEALFA(PTARE)*R5ALVCP*(_ONE_/(PTARE-R4LES)**2)+&
             &(_ONE_-FOEALFA(PTARE))*R5ALSCP*(_ONE_/(PTARE-R4IES)**2)

FOELDCPM ( PTARE ) = FOEALFA(PTARE)*RALVDCP+&
            &(_ONE_-FOEALFA(PTARE))*RALSDCP

FOELHM ( PTARE ) =&
         &FOEALFA(PTARE)*RLVTT+(_ONE_-FOEALFA(PTARE))*RLSTT
#ifdef DOC
!     ------------------------------------------------------------------
!     *****************************************************************

!           CONSIDERATION OF DIFFERENT MIXED PHASE FOR CONV

!     *****************************************************************

!     FOEALFCU is calculated to distinguish the three cases:

!                       FOEALFCU=1            water phase
!                       FOEALFCU=0            ice phase
!                       0 < FOEALFCU < 1      mixed phase

!               INPUT : PTARE = TEMPERATURE
#endif
REAL_B :: FOEALFCU 
FOEALFCU (PTARE) = MIN(_ONE_,((MAX(RTICECU,MIN(RTWAT,PTARE))&
&-RTICECU)/(RTWAT-RTICECU))**2) 

#ifdef DOC

!     Pressure of water vapour at saturation
!        INPUT : PTARE = TEMPERATURE
#endif
REAL_B :: FOEEWMCU,FOEDEMCU,FOELDCPMCU,FOELHMCU
FOEEWMCU ( PTARE ) = R2ES *&
     &(FOEALFCU(PTARE)*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))+&
  &(_ONE_-FOEALFCU(PTARE))*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES)))

FOEDEMCU ( PTARE )=FOEALFCU(PTARE)*R5ALVCP*(_ONE_/(PTARE-R4LES)**2)+&
             &(_ONE_-FOEALFCU(PTARE))*R5ALSCP*(_ONE_/(PTARE-R4IES)**2)

FOELDCPMCU ( PTARE ) = FOEALFCU(PTARE)*RALVDCP+&
            &(_ONE_-FOEALFCU(PTARE))*RALSDCP

FOELHMCU ( PTARE ) =&
         &FOEALFCU(PTARE)*RLVTT+(_ONE_-FOEALFCU(PTARE))*RLSTT
!     ------------------------------------------------------------------
#ifdef DOC

!     Pressure of water vapour at saturation
!     This one is for the WMO definition of saturation, i.e. always
!     with respect to water.
#endif
REAL_B :: FOEEWMO
FOEEWMO( PTARE ) = R2ES*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))
